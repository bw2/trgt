use crate::cli;
use crate::cli::DeepdiveArgs;
use crate::trgt::locus::create_chrom_lookup;
use crate::utils::locus::Locus;
use crate::utils::{input, open_catalog_reader, open_genome_reader, Result};
use crate::wfaligner::{AlignmentScope, MemoryModel, WFAligner, WfaOp};
use itertools::Itertools;
use rust_htslib::bam::header::{Header, HeaderRecord};
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use rust_htslib::bam::{self, Record, Writer};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};

/// Convert individual WFA ops to rust htslib
fn op_to_cigar(op: WfaOp, len: u32) -> Cigar {
    match op {
        WfaOp::Match => Cigar::Equal(len),
        WfaOp::Subst => Cigar::Diff(len),
        WfaOp::Ins => Cigar::Ins(len),
        WfaOp::Del => Cigar::Del(len),
    }
}

/// Convert the WFA Cigar strings to the ones from rust_htslib
fn wfa_cigar_to_htslib(ops: &[WfaOp]) -> Vec<Cigar> {
    ops.iter()
        .chunk_by(|&&op| op)
        .into_iter()
        .map(|(op, group)| op_to_cigar(op, group.count() as u32))
        .collect()
}

/// Given a DNA sequence, make the MM tag for CpGs
fn get_mm_tag(seq: &[u8]) -> String {
    let mut deltas = Vec::with_capacity(seq.len() / 2);
    let mut c_since_last_cpg = 0;
    for i in seq
        .iter()
        .enumerate()
        .filter_map(|(i, &b)| if b == b'C' { Some(i) } else { None })
    {
        if i + 1 < seq.len() && seq[i + 1] == b'G' {
            deltas.push(c_since_last_cpg.to_string());
            c_since_last_cpg = 0;
        } else {
            c_since_last_cpg += 1;
        }
    }
    format!("C+m?,{};", deltas.join(","))
}

/// Get the relevant info for a locus to avoid re-opening the caller
fn find_locus_line(catalog_reader: impl BufRead, tr_id: &str) -> Result<String> {
    let query = format!("ID={};", tr_id);
    for result_line in catalog_reader.lines() {
        let line = result_line.map_err(|e| format!("Error reading catalog file: {}", e))?;
        if line.contains(&query) {
            return Ok(line);
        }
    }
    Err(format!("Unable to find locus {}", tr_id))
}

/// Align .spanning.bam reads to their consensus and produce
/// outputs meant to be used on downstream tools, using the
/// consensus sequences as the reference
pub fn deepdive(args: DeepdiveArgs) -> Result<()> {
    let catalog_reader = open_catalog_reader(&args.repeats_path)?;
    let locus_line = find_locus_line(catalog_reader, &args.tr_id)?;

    let genome_reader = open_genome_reader(&args.genome_path)?;
    let chrom_lookup = create_chrom_lookup(&genome_reader)?;

    let locus_prelim = Locus::new(&genome_reader, &chrom_lookup, &locus_line, 0)
        .map_err(|e| format!("Error parsing locus line: {}", e))?;
    let reads = input::get_reads(&args.reads_path, &locus_prelim, None)?;

    // largest FL value in spanning.bam reads to be used to fetch flanks
    let flank_length: usize = reads
        .iter()
        .map(|r| std::cmp::max(r.left_flank, r.right_flank))
        .max()
        .unwrap();

    let locus = Locus::new(&genome_reader, &chrom_lookup, &locus_line, flank_length)
        .map_err(|e| format!("Error parsing locus line: {}", e))?;

    let allele_seqs = input::get_alleles(&args.bcf_path, &locus)?;

    let mut output_path_fasta = args.output_prefix.to_path_buf().into_os_string();
    output_path_fasta.push(".fasta");

    // Add program command to output BAM headers
    let command_line = env::args().collect::<Vec<String>>().join(" ");
    let mut out_header = Header::new();
    let mut pg_record = HeaderRecord::new(b"PG");
    pg_record.push_tag(b"ID", env!("CARGO_PKG_NAME"));
    pg_record.push_tag(b"PN", env!("CARGO_PKG_NAME"));
    pg_record.push_tag(b"CL", command_line);
    pg_record.push_tag(b"VN", (*cli::FULL_VERSION).to_string());
    out_header.push_record(&pg_record);

    let consensus_file = File::create(output_path_fasta)
        .map_err(|e| format!("Error creating FASTA output: {}", e))?;
    let mut consensus_writer = BufWriter::new(consensus_file);

    let mut output_path_bed = args.output_prefix.to_path_buf().into_os_string();
    output_path_bed.push(".bed");
    let annotation_file =
        File::create(output_path_bed).map_err(|e| format!("Error creating BED output: {}", e))?;
    let mut annotation_writer = BufWriter::new(annotation_file);
    for (i, seq) in allele_seqs.iter().enumerate() {
        let consensus_name = format!("{}_{}", &locus.id, i);
        let consensus_len = seq.len();
        writeln!(consensus_writer, ">{}\n{}", consensus_name, seq)
            .map_err(|e| format!("Error writing entry into consensus FASTA: {}", e))?;
        writeln!(
            annotation_writer,
            "{}\t{}\t{}\tleft_flank",
            consensus_name, 0, flank_length
        )
        .map_err(|e| format!("Error writing left flank region into annotation BED: {}", e))?;

        writeln!(
            annotation_writer,
            "{}\t{}\t{}\ttandem_repeat",
            consensus_name,
            flank_length,
            consensus_len - flank_length
        )
        .map_err(|e| {
            format!(
                "Error writing tandem repeat region into annotation BED: {}",
                e
            )
        })?;
        writeln!(
            annotation_writer,
            "{}\t{}\t{}\tright_flank",
            consensus_name,
            consensus_len - flank_length,
            consensus_len
        )
        .map_err(|e| {
            format!(
                "Error writing right flank region into annotation BED: {}",
                e
            )
        })?;
        out_header.push_record(&HeaderRecord::new(
            format!("SQ\tSN:{}\tLN:{}", consensus_name, consensus_len).as_bytes(),
        ));
    }
    consensus_writer
        .flush()
        .map_err(|e| format!("Error flushing FASTA writer: {}", e))?;
    annotation_writer
        .flush()
        .map_err(|e| format!("Error flushing BED writer: {}", e))?;

    let mut output_path_bam = args.output_prefix.to_path_buf().into_os_string();
    output_path_bam.push(".bam");

    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(2, 5, 1)
        .build();

    let left_flank = locus.left_flank.as_bytes();
    let right_flank = locus.right_flank.as_bytes();

    // NB: allele sequences have flanks as well
    let allele_bytes = allele_seqs
        .iter()
        .map(|x| &x.as_bytes()[left_flank.len()..x.len() - right_flank.len()])
        .collect::<Vec<&[u8]>>();

    let mut bam_records = reads
        .iter()
        .map(|read| {
            let ind = read.allele as usize;
            let seq = read.seq.as_bytes();

            // locus slices
            let align_soft_clip_start = left_flank.len() - read.left_flank;
            let locus_left_flank_right_side = &left_flank[align_soft_clip_start..];
            let locus_right_flank_left_side = &right_flank[..read.right_flank];
            let align_soft_clip_end = right_flank.len() - read.right_flank;

            // read slices
            let read_left_flank = &seq[..read.left_flank];
            let read_tr = &seq[read.left_flank..seq.len() - read.right_flank];
            let read_right_flank = &seq[seq.len() - read.right_flank..];

            // block alignments
            let _status = aligner.align_end_to_end(locus_left_flank_right_side, read_left_flank);
            let align_lf = aligner.get_alignment();
            let mut cigar_wfa = align_lf.operations;

            let _status = aligner.align_end_to_end(allele_bytes[ind], read_tr);
            let align_tr = aligner.get_alignment();
            cigar_wfa.extend(align_tr.operations);

            let _status = aligner.align_end_to_end(locus_right_flank_left_side, read_right_flank);
            let align_rf = aligner.get_alignment();
            cigar_wfa.extend(align_rf.operations);

            let mut cigar = if align_soft_clip_start == 0 {
                Vec::new()
            } else {
                vec![Cigar::SoftClip(align_soft_clip_start as u32)]
            };

            for op in wfa_cigar_to_htslib(&cigar_wfa) {
                cigar.push(op);
            }

            if align_soft_clip_end != 0 {
                cigar.push(Cigar::SoftClip(align_soft_clip_end as u32));
            }

            let mut rec = Record::new();
            let quals = vec![255_u8; seq.len()];
            rec.set(
                read.read_name.as_bytes(),
                Some(&CigarString(cigar)),
                seq,
                &quals,
            );
            rec.set_tid(ind as i32);
            rec.set_mapq(255_u8);

            let ml_tag = read
                .betas
                .iter()
                .map(|x| (255.0 * x.value).round() as u8)
                .collect::<Vec<u8>>();
            let mm_tag = get_mm_tag(seq);
            if !ml_tag.is_empty() {
                rec.push_aux(b"MM", Aux::String(mm_tag.as_str())).unwrap();
                rec.push_aux(b"ML", Aux::ArrayU8((&ml_tag).into())).unwrap();
            }
            rec.push_aux(b"HP", Aux::I32(ind as i32 + 1)).unwrap();
            rec
        })
        .collect::<Vec<Record>>();
    bam_records.sort_by_key(|r| (r.tid(), r.pos()));
    {
        let mut bam_writer = Writer::from_path(&output_path_bam, &out_header, bam::Format::Bam)
            .map_err(|e| format!("Error creating BAM output: {}", e))?;
        for rec in &bam_records {
            bam_writer
                .write(rec)
                .map_err(|e| format!("Error writing BAM entry to output: {}", e))?;
        }
    }
    rust_htslib::bam::index::build(
        &output_path_bam,
        None,
        rust_htslib::bam::index::Type::Bai,
        1,
    )
    .map_err(|e| format!("Failed to index BAM file: {}", e))?;
    Ok(())
}
