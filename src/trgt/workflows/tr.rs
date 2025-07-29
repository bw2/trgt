use super::{Allele, Genotype, LocusResult};
use crate::hmm::{
    build_hmm, calc_purity, collapse_labels, count_motifs, remove_imperfect_motifs,
    replace_invalid_bases, Annotation, Hmm,
};
use crate::trgt::{
    genotype::{find_tr_spans, genotype_cluster, genotype_flank, genotype_size, Gt},
    locus::Locus,
    reads::HiFiRead,
};
use crate::utils::{Genotyper, Ploidy, Result};
use itertools::{izip, Itertools};
use std::{mem, vec};

pub struct Params {
    pub min_flank_id_frac: f64,
    pub min_read_qual: f32,
    pub search_flank_len: usize,
    pub max_depth: usize,
}

pub fn analyze(mut locus: Locus, params: &Params) -> (Locus, Result<LocusResult>) {
    if locus.ploidy == Ploidy::Zero {
        return (locus, Ok(LocusResult::empty()));
    }

    let reads = mem::take(&mut locus.reads);
    let (reads, spans) = get_spanning_reads(&locus, params, reads);

    const MIN_RQ_FOR_PURITY: f32 = 0.9;
    let (reads, spans) = if params.min_read_qual < MIN_RQ_FOR_PURITY {
        let ret = filter_impure_trs(&locus, &reads, &spans, MIN_RQ_FOR_PURITY);
        if ret.0.len() < reads.len() {
            log::warn!(
                "{}: Filtered out {} impure reads",
                locus.id,
                reads.len() - ret.0.len()
            );
        }
        ret
    } else {
        (reads, spans)
    };

    if reads.is_empty() {
        return (locus, Ok(LocusResult::empty()));
    }

    let trs = reads
        .iter()
        .zip(spans.iter())
        .map(|(r, s)| std::str::from_utf8(&r.bases[s.0..s.1]).unwrap())
        .collect_vec();

    let (mut gt, mut allele_seqs, mut classification) = match locus.genotyper {
        Genotyper::Size => genotype_size::genotype(locus.ploidy, &trs),
        Genotyper::Cluster => genotype_cluster::genotype(locus.ploidy, &trs),
    };

    // Attempt flank re-genotyping only if alleles have similar length
    if gt.len() == 2 && gt[0].size.abs_diff(gt[1].size) <= 10 {
        let snp_result = genotype_flank::genotype(&reads, &trs, &locus.region);
        if let Some((snp_gt, snp_alleles, snp_assignment)) = snp_result {
            (gt, allele_seqs, classification) = (snp_gt, snp_alleles, snp_assignment);
        }
    }

    let annotations = label_with_hmm(&locus, &allele_seqs);

    let spanning_by_hap = [
        classification.iter().filter(|&x| *x == 0).count(),
        classification.iter().filter(|&x| *x == 1).count(),
    ];
    let meth_by_hap = get_meth(&gt, &reads, &spans);
    let mut genotype = Genotype::new();
    for allele_index in 0..gt.len() {
        genotype.push(Allele {
            seq: allele_seqs[allele_index].clone(),
            annotation: annotations[allele_index].clone(),
            ci: gt[allele_index].ci,
            num_spanning: spanning_by_hap[allele_index],
            meth: meth_by_hap[allele_index],
        });
    }

    // Put reference allele first
    if genotype.len() != 1 && genotype[0].seq != locus.tr && genotype[1].seq == locus.tr {
        genotype.swap(0, 1);
        for c in classification.iter_mut() {
            *c = 1 - *c;
        }
    }

    (
        locus,
        Ok(LocusResult {
            genotype,
            reads,
            tr_spans: spans,
            classification,
        }),
    )
}

fn get_spanning_reads(
    locus: &Locus,
    params: &Params,
    reads: Vec<HiFiRead>,
) -> (Vec<HiFiRead>, Vec<(usize, usize)>) {
    let tr_spans = find_tr_spans(
        locus.left_flank.as_bytes(),
        locus.right_flank.as_bytes(),
        &reads,
        params,
    );

    let reads_and_spans = reads
        .into_iter()
        .zip(tr_spans)
        .filter(|(_r, s)| s.is_some())
        .map(|(r, s)| (r, s.unwrap()))
        .collect_vec();
    log::debug!(
        "{}: Found {} spanning reads",
        locus.id,
        reads_and_spans.len()
    );

    if reads_and_spans.is_empty() {
        return (Vec::new(), Vec::new());
    }

    // Check flanks
    let mut reads_and_spans = reads_and_spans
        .into_iter()
        .filter(|(r, s)| {
            s.0 >= params.search_flank_len && r.bases.len() - s.1 >= params.search_flank_len
        })
        .collect_vec();
    log::debug!(
        "{}: {} spanning reads had sufficiently long flanks",
        locus.id,
        reads_and_spans.len()
    );

    if reads_and_spans.is_empty() {
        return (Vec::new(), Vec::new());
    }

    // Sort and (possibly) downsample
    reads_and_spans.sort_by(|(ra, sa), (rb, sb)| {
        (sa.1 - sa.0)
            .cmp(&(sb.1 - sb.0))
            .then_with(|| ra.id.cmp(&rb.id))
    });
    if reads_and_spans.len() > params.max_depth {
        uniform_downsample(&mut reads_and_spans, params.max_depth);
        log::debug!(
            "{}: downsampled to {} reads",
            locus.id,
            reads_and_spans.len()
        );
    }

    let (reads, spans) = reads_and_spans.into_iter().unzip();

    (reads, spans)
}

fn uniform_downsample(reads_and_spans: &mut Vec<(HiFiRead, (usize, usize))>, output_length: usize) {
    let num_reads = reads_and_spans.len() as f64;
    let mut fast: f64 = 0.0;
    let step = num_reads / (output_length as f64);
    for i in 0..output_length {
        let ind = fast.floor() as usize;
        if ind != i {
            reads_and_spans.swap(i, ind);
        }
        fast += step;
    }
    reads_and_spans.truncate(output_length);
}

fn get_meth(gt: &Gt, reads: &[HiFiRead], spans: &[(usize, usize)]) -> Vec<Option<f64>> {
    let mut meths_1 = Vec::new();
    let mut meths_2 = Vec::new();

    for (read, span) in reads.iter().zip(spans.iter()) {
        if let Some(level) = read.meth.as_ref().and_then(|_| get_tr_meth(read, span)) {
            match assign_read(gt, span.1 - span.0) {
                Assignment::First => meths_1.push(level),
                Assignment::Second => meths_2.push(level),
                Assignment::Both => {
                    meths_1.push(level);
                    meths_2.push(level);
                }
                Assignment::None => {}
            };
        }
    }

    let calculate_average = |meths: &Vec<f64>| {
        meths
            .first()
            .map(|_| meths.iter().sum::<f64>() / meths.len() as f64)
    };

    let meth_1 = calculate_average(&meths_1);
    let meth_2 = calculate_average(&meths_2);

    if gt.len() == 2 {
        vec![meth_1, meth_2]
    } else {
        vec![meth_1]
    }
}

enum Assignment {
    First,
    Second,
    Both,
    None,
}

fn assign_read(gt: &Gt, tr_len: usize) -> Assignment {
    if gt.len() == 1 {
        return Assignment::First;
    }

    let hap1_len = gt[0].size;
    let hap2_len = gt[1].size;

    let spans_1 = gt[0].ci.0 <= tr_len && tr_len <= gt[0].ci.1;
    let spans_2 = gt[1].ci.0 <= tr_len && tr_len <= gt[1].ci.1;

    let dist_1 = tr_len.abs_diff(gt[0].size);
    let dist_2 = tr_len.abs_diff(gt[1].size);

    if dist_1 < dist_2 && spans_1 {
        return Assignment::First;
    }

    if dist_2 < dist_1 && spans_2 {
        return Assignment::Second;
    }

    if hap1_len == hap2_len && spans_1 {
        return Assignment::Both;
    }

    Assignment::None
}

fn get_tr_meth(read: &HiFiRead, span: &(usize, usize)) -> Option<f64> {
    if read.meth.is_none() || read.meth.as_ref().unwrap().is_empty() {
        return None;
    }

    let meth = read.meth.as_ref().unwrap();

    let mut total_meth = 0.0;
    let mut cpg_count = 0;
    let mut cpg_index = 0;
    for pos in 0..read.bases.len() - 1 {
        let dinuc = &read.bases[pos..pos + 2];

        if dinuc == b"CG" {
            if span.0 <= pos && pos < span.1 {
                cpg_count += 1;
                total_meth += match meth.get(cpg_index) {
                    Some(value) => *value as f64 / 255.0,
                    None => {
                        log::error!("Read {} has malformed methylation profile", read.id);
                        std::process::exit(1);
                    }
                };
            }

            cpg_index += 1;
        }
    }

    if cpg_count != 0 {
        let mean_meth = total_meth / cpg_count as f64;
        Some(mean_meth)
    } else {
        None
    }
}

fn filter_impure_trs(
    locus: &Locus,
    reads: &[HiFiRead],
    spans: &[(usize, usize)],
    rq_cutoff: f32,
) -> (Vec<HiFiRead>, Vec<(usize, usize)>) {
    let max_filter = std::cmp::max(1_usize, (0.1 * (reads.len() as f64)).round() as usize);
    let mut num_filtered = 0;
    let mut hmm = Hmm::new(0);
    let mut motifs = Vec::new();
    const PURITY_CUTOFF: f64 = 0.9;
    let mut ret = izip!(reads, spans)
        .map(|(read, span)| {
            if let Some(rq) = read.read_qual {
                if rq >= rq_cutoff {
                    return (read, span, 1.0);
                }
            }
            // since HMM building is costly, we will do a "lazy build", i.e., only
            // run the constructor if we find a read with low rq
            if hmm.num_states == 0 {
                motifs = locus
                    .motifs
                    .iter()
                    .map(|m| replace_invalid_bases(m, &['A', 'T', 'C', 'G', 'N']))
                    .map(|m| m.as_bytes().to_vec())
                    .collect_vec();

                hmm = build_hmm(&motifs);
            }
            let seq = std::str::from_utf8(&read.bases[span.0..span.1]).unwrap();
            let seq = replace_invalid_bases(seq, &['A', 'T', 'C', 'G']);
            let labels = hmm.label(&seq);
            let purity = calc_purity(seq.as_bytes(), &hmm, &motifs, &labels);
            (read, span, purity)
        })
        .collect_vec();

    ret.sort_by(|(_, _, pa), (_, _, pb)| f64::total_cmp(pa, pb));
    let ret = ret
        .iter()
        .filter(|(_, _, purity)| {
            if *purity >= PURITY_CUTOFF || num_filtered >= max_filter {
                return true;
            }
            num_filtered += 1;
            false
        })
        .collect_vec();
    let reads = ret.iter().map(|r| (r.0).clone()).collect();
    let spans = ret.iter().map(|r| *(r.1)).collect();
    (reads, spans)
}

fn label_with_hmm(locus: &Locus, seqs: &Vec<String>) -> Vec<Annotation> {
    let motifs = locus
        .motifs
        .iter()
        .map(|m| replace_invalid_bases(m, &['A', 'T', 'C', 'G', 'N']))
        .map(|m| m.as_bytes().to_vec())
        .collect_vec();
    let hmm = build_hmm(&motifs);

    let mut annotations = Vec::new();
    for seq in seqs {
        let seq = replace_invalid_bases(seq, &['A', 'T', 'C', 'G']);
        let labels = hmm.label(&seq);
        let purity = calc_purity(seq.as_bytes(), &hmm, &motifs, &labels);
        let labels = remove_imperfect_motifs(&hmm, &motifs, &labels, seq.as_bytes(), 6);
        let labels = hmm.label_motifs(&labels);
        // Remove labels corresponding to the skip state
        let labels = labels
            .into_iter()
            .filter(|rec| rec.motif_index < motifs.len())
            .collect_vec();
        let motif_counts = count_motifs(&locus.motifs, &labels);
        let labels = collapse_labels(labels);
        // TODO: Consider using empty labels instead of Nones
        let labels = if !labels.is_empty() {
            Some(labels)
        } else {
            None
        };

        annotations.push(Annotation {
            labels,
            motif_counts,
            purity,
        });
    }

    annotations
}
