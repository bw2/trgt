use super::align::Align;
use super::params::{base_color, get_meth_colors, Color, ColorMap, PlotParams};
use super::scale::get_multi_scale;
use crate::trvz::{
    align::{AlignOp, AlignSeg, SegType},
    align_consensus::align_motifs_with_labels,
};
use crate::utils::{
    locus::InputLocus,
    read::{project_betas, Beta, Betas, Read},
};
use crate::wfaligner::{AlignmentScope, MemoryModel, WFAligner, WfaAlign, WfaOp};
use itertools::Itertools;
use pipeplot::{Band, FontConfig, Legend, Pipe, PipePlot, Seg, Shape, TextLabel};

pub fn plot_waterfall(
    locus: &InputLocus,
    what_to_show: &str,
    reads: &[Read],
    params: &PlotParams,
) -> PipePlot {
    let reads = reads
        .iter()
        .sorted_by(|r1, r2| r1.seq.len().cmp(&r2.seq.len()))
        .collect_vec();

    let longest_read = reads.iter().map(|r| r.seq.len()).max().unwrap();
    let aligned_reads = reads
        .iter()
        .map(|r| align(locus, longest_read, r))
        .collect_vec();

    plot(locus, what_to_show, &aligned_reads, params)
}

fn align(locus: &InputLocus, longest_read: usize, read: &Read) -> (Align, Vec<Beta>, Vec<TextLabel>) {
    let lf_ref = &locus.left_flank;
    let rf_ref = &locus.right_flank;

    let lf_read = &read.seq[..lf_ref.len()];
    let rf_read = &read.seq[read.seq.len() - locus.right_flank.len()..];

    // Align left flank and get mismatch labels
    let lf_wfa_align = get_flank_align(lf_ref, lf_read);
    let mut align = convert(&lf_wfa_align, SegType::LeftFlank);
    let mut labels = get_flank_mismatch_labels(lf_read, &lf_wfa_align, 0);

    // Calculate pipe position after left flank alignment
    let lf_pipe_len: usize = align.iter().map(|seg| seg.width).sum();

    // Align TR region and get mismatch labels
    let tr = &read.seq[locus.left_flank.len()..read.seq.len() - locus.right_flank.len()];
    let (motif_aligns, tr_labels) = align_motifs_with_labels(&locus.motifs, tr, lf_pipe_len);
    align.extend(motif_aligns);
    labels.extend(tr_labels);

    // Calculate pipe position after TR alignment
    let tr_pipe_len: usize = align.iter().map(|seg| seg.width).sum();

    // Add deletion that lines up right flanks
    let deletion_width = longest_read.saturating_sub(read.seq.len());
    // prevents adding a 0-width segment
    if deletion_width > 0 {
        align.push(AlignSeg {
            width: deletion_width,
            op: AlignOp::Del,
            seg_type: SegType::RightFlank,
            insertion_size: 0,
        });
    }

    // Calculate offset for right flank labels (after TR + deletion)
    let rf_offset = tr_pipe_len + deletion_width;

    // Align right flank and get mismatch labels
    let rf_wfa_align = get_flank_align(rf_ref, rf_read);
    align.extend(convert(&rf_wfa_align, SegType::RightFlank));
    labels.extend(get_flank_mismatch_labels(rf_read, &rf_wfa_align, rf_offset));

    let mut proj_betas = Vec::new();

    let lf_betas = &read
        .betas
        .iter()
        .filter(|beta| beta.pos < lf_read.len())
        .cloned()
        .collect_vec();
    proj_betas.extend(project_betas(&lf_wfa_align, lf_betas));

    let tr_betas = &read
        .betas
        .iter()
        .filter(|beta| lf_read.len() <= beta.pos && beta.pos < lf_read.len() + tr.len())
        .map(|beta| Beta {
            pos: beta.pos - lf_read.len(),
            value: beta.value,
        })
        .collect_vec();
    proj_betas.extend(tr_betas.iter().map(|beta| Beta {
        value: beta.value,
        pos: beta.pos + lf_read.len(),
    }));

    let rf_betas = &read
        .betas
        .iter()
        .filter(|beta| lf_read.len() + tr.len() <= beta.pos)
        .map(|beta| Beta {
            pos: beta.pos - lf_read.len() - tr.len(),
            value: beta.value,
        })
        .collect_vec();

    proj_betas.extend(
        project_betas(&rf_wfa_align, rf_betas)
            .iter()
            .map(|beta| Beta {
                pos: beta.pos + lf_read.len() + tr.len() + longest_read - read.seq.len(),
                value: beta.value,
            }),
    );

    assert_eq!(
        read.betas.len(),
        lf_betas.len() + tr_betas.len() + rf_betas.len()
    );

    (align, proj_betas, labels)
}

/// Generate mismatch labels for flank alignment
fn get_flank_mismatch_labels(
    read_seq: &[u8],
    wfa_align: &WfaAlign,
    offset: usize,
) -> Vec<TextLabel> {
    let mut labels = Vec::new();
    let mut read_pos = 0usize;
    let mut pipe_pos = offset;

    for op in &wfa_align.operations {
        match op {
            WfaOp::Match => {
                // Match - no label needed
                read_pos += 1;
                pipe_pos += 1;
            }
            WfaOp::Subst => {
                // Mismatch - emit label
                let base = read_seq[read_pos];
                labels.push(TextLabel {
                    pos: pipe_pos as u32,
                    text: base as char,
                    color: base_color(base).to_string(),
                });
                read_pos += 1;
                pipe_pos += 1;
            }
            WfaOp::Del => {
                // Deletion in read (gap in read) - advances pipe but not read
                pipe_pos += 1;
            }
            WfaOp::Ins => {
                // Insertion in read (extra base) - advances read but not pipe
                read_pos += 1;
            }
        }
    }

    labels
}

fn get_flank_align(ref_seq: &[u8], read_seq: &[u8]) -> WfaAlign {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(2, 5, 1)
        .build();
    let _status = aligner.align_end_to_end(ref_seq, read_seq);
    aligner.get_alignment()
}

/// Convert a rust-wfa alignment into an internal alignments
fn convert(wfa_align: &WfaAlign, seg_type: SegType) -> Align {
    assert!(wfa_align.xstart == 0, "WFA alignment xstart should be 0");
    assert!(wfa_align.ystart == 0, "WFA alignment ystart should be 0");

    let mut align = Vec::new();
    let mut ref_pos = 0;
    let mut seq_pos = 0;

    for (wfa_op, group) in &wfa_align.operations.iter().chunk_by(|op| *op) {
        let run_len = group.count();

        let align_seg = match *wfa_op {
            WfaOp::Match => AlignSeg {
                width: run_len,
                op: AlignOp::Match,
                seg_type,
                insertion_size: 0,
            },
            WfaOp::Subst => AlignSeg {
                width: run_len,
                op: AlignOp::Subst,
                seg_type,
                insertion_size: 0,
            },
            WfaOp::Del => AlignSeg {
                width: run_len,
                op: AlignOp::Del,
                seg_type,
                insertion_size: 0,
            },
            WfaOp::Ins => AlignSeg {
                width: 0,
                op: AlignOp::Ins,
                seg_type,
                insertion_size: run_len,
            },
        };

        align.push(align_seg);

        ref_pos += match *wfa_op {
            WfaOp::Match | WfaOp::Subst | WfaOp::Del => run_len,
            WfaOp::Ins => 0,
        };

        seq_pos += match *wfa_op {
            WfaOp::Match | WfaOp::Subst | WfaOp::Ins => run_len,
            WfaOp::Del => 0,
        };
    }

    assert_eq!(
        wfa_align.ylen, seq_pos,
        "Sequence length mismatch after conversion"
    );
    assert_eq!(
        wfa_align.xlen, ref_pos,
        "Reference length mismatch after conversion"
    );

    align
}

pub fn plot(
    locus: &InputLocus,
    what_to_show: &str,
    reads: &[(Align, Vec<Beta>, Vec<TextLabel>)],
    params: &PlotParams,
) -> PipePlot {
    let xpos = 0;
    let mut ypos = 0;
    let mut pipes = Vec::new();
    let scale_pipes = get_multi_scale(xpos, ypos, params.pipe_height, &reads.last().unwrap().0, locus);
    pipes.extend(scale_pipes);
    ypos += 4;

    for (align, betas, mismatch_labels) in reads {
        let (colors, betas, labels) = if what_to_show == "meth" {
            (get_meth_colors(&locus.motifs), betas.clone(), Vec::new())
        } else {
            (params.colors.clone(), Vec::new(), mismatch_labels.clone())
        };
        let pipe = get_pipe(xpos, ypos, params.pipe_height, align, &betas, &colors, labels);
        pipes.push(pipe);
        ypos += params.pipe_height + params.pipe_pad;
    }

    let mut labels = Vec::new();
    if what_to_show == "motifs" {
        for (index, motif) in locus.motifs.iter().enumerate() {
            let color = params.colors.get(&SegType::Tr(index)).unwrap().to_string();
            labels.push((std::str::from_utf8(motif).unwrap().to_string(), color));
        }
    } else {
        labels = vec![
            ("Methylated".to_string(), Color::Grad(1.0).to_string()),
            ("Unmethylated".to_string(), Color::Grad(0.0).to_string()),
        ];
    }
    ypos += 1;

    let legend = Legend {
        xpos,
        ypos,
        height: 4,
        labels,
    };

    PipePlot {
        pipes,
        legend,
        font: FontConfig::default(),
        grid_lines: Vec::new(),
    }
}

fn get_pipe(
    xpos: u32,
    ypos: u32,
    height: u32,
    align: &Align,
    betas: &Betas,
    colors: &ColorMap,
    labels: Vec<TextLabel>,
) -> Pipe {
    let segs = align
        .iter()
        .map(|align_seg| {
            let shape = match align_seg.op {
                AlignOp::Del => Shape::HLine,
                AlignOp::Ins => Shape::VLine(align_seg.insertion_size as u32),
                AlignOp::Match | AlignOp::Subst => Shape::Rect,
            };
            let color = match align_seg.op {
                AlignOp::Match => colors.get(&align_seg.seg_type).unwrap(),
                AlignOp::Subst => &Color::Gray,
                AlignOp::Del => &Color::LightGray,
                _ => &Color::Black,
            };
            Seg {
                width: align_seg.width as u32,
                color: color.to_string(),
                shape,
                dashed: false,
            }
        })
        .collect();

    let bands = betas
        .iter()
        .map(|beta| Band {
            pos: beta.pos as u32,
            width: 2,
            color: Color::Grad(beta.value).to_string(),
        })
        .collect_vec();

    Pipe {
        xpos,
        ypos,
        height,
        segs,
        bands,
        outline: false,
        labels,
    }
}
