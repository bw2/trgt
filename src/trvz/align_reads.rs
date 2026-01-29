use super::align::{Align, AlignOp, AlignSeg};
use super::params::base_color;
use crate::utils::read::{project_betas, Betas, Read};
use crate::wfaligner::{AlignmentScope, MemoryModel, WFAligner, WfaAlign, WfaOp};
use itertools::Itertools;
use pipeplot::TextLabel;

/// Align reads to the consensus sequence
pub fn align_reads(
    consensus: &[u8],
    consensus_align: &[AlignSeg],
    reads: &[&Read],
) -> Vec<(Align, Betas, Vec<TextLabel>)> {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(2, 5, 1)
        .build();

    let mut ret: Vec<(Align, Betas, Vec<TextLabel>, i32, usize)> = reads
        .iter()
        .map(|read| {
            let _status = aligner.align_end_to_end(consensus, &read.seq);
            let wfa_align = aligner.get_alignment();
            let align = convert(consensus_align, &wfa_align);
            let betas = project_betas(&wfa_align, &read.betas);
            let labels = get_mismatch_labels(&read.seq, &wfa_align);
            (align, betas, labels, wfa_align.score, read.seq.len())
        })
        .collect();
    ret.sort_by_key(|r| (r.4, std::cmp::Reverse(r.3)));
    ret.iter()
        .map(|r| (r.0.clone(), r.1.clone(), r.2.clone()))
        .collect()
}

/// Convert a WFA alignment into an internal alignment
fn convert(consensus_align: &[AlignSeg], wfa_align: &WfaAlign) -> Align {
    assert!(wfa_align.xstart == 0, "WFA alignment xstart should be 0");
    assert!(wfa_align.ystart == 0, "WFA alignment ystart should be 0");

    let mut seg_type_by_ref = Vec::new();
    for align_seg in consensus_align {
        let seg_type = align_seg.seg_type;
        seg_type_by_ref.extend(match align_seg.op {
            AlignOp::Del | AlignOp::Match | AlignOp::Subst => vec![seg_type; align_seg.width],
            AlignOp::Ins => vec![],
        });
    }

    let mut ref_pos = 0;
    let mut ops_and_segs = Vec::with_capacity(wfa_align.operations.len());
    for op in &wfa_align.operations {
        let seg_type = if ref_pos == seg_type_by_ref.len() {
            // Handle trailing insertion
            assert_eq!(*op, WfaOp::Ins);
            seg_type_by_ref[ref_pos - 1]
        } else {
            seg_type_by_ref[ref_pos]
        };
        ops_and_segs.push((*op, seg_type));
        ref_pos += match *op {
            WfaOp::Match | WfaOp::Subst | WfaOp::Del => 1,
            WfaOp::Ins => 0,
        };
    }

    let mut align = Vec::new();
    let mut ref_pos = 0;
    let mut seq_pos = 0;
    for ((wfa_op, seg_type), group) in &ops_and_segs.iter().chunk_by(|rec| *rec) {
        let seg_type = *seg_type;
        let run_len = group.count();

        let align_seg = match *wfa_op {
            WfaOp::Match => AlignSeg {
                width: run_len,
                op: AlignOp::Match,
                seg_type,
            },
            WfaOp::Subst => AlignSeg {
                width: run_len,
                op: AlignOp::Subst,
                seg_type,
            },
            WfaOp::Del => AlignSeg {
                width: run_len,
                op: AlignOp::Del,
                seg_type,
            },
            WfaOp::Ins => AlignSeg {
                width: 0,
                op: AlignOp::Ins,
                seg_type,
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

/// Generate text labels for bases that differ from consensus
fn get_mismatch_labels(read_seq: &[u8], wfa_align: &WfaAlign) -> Vec<TextLabel> {
    let mut labels = Vec::new();
    let mut cons_pos = 0usize;
    let mut read_pos = 0usize;

    for op in &wfa_align.operations {
        match op {
            WfaOp::Match => {
                cons_pos += 1;
                read_pos += 1;
            }
            WfaOp::Subst => {
                let base = read_seq[read_pos];
                labels.push(TextLabel {
                    pos: cons_pos as u32,
                    text: base as char,
                    color: base_color(base).to_string(),
                });
                cons_pos += 1;
                read_pos += 1;
            }
            WfaOp::Del => {
                cons_pos += 1;
            }
            WfaOp::Ins => {
                read_pos += 1;
            }
        }
    }

    labels
}
