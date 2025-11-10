use crate::commands::genotype::THREAD_WFA_FLANK;
use crate::trgt::{
    reads::{LocusRead, Span},
    workflows::Params,
};
use rayon::prelude::*;

fn find_spans(piece: &[u8], seqs: &[&[u8]], min_n_matches: f64) -> Vec<Option<Span>> {
    seqs.par_iter()
        .map(|s| {
            s.windows(piece.len())
                .position(|window| window == piece)
                .map(|start| (start, start + piece.len()))
                .or_else(|| {
                    THREAD_WFA_FLANK.with(|aligner_cell| {
                        let mut aligner = aligner_cell.borrow_mut();
                        let _status =
                            aligner.align_ends_free(piece, 0, 0, s, s.len() as i32, s.len() as i32);
                        let n_matches = aligner.count_matches() as f64;
                        if n_matches >= min_n_matches {
                            let ((_pattern_start, _pattern_end), (text_start, text_end)) =
                                aligner.get_alignment_span();
                            Some((text_start, text_end))
                        } else {
                            None
                        }
                    })
                })
        })
        .collect()
}

pub fn find_tr_spans(
    lf: &[u8],
    rf: &[u8],
    reads: &[LocusRead],
    params: &Params,
) -> Vec<Option<Span>> {
    let lf_piece = &lf[lf.len() - params.search_flank_len..];
    let rf_piece = &rf[..params.search_flank_len];

    let seqs = reads
        .iter()
        .map(|r| r.bases.as_slice())
        .collect::<Vec<&[u8]>>();

    let min_n_matches = (params.search_flank_len as f64) * params.min_flank_id_frac;
    let (lf_spans, rf_spans) = rayon::join(
        || find_spans(lf_piece, &seqs, min_n_matches),
        || find_spans(rf_piece, &seqs, min_n_matches),
    );

    lf_spans
        .into_iter()
        .zip(rf_spans)
        .map(|(lf_span, rf_span)| match (lf_span, rf_span) {
            (None, None) => None,      // No left or right span
            (Some(_lf), None) => None, // Left flanking only
            (None, Some(_rf)) => None, // Right flanking only
            (Some(lf), Some(rf)) => {
                if lf.1 <= rf.0 {
                    Some((lf.1, rf.0))
                } else {
                    None // Discordant flanks (overlap or wrong order)
                }
            }
        })
        .collect()
}
