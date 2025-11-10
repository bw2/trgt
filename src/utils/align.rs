use crate::{
    commands::genotype::THREAD_WFA_CONSENSUS,
    wfaligner::{CigarOp, WFAligner},
};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

#[derive(Debug, Clone, Copy)]
pub struct TrgtScoring {
    pub mism_scr: i32,
    pub gapo_scr: i32,
    pub gape_scr: i32,
}

pub fn align(backbone: &[u8], seqs: &[&[u8]]) -> Vec<Vec<CigarOp>> {
    let alignments: Vec<Vec<CigarOp>> = seqs
        .par_iter()
        .map(|seq| {
            THREAD_WFA_CONSENSUS.with(|aligner_cell| {
                let mut aligner = aligner_cell.borrow_mut();
                let _status = aligner.align_end_to_end(backbone, seq);
                WFAligner::decode_sam_cigar(&aligner.get_sam_cigar(true))
            })
        })
        .collect();
    alignments
}
