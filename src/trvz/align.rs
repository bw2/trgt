use crate::utils::read::Betas;
use pipeplot::TextLabel;

#[derive(Debug, Clone)]
pub struct AlleleAlign {
    pub seq: Align,
    pub reads: Vec<(Align, Betas, Vec<TextLabel>)>,
}

pub type Align = Vec<AlignSeg>;

#[derive(Debug, Clone)]
pub struct AlignSeg {
    pub width: usize,
    pub op: AlignOp,
    pub seg_type: SegType,
    pub insertion_size: usize,  // For Ins operations, stores the actual insertion size
}

#[derive(Debug, Clone, PartialEq)]
pub enum AlignOp {
    Match,
    Subst,
    Ins,
    Del,
}

#[derive(Debug, Clone, Copy, PartialEq, Hash, Eq)]
pub enum SegType {
    // The integer is the motif index or motifs.len() for unsegmented regions
    Tr(usize),
    LeftFlank,
    RightFlank,
}
