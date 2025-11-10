use crate::hmm::spans::Annotation;
use crate::trgt::reads::SpanningRead;
use arrayvec::ArrayVec;

#[derive(Debug)]
pub struct Allele {
    pub seq: Vec<u8>,
    pub annotation: Annotation,
    pub ci: (usize, usize),
    pub num_spanning: usize,
    pub meth: Option<f64>,
}

pub type Genotype = ArrayVec<Allele, 2>;

#[derive(Debug)]
pub struct LocusResult {
    pub genotype: Genotype,
    pub spanning_reads: Vec<SpanningRead>,
}

impl LocusResult {
    pub fn empty() -> LocusResult {
        LocusResult {
            genotype: Genotype::new(),
            spanning_reads: Vec::new(),
        }
    }
}
