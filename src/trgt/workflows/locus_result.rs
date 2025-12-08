use crate::hmm::spans::Annotation;
use crate::trgt::reads::SpanningRead;
use crate::trgt::workflows::tr::PhaseResult;
use arrayvec::ArrayVec;

#[derive(Debug)]
pub struct Allele {
    pub seq: Vec<u8>,
    pub annotation: Annotation,
    pub ci: (usize, usize),
    pub num_spanning: usize,
    pub meth: Option<f64>,
}

impl Allele {
    #[cfg(test)]
    pub fn test_base(seq: &[u8]) -> Self {
        Self {
            seq: seq.to_vec(),
            annotation: Annotation::base(),
            ci: (0, 0),
            num_spanning: 0,
            meth: None,
        }
    }
}

pub type Genotype = ArrayVec<Allele, 2>;

#[derive(Debug)]
pub struct LocusResult {
    pub genotype: Genotype,
    pub spanning_reads: Vec<SpanningRead>,
    pub phase_result: Option<PhaseResult>,
}

impl LocusResult {
    pub fn empty() -> LocusResult {
        LocusResult {
            genotype: Genotype::new(),
            spanning_reads: Vec::new(),
            phase_result: None,
        }
    }
}
