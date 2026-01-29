use super::align::AlleleAlign;
use super::align_consensus::align_consensus;
use super::align_reads::align_reads;
use crate::utils::{locus::InputLocus, read::Read};

pub fn get_allele_align(locus: &InputLocus, consensus: &[u8], reads: &[&Read]) -> AlleleAlign {
    let consensus_align = align_consensus(locus, consensus);
    let reads = align_reads(consensus, &consensus_align, reads);
    AlleleAlign {
        seq: consensus_align,
        reads,
    }
}
