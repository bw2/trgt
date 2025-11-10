use super::{cigar::Cigar, HiFiRead};
use rust_htslib::bam::record::CigarString;

pub fn make_read(bases: &str, meths: Vec<u8>, cigar: Cigar) -> HiFiRead {
    HiFiRead {
        id: "test_read".to_string(),
        is_reverse: false,
        bases: bases.as_bytes().to_vec(),
        quals: "(".repeat(bases.len()).as_bytes().to_vec(),
        meth: Some(meths),
        read_qual: None,
        mismatch_positions: None,
        cigar: Some(cigar),
        hp_tag: None,
        mapq: 60,
        ref_start: 0,
        ref_end: 0,
    }
}

pub fn make_cigar(ref_pos: i64, encoding: &str) -> Cigar {
    let ops = CigarString::try_from(encoding).unwrap().to_vec();
    Cigar { ref_pos, ops }
}
