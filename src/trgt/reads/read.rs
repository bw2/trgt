//! Module for representing and building read information from alignment records.
//!

use super::cigar::Cigar;
use crate::trgt::reads::snp;
use rust_htslib::bam::{self, ext::BamRecordExtensions, record::Aux};
use std::str;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AlleleAssign {
    A0,
    A1,
    Both,
    None,
}

pub type Span = (usize, usize);

#[derive(Debug)]
pub struct SpanningRead {
    pub read: LocusRead,
    pub span: Span,
    pub assign: AlleleAssign,
}

impl SpanningRead {
    #[inline]
    pub fn tr_slice(&self) -> &[u8] {
        &self.read.bases[self.span.0..self.span.1]
    }
    #[inline]
    pub fn tr_len(&self) -> usize {
        self.span.1 - self.span.0
    }
}

/// Represents a single read from an alignment record after clipping to a locus.
pub type LocusRead = HiFiRead;

// NOTE: at construction this encapsulates the entire read, when it is clipped, it is still called a HiFiRead even though it is not the entire read anymore
/// Represents a single read from an alignment record.
#[derive(PartialEq)]
pub struct HiFiRead {
    /// Unique identifier for the read.
    pub id: String,
    /// Flag indicating if the read is from the reverse strand.
    pub is_reverse: bool,
    /// Vector of bases (nucleotides) in the read.
    pub bases: Vec<u8>,
    /// Vector of quality scores for the bases.
    pub quals: Vec<u8>,
    /// Optional vector of methylation calls.
    pub meth: Option<Vec<u8>>,
    /// Optional overall quality score for the read.
    pub read_qual: Option<f32>,
    /// Optional vector of positions of mismatches
    pub mismatch_positions: Option<Vec<u32>>,
    /// Optional CIGAR string representing the alignment.
    pub cigar: Option<Cigar>,
    /// Optional haplotype tag.
    pub hp_tag: Option<u8>,
    /// Optional phase block id tag.
    pub ps_tag: Option<u32>,
    /// Mapping quality score.
    pub mapq: u8,
    /// Reference start position.
    pub ref_start: i64,
    /// Reference end position.
    pub ref_end: i64,
}

impl std::fmt::Debug for HiFiRead {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let meth = self
            .meth
            .as_deref()
            .map(|v| v.iter().map(u8::to_string).collect::<Vec<_>>().join(","))
            .unwrap_or_else(|| "NA".into());

        f.debug_struct("Read")
            .field("id", &self.id)
            .field("bases", &String::from_utf8_lossy(&self.bases))
            .field("meth", &meth)
            .field("cigar", &self.cigar)
            .finish()
    }
}

fn get_meth(rec: &bam::Record, bases: &[u8]) -> Option<Vec<u8>> {
    let reverse = rec.is_reverse();
    let shift = usize::from(reverse); // 1 on reverse, 0 otherwise

    let cpg_indices: Vec<usize> = bases
        .windows(2)
        .enumerate()
        .filter_map(|(i, w)| (w == b"CG").then_some(i + shift))
        .collect();

    let num_cpgs = cpg_indices.len();
    let mut ans = vec![0u8; cpg_indices.len()];
    let mut ind = 0;
    if let Ok(mods) = rec.basemods_iter() {
        for (pos, m) in mods.flatten() {
            if m.canonical_base as u8 == b'C' {
                let pos = pos as usize;
                while ind < num_cpgs && cpg_indices[ind] < pos {
                    ind += 1;
                }
                if ind < num_cpgs && pos == cpg_indices[ind] {
                    ans[ind] = m.qual as u8;
                    ind += 1;
                }
            }
        }
        if ind == 0 {
            // empty MM/ML
            return None;
        }
        if reverse {
            ans.reverse();
        }
        return Some(ans);
    }
    None
}

impl HiFiRead {
    /// Creates a `HiFiRead` from an HTSlib record.
    ///
    /// # Arguments
    /// * `rec` - A BAM record from HTSlib.
    ///
    /// # Returns
    /// Returns a `HiFiRead` populated with data extracted from the BAM record.
    pub fn from_hts_rec(rec: &bam::Record) -> HiFiRead {
        let id = str::from_utf8(rec.qname()).unwrap().to_owned();
        let is_reverse = rec.is_reverse();
        let bases = rec.seq().as_bytes();
        let quals = rec.qual().to_vec();

        let meth = get_meth(rec, &bases);

        let mapq = rec.mapq();
        let hp_tag = get_hp_tag(rec);
        let ps_tag = get_ps_tag(rec);
        let read_qual = get_rq_tag(rec);

        let cigar = if !rec.is_unmapped() {
            Some(Cigar {
                ref_pos: rec.reference_start(),
                ops: rec.cigar().take().to_vec(),
            })
        } else {
            None
        };

        let ref_start = rec.reference_start();
        let ref_end = rec.reference_end();

        let mismatch_positions = cigar.as_ref().map(snp::extract_mismatch_positions);

        HiFiRead {
            id,
            is_reverse,
            bases,
            quals,
            meth,
            read_qual,
            mismatch_positions,
            cigar,
            hp_tag,
            ps_tag,
            mapq,
            ref_start,
            ref_end,
        }
    }
}

/// Retrieves the RQ (read quality) tag from a BAM record.
///
/// # Arguments
/// * `rec` - A reference to the BAM record.
///
/// # Returns
/// Returns an `Option<f32>` which is `Some` if the RQ tag is present and can be parsed as a float, otherwise `None`.
pub fn get_rq_tag(rec: &bam::Record) -> Option<f32> {
    match rec.aux(b"rq") {
        Ok(Aux::Float(value)) => Some(value),
        _ => None,
    }
}

/// Retrieves the HP (haplotype) tag from a BAM record.
///
/// # Arguments
/// * `rec` - A reference to the BAM record.
///
/// # Returns
/// Returns an `Option<u8>` which is `Some` if the HP tag is present and can be parsed as a byte (1 or 2), otherwise `None`.
fn get_hp_tag(rec: &bam::Record) -> Option<u8> {
    match rec.aux(b"HP") {
        Ok(Aux::U8(value)) if value == 1 || value == 2 => Some(value),
        _ => None,
    }
}

/// Retrieves the PS (phase set) tag from a BAM record.
///
/// # Arguments
/// * `rec` - A reference to the BAM record.
///
/// # Returns
/// Returns an `Option<u32>` which is `Some` if the PS tag is present and can be parsed as a u32, otherwise `None`.
fn get_ps_tag(rec: &bam::Record) -> Option<u32> {
    let aux = rec.aux(b"PS").ok()?;
    let val = match aux {
        Aux::I8(v) => v as i64,
        Aux::I16(v) => v as i64,
        Aux::I32(v) => v as i64,
        Aux::U8(v) => v as i64,
        Aux::U16(v) => v as i64,
        Aux::U32(v) => v as i64,
        _ => return None,
    };
    u32::try_from(val).ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::{record::Aux, Record};

    fn create_record(bases: &[u8], mm: &str, ml: &[u8], reverse: bool) -> Record {
        let mut rec = Record::new();
        let qual = vec![40; bases.len()];
        if reverse {
            rec.set_flags(0x10);
        }
        rec.set(b"test", None, bases, &qual);
        rec.push_aux(b"MM", Aux::String(mm)).unwrap();
        rec.push_aux(b"ML", Aux::ArrayU8(ml.into())).unwrap();
        rec
    }

    #[test]
    fn test_basemods_error() {
        let bases = b"ACGTCG";
        let rec = create_record(bases, "no", &[42], false);
        let res = get_meth(&rec, bases);
        assert!(res.is_none());
    }

    #[test]
    fn test_hp_tag_one_and_two_valid() {
        let mut rec1 = create_record(b"N", "C+m,1,0;", &[42], false);
        rec1.push_aux(b"HP", Aux::U8(1)).unwrap();
        assert_eq!(get_hp_tag(&rec1), Some(1));
        let mut rec2 = create_record(b"N", "C+m,1,0;", &[42], false);
        rec2.push_aux(b"HP", Aux::U8(2)).unwrap();
        assert_eq!(get_hp_tag(&rec2), Some(2));
    }

    #[test]
    fn test_hp_tag_non_1_or_2_invalid() {
        let mut rec = create_record(b"N", "C+m,1,0;", &[42], false);
        rec.push_aux(b"HP", Aux::U8(3)).unwrap();
        assert_eq!(get_hp_tag(&rec), None);
    }

    #[test]
    fn test_matching_modifications() {
        let bases = b"AGTCTAGACTCCGTAATTACTCGCCTAG";
        let mm = "C+m,3,1;";
        let ml = [249, 4];
        let rec = create_record(bases, mm, &ml, false);
        let res = get_meth(&rec, bases);
        assert_eq!(res, Some(vec![249, 4]));
    }

    #[test]
    fn test_reverse_matching_modifications() {
        let a = true;
        let b = false;
        let shift_a = usize::from(a);
        let shift_b = usize::from(b);
        assert_eq!(shift_a, 1);
        assert_eq!(shift_b, 0);
    }
}
