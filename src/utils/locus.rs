use crate::{
    trgt::locus::{check_region_bounds, decode_fields, get_field, get_tr_and_flanks},
    utils::{GenomicRegion, Result},
};
use rust_htslib::faidx;
use std::collections::HashMap;

#[derive(Debug, PartialEq, Clone)]
pub enum RegionLabel {
    Flank(usize, usize),       // Coordinates
    Tr(usize, usize, Vec<u8>), // Coordinates, Motif
    Seq(usize, usize),
    Other(usize, usize),
}

#[derive(Debug)]
pub struct Allele {
    pub seq: Vec<u8>,
    pub region_labels: Vec<RegionLabel>,
    pub flank_labels: Vec<RegionLabel>,
    //pub base_labels: Vec<BaseLabel>,
}

#[derive(Debug)]
pub struct InputLocus {
    pub id: String,
    pub struc: String,
    pub motifs: Vec<Vec<u8>>,
    pub left_flank: Vec<u8>,
    pub right_flank: Vec<u8>,
    pub region: GenomicRegion,
}

impl InputLocus {
    #[cfg(test)]
    pub fn builder() -> InputLocusBuilder {
        InputLocusBuilder::default()
    }

    pub fn new(
        genome_reader: &faidx::Reader,
        chrom_lookup: &HashMap<String, u32>,
        line: &str,
        flank_len: usize,
    ) -> Result<Self> {
        const EXPECTED_FIELD_COUNT: usize = 4;
        let split_line: Vec<&str> = line.split_whitespace().collect();
        if split_line.len() != EXPECTED_FIELD_COUNT {
            return Err(format!(
                "Expected {} fields in the format 'chrom start end info', found {}: {}",
                EXPECTED_FIELD_COUNT,
                split_line.len(),
                line
            ));
        }

        let (chrom, start, end, info_fields) = match &split_line[..] {
            [chrom, start, end, info_fields] => (*chrom, *start, *end, *info_fields),
            _ => unreachable!(),
        };

        let region = GenomicRegion::from_string(&format!("{}:{}-{}", chrom, start, end))?;
        check_region_bounds(&region, flank_len, chrom_lookup)?;

        let fields = decode_fields(info_fields)?;
        let id = get_field(&fields, "ID")?;

        let motifs = get_field(&fields, "MOTIFS")?
            .split(',')
            .map(|s| s.as_bytes().to_vec())
            .collect();
        let struc = get_field(&fields, "STRUC")?;

        let (left_flank, _, right_flank) = get_tr_and_flanks(genome_reader, &region, flank_len)?;

        Ok(InputLocus {
            id,
            struc,
            motifs,
            left_flank,
            right_flank,
            region,
        })
    }
}

#[cfg(test)]
pub struct InputLocusBuilder {
    id: String,
    struc: String,
    motifs: Vec<Vec<u8>>,
    left_flank: Vec<u8>,
    right_flank: Vec<u8>,
    region: GenomicRegion,
}

#[cfg(test)]
impl InputLocusBuilder {
    pub fn new() -> Self {
        Self {
            id: "tr1".to_string(),
            struc: "".to_string(),
            motifs: vec![],
            left_flank: b"L".to_vec(),
            right_flank: b"R".to_vec(),
            region: GenomicRegion::new("chr1", 0, 10).unwrap(),
        }
    }

    pub fn id(mut self, id: &str) -> Self {
        self.id = id.to_string();
        self
    }

    pub fn struc(mut self, struc: &str) -> Self {
        self.struc = struc.to_string();
        self
    }

    pub fn motifs(mut self, motifs: Vec<Vec<u8>>) -> Self {
        self.motifs = motifs;
        self
    }

    pub fn left_flank(mut self, flank: Vec<u8>) -> Self {
        self.left_flank = flank;
        self
    }

    pub fn right_flank(mut self, flank: Vec<u8>) -> Self {
        self.right_flank = flank;
        self
    }

    pub fn region(mut self, region: GenomicRegion) -> Self {
        self.region = region;
        self
    }

    pub fn build(self) -> InputLocus {
        InputLocus {
            id: self.id,
            struc: self.struc,
            motifs: self.motifs,
            left_flank: self.left_flank,
            right_flank: self.right_flank,
            region: self.region,
        }
    }
}

#[cfg(test)]
impl Default for InputLocusBuilder {
    fn default() -> Self {
        Self::new()
    }
}
