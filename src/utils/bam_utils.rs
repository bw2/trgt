use crate::utils::{open_bam_reader, InputSource, Result};
use rust_htslib::bam::{self, Read};
use std::collections::HashSet;

pub fn get_bam_header(src: &InputSource) -> Result<bam::Header> {
    let bam = open_bam_reader(src, 1)?;
    Ok(bam::Header::from_template(bam.header()))
}

pub fn is_bam_mapped(bam_header: &bam::Header) -> bool {
    // input is already sorted because it fails an index.
    // If it is mapped, the index needs the SQ tags to fetch data.
    for line in String::from_utf8(bam_header.to_bytes()).unwrap().lines() {
        if line.starts_with("@SQ") {
            return true;
        }
    }
    false
}

pub fn get_sample_name(reads_path: &InputSource, bam_header: &bam::Header) -> Result<String> {
    let header_hashmap = bam_header.to_hashmap();
    let mut sample_names = HashSet::new();

    if let Some(rg_fields) = header_hashmap.get("RG") {
        for rg_field in rg_fields {
            if let Some(sample_name) = rg_field.get("SM") {
                sample_names.insert(sample_name.to_owned());
            }
        }
    }

    match sample_names.len() {
        1 => return Ok(sample_names.into_iter().next().unwrap()),
        0 => log::warn!("No sample names found"),
        _ => log::warn!("Multiple sample names found"),
    };

    let sample = reads_path
        .file_stem()
        .ok_or_else(|| "Invalid reads input; cannot infer sample name from path/URL".to_string())?;

    Ok(sample)
}
