use crate::utils::{readers::check_missing_faidx, InputSource, Result};
use rust_htslib::bam::{self, Read};
use std::path::Path;

pub enum ReferencePolicy {
    Never,
    Always,
    AutoForCram,
}

pub struct BamOpenOptions<'a> {
    pub threads: Option<usize>,
    pub reference: Option<&'a InputSource>,
    pub ref_policy: ReferencePolicy,
}

impl Default for BamOpenOptions<'_> {
    fn default() -> Self {
        Self {
            threads: None,
            reference: None,
            ref_policy: ReferencePolicy::Never,
        }
    }
}

fn is_cram_source(src: &InputSource) -> bool {
    match src {
        InputSource::Local(p) => has_ext(p, "cram"),
        InputSource::Remote(r) => r.url().path().to_ascii_lowercase().ends_with(".cram"),
    }
}

fn has_ext(p: &Path, ext: &str) -> bool {
    p.extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case(ext))
        .unwrap_or(false)
}

pub fn create_bam_reader_with_options(
    reads_src: &InputSource,
    opts: &BamOpenOptions<'_>,
) -> Result<bam::IndexedReader> {
    reads_src.preflight_checks()?;
    let mut reader = match reads_src {
        InputSource::Local(p) => bam::IndexedReader::from_path(p),
        InputSource::Remote(r) => bam::IndexedReader::from_url(r.url()),
    }
    .map_err(|e| reads_src.format_error("Failed to create BAM reader from", e))?;

    if let Some(n) = opts.threads {
        if let Err(e) = reader.set_threads(n) {
            log::warn!("Failed to set decompression threads: {e}");
        }
    }

    let need_ref = match opts.ref_policy {
        ReferencePolicy::Never => false,
        ReferencePolicy::Always => true,
        ReferencePolicy::AutoForCram => is_cram_source(reads_src),
    };
    if need_ref {
        if let Some(ref_src) = opts.reference {
            if let InputSource::Local(p) = ref_src {
                check_missing_faidx(p)?;
                if let Err(e) = reader.set_reference(p) {
                    log::warn!("Failed to set reference {}: {e}", p.display());
                }
            } else {
                log::warn!("Reference is remote; skipping set_reference()");
            }
        }
    }
    Ok(reader)
}

pub fn open_genotyper_bam_reader(
    reads_src: &InputSource,
    genome_src: &InputSource,
    decompression_threads: usize,
) -> Result<bam::IndexedReader> {
    let opts = BamOpenOptions {
        threads: Some(decompression_threads),
        reference: Some(genome_src),
        ref_policy: ReferencePolicy::AutoForCram,
    };
    create_bam_reader_with_options(reads_src, &opts)
}

pub fn open_bam_reader(reads_src: &InputSource, threads: usize) -> Result<bam::IndexedReader> {
    let opts = BamOpenOptions {
        threads: Some(threads),
        ..Default::default()
    };
    create_bam_reader_with_options(reads_src, &opts)
}
