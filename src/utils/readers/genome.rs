use crate::utils::{InputSource, Result};
use rust_htslib::faidx;
use std::path::Path;

pub fn open_genome_reader(genome_src: &InputSource) -> Result<faidx::Reader> {
    genome_src.preflight_checks()?;
    if let InputSource::Local(p) = genome_src {
        check_missing_faidx(p)?;
    }

    let reader = match genome_src {
        InputSource::Local(p) => faidx::Reader::from_path(p),
        InputSource::Remote(r) => faidx::Reader::from_url(r.url()),
    }
    .map_err(|e| genome_src.format_error("Failed to open genome file from", e))?;
    Ok(reader)
}

pub fn check_missing_faidx(fasta: &Path) -> Result<()> {
    let ext = fasta.extension().and_then(|s| s.to_str()).unwrap_or("");
    let fai = fasta.with_extension(format!("{ext}.fai"));
    if !fai.exists() {
        return Err(format!(
            "Reference index not found: {}. Create it with 'samtools faidx {}'",
            fai.display(),
            fasta.display()
        ));
    }
    let lower = fasta.to_string_lossy().to_ascii_lowercase();
    if lower.ends_with(".gz") || lower.ends_with(".bgz") || lower.ends_with(".bgzip") {
        let gzi = fasta.with_extension(format!("{ext}.gzi"));
        if !gzi.exists() {
            return Err(format!(
                "Compressed FASTA appears to be missing its .gzi index: {}. \
                 Create it with 'samtools faidx {}'",
                gzi.display(),
                fasta.display()
            ));
        }
    }
    Ok(())
}
