use crate::utils::{InputSource, Result};
use rust_htslib::bcf;

pub fn open_vcf_reader(bcf_src: &InputSource) -> Result<bcf::Reader> {
    bcf_src.preflight_checks()?;
    let reader = match bcf_src {
        InputSource::Local(p) => bcf::Reader::from_path(p),
        InputSource::Remote(r) => bcf::Reader::from_url(r.url()),
    }
    .map_err(|e| bcf_src.format_error("Failed to open VCF file from", e))?;
    Ok(reader)
}
