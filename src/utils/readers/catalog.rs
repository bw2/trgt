use crate::utils::{InputSource, Result};
use rust_htslib::{self, bgzf};
use std::io::BufReader;

pub type CatalogReader = BufReader<bgzf::Reader>;
const BUFFER_CAPACITY: usize = 128 * 1024;

pub fn open_catalog_reader(src: &InputSource) -> Result<CatalogReader> {
    src.preflight_checks()?;
    let inner = match src {
        InputSource::Local(p) => bgzf::Reader::from_path(p),
        InputSource::Remote(r) => bgzf::Reader::from_url(r.url()),
    }
    .map_err(|e| src.format_error("Failed to open catalog from", e))?;
    Ok(BufReader::with_capacity(BUFFER_CAPACITY, inner))
}
