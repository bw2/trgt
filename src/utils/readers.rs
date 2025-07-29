use super::Result;
use flate2::read::MultiGzDecoder;
use rust_htslib::{
    bam::{self, Read},
    faidx,
};
use std::{
    fs::File,
    io::{self, BufRead, BufReader, Read as IoRead},
    path::Path,
};

enum CatalogReader {
    Gzipped(MultiGzDecoder<File>),
    Plain(File),
}

impl IoRead for CatalogReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            CatalogReader::Gzipped(reader) => reader.read(buf),
            CatalogReader::Plain(reader) => reader.read(buf),
        }
    }
}

pub fn open_catalog_reader(path: &Path) -> Result<impl BufRead> {
    fn is_gzipped(path: &Path) -> bool {
        let path_str = path.to_string_lossy().to_lowercase();
        path_str.ends_with(".gz") || path_str.ends_with(".gzip")
    }
    let file = File::open(path).map_err(|e| e.to_string())?;

    let reader: CatalogReader = if is_gzipped(path) {
        let gz_decoder = MultiGzDecoder::new(file);
        if gz_decoder.header().is_some() {
            CatalogReader::Gzipped(gz_decoder)
        } else {
            return Err(format!("Invalid gzip header: {}", path.to_string_lossy()));
        }
    } else {
        CatalogReader::Plain(file)
    };
    const BUFFER_CAPACITY: usize = 128 * 1024;
    Ok(BufReader::with_capacity(BUFFER_CAPACITY, reader))
}

pub fn open_genome_reader(path: &Path) -> Result<faidx::Reader> {
    let extension = path.extension().unwrap().to_str().unwrap();
    let fai_path = path.with_extension(extension.to_owned() + ".fai");
    if !fai_path.exists() {
        return Err(format!(
            "Reference index file not found: {}. Create it using 'samtools faidx {}'",
            fai_path.display(),
            path.display()
        ));
    }
    faidx::Reader::from_path(path).map_err(|e| e.to_string())
}

pub fn create_bam_reader(
    reads_path: &Path,
    genome_path: &Path,
    decompression_threads: usize,
) -> Result<bam::IndexedReader> {
    let mut reader = bam::IndexedReader::from_path(reads_path)
        .map_err(|e| format!("Failed to create BAM reader: {}", e))?;

    if let Err(e) = reader.set_threads(decompression_threads) {
        log::warn!("Failed to set decompression threads: {}", e);
    }

    if let Err(e) = reader.set_reference(genome_path) {
        log::warn!("Failed to set reference: {}", e);
    }

    Ok(reader)
}
