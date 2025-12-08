//! A thread pool for the rust_htslib::bcf module
//!

// NOTE: There's https://github.com/rust-bio/rust-htslib/blob/master/src/tpool.rs but there's no methods that can use it in the bcf module (bgzf::Reader, bgzf::Writer, bam::Reader, bam::IndexedReader all support it)
//       probably worth extending to cover bcf so we do not have to do this work below:

use crate::utils::Result;
use rust_htslib::{bcf, htslib};
use std::ptr::NonNull;

/// Thread pool for parallel BGZF decompression across multiple VCF readers.
/// This allows sharing decompression threads among multiple BCF IndexedReader's, Reader's, and Writer's.
///
/// # Safety
/// 1. Assumes rust-htslib IndexedReader, Reader, Writer layout as of 0.50.x.
/// 2. The pool must outlive all future uses of Reader, Writer, IndexedReader.
///
/// Note: The const assertions below do not fully guarantee layout compatibility; they only detect some incompatible changes in rust-htslib updates.
pub struct HtsThreadPool {
    pool: NonNull<htslib::htsThreadPool>,
}

impl HtsThreadPool {
    pub fn new(n_threads: i32) -> Result<Self> {
        unsafe {
            let mut inner = Box::new(htslib::htsThreadPool {
                pool: std::ptr::null_mut(),
                qsize: n_threads * 2,
            });

            inner.pool = htslib::hts_tpool_init(n_threads);
            if inner.pool.is_null() {
                return Err("Failed to initialize thread pool".to_string());
            }

            Ok(Self {
                pool: NonNull::new(Box::into_raw(inner)).expect("Box::into_raw never returns null"),
            })
        }
    }

    pub unsafe fn attach_to_indexed_reader(&self, reader: &mut bcf::IndexedReader) -> Result<()> {
        unsafe {
            #[repr(C)]
            struct IndexedReaderLayout {
                inner: *mut htslib::bcf_srs_t,
                header: std::rc::Rc<bcf::header::HeaderView>,
                current_region: Option<(u32, u64, Option<u64>)>,
            }

            const _: () = {
                assert!(
                    std::mem::size_of::<bcf::IndexedReader>()
                        == std::mem::size_of::<IndexedReaderLayout>()
                );
                assert!(
                    std::mem::align_of::<bcf::IndexedReader>()
                        == std::mem::align_of::<IndexedReaderLayout>()
                );
            };

            let reader_ptr = reader as *mut bcf::IndexedReader as *mut IndexedReaderLayout;

            debug_assert!(
                reader_ptr.is_aligned(),
                "IndexedReader pointer is not properly aligned for IndexedReaderLayout cast"
            );

            let bcf_srs_ptr = (*reader_ptr).inner;
            if bcf_srs_ptr.is_null() {
                return Err("IndexedReader inner pointer is null".to_string());
            }

            if (*bcf_srs_ptr).nreaders != 1 {
                let reader_count = (*bcf_srs_ptr).nreaders;
                return Err(format!(
                    "Expected exactly 1 reader in bcf_srs_t, got {reader_count}"
                ));
            }

            let bcf_sr_ptr = (*bcf_srs_ptr).readers;
            if bcf_sr_ptr.is_null() {
                return Err("Readers array is null".to_string());
            }

            let hts_file_ptr = (*bcf_sr_ptr).file;
            self.attach_to_hts_file(hts_file_ptr, "bcf::IndexedReader")?;
        }

        Ok(())
    }

    pub unsafe fn attach_to_reader(&self, reader: &mut bcf::Reader) -> Result<()> {
        #[repr(C)]
        struct ReaderLayout {
            inner: *mut htslib::htsFile,
            header: std::rc::Rc<bcf::header::HeaderView>,
        }

        const _: () = {
            assert!(std::mem::size_of::<bcf::Reader>() == std::mem::size_of::<ReaderLayout>());
            assert!(std::mem::align_of::<bcf::Reader>() == std::mem::align_of::<ReaderLayout>());
        };

        let reader_ptr = reader as *mut bcf::Reader as *mut ReaderLayout;

        debug_assert!(
            reader_ptr.is_aligned(),
            "Reader pointer is not properly aligned for ReaderLayout cast"
        );

        let hts_file_ptr = (*reader_ptr).inner;

        self.attach_to_hts_file(hts_file_ptr, "bcf::Reader")
    }

    pub unsafe fn attach_to_writer(&self, writer: &mut bcf::Writer) -> Result<()> {
        #[repr(C)]
        struct WriterLayout {
            subset: Option<bcf::header::SampleSubset>,
            header: std::rc::Rc<bcf::header::HeaderView>,
            inner: *mut htslib::htsFile,
        }

        const _: () = {
            assert!(std::mem::size_of::<bcf::Writer>() == std::mem::size_of::<WriterLayout>());
            assert!(std::mem::align_of::<bcf::Writer>() == std::mem::align_of::<WriterLayout>());
        };

        let writer_ptr = writer as *mut bcf::Writer as *mut WriterLayout;

        debug_assert!(
            writer_ptr.is_aligned(),
            "Writer pointer is not properly aligned for WriterLayout cast"
        );

        debug_assert!(
            std::rc::Rc::as_ptr(&(*writer_ptr).header) == writer.header() as *const _,
            "Writer header pointer layout mismatch"
        );

        let hts_file_ptr = (*writer_ptr).inner;

        self.attach_to_hts_file(hts_file_ptr, "bcf::Writer")
    }

    unsafe fn attach_to_hts_file(
        &self,
        hts_file_ptr: *mut htslib::htsFile,
        target: &str,
    ) -> Result<()> {
        if hts_file_ptr.is_null() {
            return Err(format!("{target} htsFile pointer is null"));
        }

        let ret = htslib::hts_set_thread_pool(hts_file_ptr, self.pool.as_ptr());
        if ret != 0 {
            return Err(format!(
                "Failed to set thread pool on {target} (error code: {ret})"
            ));
        }

        let pool_ref = self.pool.as_ref();
        let thread_count = if pool_ref.pool.is_null() {
            0
        } else {
            htslib::hts_tpool_size(pool_ref.pool)
        };

        log::trace!("Attached thread pool with {thread_count} threads to {target}");

        Ok(())
    }
}

impl Drop for HtsThreadPool {
    fn drop(&mut self) {
        unsafe {
            let raw = self.pool.as_ptr();
            if !(*raw).pool.is_null() {
                htslib::hts_tpool_destroy((*raw).pool);
            }
            drop(Box::from_raw(raw));
        }
    }
}

unsafe impl Send for HtsThreadPool {}
unsafe impl Sync for HtsThreadPool {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::test_util::TestVcfBuilder;
    use rust_htslib::bcf::{self, record::GenotypeAllele, Read};
    use tempfile::tempdir;

    #[test]
    fn test_attach_to_reader_is_callable() -> Result<()> {
        let pool = HtsThreadPool::new(2)?;
        let dir = tempdir().map_err(|e| e.to_string())?;
        let path = dir.path().join("reader_attach.vcf.gz");

        let header = TestVcfBuilder::new()
            .contig("chr1", 16)
            .sample("sample")
            .build_header();

        {
            let mut writer = bcf::Writer::from_path(&path, &header, false, bcf::Format::Vcf)
                .map_err(|e| {
                    let display_path = path.to_string_lossy();
                    format!("failed to create writer at {display_path}: {e}")
                })?;
            let mut record = writer.empty_record();
            let rid = writer
                .header()
                .name2rid(b"chr1")
                .map_err(|e| e.to_string())?;
            record.set_rid(Some(rid));
            record.set_pos(0);
            record
                .set_alleles(&[b"A", b"T"])
                .map_err(|e| e.to_string())?;
            record
                .push_genotypes(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
                .map_err(|e| e.to_string())?;
            writer.write(&record).map_err(|e| e.to_string())?;
        }

        bcf::index::build(&path, None::<&std::path::PathBuf>, 1, bcf::index::Type::Tbx)
            .map_err(|e| e.msg)?;

        let mut reader = bcf::Reader::from_path(&path).map_err(|e| e.to_string())?;
        unsafe { pool.attach_to_reader(&mut reader)? };
        let mut record = reader.empty_record();
        assert!(matches!(reader.read(&mut record), Some(Ok(()))));

        let mut indexed_reader = bcf::IndexedReader::from_path(&path).map_err(|e| e.to_string())?;
        unsafe { pool.attach_to_indexed_reader(&mut indexed_reader)? };
        let mut indexed_record = indexed_reader.empty_record();
        assert!(matches!(
            indexed_reader.read(&mut indexed_record),
            Some(Ok(()))
        ));

        Ok(())
    }

    #[test]
    fn attach_to_writer_is_callable() -> Result<()> {
        let pool = HtsThreadPool::new(2)?;
        let header = TestVcfBuilder::new()
            .contig("chr1", 16)
            .sample("sample")
            .build_header();

        let dir = tempdir().map_err(|e| e.to_string())?;
        let path = dir.path().join("writer_attach.vcf.gz");
        let mut writer =
            bcf::Writer::from_path(&path, &header, false, bcf::Format::Vcf).map_err(|e| {
                let display_path = path.to_string_lossy();
                format!("failed to create writer at {display_path}: {e}")
            })?;
        unsafe { pool.attach_to_writer(&mut writer)? };

        let mut record = writer.empty_record();
        let rid = writer
            .header()
            .name2rid(b"chr1")
            .map_err(|e| e.to_string())?;
        record.set_rid(Some(rid));
        record.set_pos(0);
        record
            .set_alleles(&[b"A", b"T"])
            .map_err(|e| e.to_string())?;
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
            .map_err(|e| e.to_string())?;
        writer.write(&record).map_err(|e| e.to_string())?;
        drop(writer);

        let mut reader = bcf::Reader::from_path(&path).map_err(|e| e.to_string())?;
        let mut read_record = reader.empty_record();
        assert!(matches!(reader.read(&mut read_record), Some(Ok(()))));

        Ok(())
    }

    #[test]
    fn test_attach_to_all() -> Result<()> {
        let pool = HtsThreadPool::new(4)?;
        let dir = tempdir().map_err(|e| e.to_string())?;
        let path = dir.path().join("all_attach.vcf.gz");

        let header = TestVcfBuilder::new()
            .contig("chr1", 1000)
            .sample("sample")
            .build_header();

        // Scoped to make sure the file is closed before the index is built
        {
            let mut writer = bcf::Writer::from_path(&path, &header, false, bcf::Format::Vcf)
                .map_err(|e| {
                    let display_path = path.to_string_lossy();
                    format!("failed to create writer at {display_path}: {e}")
                })?;
            unsafe { pool.attach_to_writer(&mut writer)? };

            let rid = writer
                .header()
                .name2rid(b"chr1")
                .map_err(|e| e.to_string())?;

            for pos in [0, 100, 200] {
                let mut record = writer.empty_record();
                record.set_rid(Some(rid));
                record.set_pos(pos);
                record
                    .set_alleles(&[b"A", b"T"])
                    .map_err(|e| e.to_string())?;
                record
                    .push_genotypes(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
                    .map_err(|e| e.to_string())?;
                writer.write(&record).map_err(|e| e.to_string())?;
            }
        }

        bcf::index::build(&path, None::<&std::path::PathBuf>, 1, bcf::index::Type::Tbx)
            .map_err(|e| e.msg)?;

        // Attach to readers and read all records
        let mut reader = bcf::Reader::from_path(&path).map_err(|e| e.to_string())?;
        unsafe { pool.attach_to_reader(&mut reader)? };

        let mut indexed_reader = bcf::IndexedReader::from_path(&path).map_err(|e| e.to_string())?;
        unsafe { pool.attach_to_indexed_reader(&mut indexed_reader)? };

        let mut record_count = 0;
        let mut record = reader.empty_record();
        while let Some(Ok(())) = reader.read(&mut record) {
            record_count += 1;
        }
        assert_eq!(record_count, 3);

        indexed_reader
            .fetch(indexed_reader.header().name2rid(b"chr1").unwrap(), 0, None)
            .map_err(|e| e.to_string())?;

        let mut indexed_record_count = 0;
        let mut indexed_record = indexed_reader.empty_record();
        while let Some(Ok(())) = indexed_reader.read(&mut indexed_record) {
            indexed_record_count += 1;
        }
        assert_eq!(indexed_record_count, 3);

        Ok(())
    }
}
