mod cigar;
mod clip_bases;
mod clip_region;
mod snp;
pub use snp::extract_mismatch_offsets;

mod read;
pub use read::get_rq_tag;
pub use read::HiFiRead;
