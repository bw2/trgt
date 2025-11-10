mod consensus;
mod diploid;
pub mod genotype_cluster;
pub mod genotype_flank;
pub mod genotype_size;
mod gt;
mod haploid;
mod span_locater;

pub use gt::{Gt, TrSize};
pub use span_locater::find_tr_spans;
