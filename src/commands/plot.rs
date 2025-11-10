use crate::cli::PlotArgs;
use crate::trvz::allele_plot::plot_alleles;
use crate::trvz::params::pick_params;
use crate::trvz::waterfall_plot::plot_waterfall;
use crate::utils::{
    input::{get_alleles, get_locus, get_reads},
    open_catalog_reader, open_genome_reader, Result,
};
use pipeplot::generate_image;

pub fn trvz(args: PlotArgs) -> Result<()> {
    let catalog_reader = open_catalog_reader(&args.repeats_src)?;
    let genome_reader = open_genome_reader(&args.genome_src)?;
    let locus = get_locus(catalog_reader, genome_reader, &args.tr_id, args.flank_len)?;
    let reads = get_reads(&args.reads_src, &locus, args.max_allele_reads)?;
    let params = pick_params(&locus.motifs, args.is_squished);
    let mut pipe_plot = if args.plot_type == "allele" {
        let allele_seqs = get_alleles(&args.bcf_src, &locus)?;
        plot_alleles(&locus, &args.what_to_show, &allele_seqs, &reads, params)
    } else {
        plot_waterfall(&locus, &args.what_to_show, &reads, &params)
    };

    if let Some(font_family) = args.font_family {
        pipe_plot.set_font_family(&font_family);
    }

    generate_image(&pipe_plot, &args.output_path)?;
    Ok(())
}
