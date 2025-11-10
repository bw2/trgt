use clap::Parser;
use std::{env, time};
use trgt::{
    cli::{init_verbose, Cli, Command, FULL_VERSION},
    commands::{deepdive, genotype, merge, plot, validate},
    utils::{handle_error_and_exit, Result},
};

fn disable_htslib_logging() {
    if env::var_os("TRGT_ENABLE_HTSLIB_LOGGING").is_some() {
        log::debug!("TRGT_ENABLE_HTSLIB_LOGGING is set, keeping htslib logging enabled");
        return;
    }
    unsafe {
        use rust_htslib::htslib::{htsLogLevel_HTS_LOG_OFF, hts_set_log_level};
        hts_set_log_level(htsLogLevel_HTS_LOG_OFF);
    }
}

fn runner() -> Result<()> {
    let cli = Cli::parse();
    init_verbose(&cli);
    disable_htslib_logging();
    log::info!(
        "Running {}-{} [{}]",
        env!("CARGO_PKG_NAME"),
        FULL_VERSION,
        cli.command.name()
    );

    let start_timer = time::Instant::now();
    match cli.command {
        Command::Genotype(args) => {
            log::trace!("Genotype arguments: {:#?}", args);
            args.preflight()?;
            genotype::trgt(args)?
        }
        Command::Plot(args) => {
            log::trace!("Plot arguments: {:#?}", args);
            args.preflight()?;
            plot::trvz(args)?
        }
        Command::Validate(args) => {
            log::trace!("Validate arguments: {:#?}", args);
            args.preflight()?;
            validate::validate(args)?
        }
        Command::Merge(args) => {
            log::trace!("Merge arguments: {:#?}", args);
            args.preflight()?;
            merge::merge(args)?
        }
        Command::Deepdive(args) => {
            log::trace!("Deep dive arguments: {:#?}", args);
            args.preflight()?;
            deepdive::deepdive(args)?
        }
    }

    log::info!("Total execution time: {:.2?}", start_timer.elapsed());
    log::info!("{} end", env!("CARGO_PKG_NAME"));
    Ok(())
}

fn main() {
    if let Err(e) = runner() {
        handle_error_and_exit(e);
    }
}
