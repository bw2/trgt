use crate::cli::GenotypeArgs;
use crate::trgt::{
    locus::{stream_locus_groups_into_channel, Locus},
    locus_group::{extract_reads_for_group, LocusGroup},
    workflows::{self, analyze_tr, LocusResult, Params},
    writers::{BamWriter, VcfWriter},
};
use crate::utils::{
    create_writer, get_bam_header, get_sample_name, is_bam_mapped, readers::create_bam_reader,
    Karyotype, Result, TrgtScoring,
};
use crate::wfaligner::{AlignmentScope, Heuristic, MemoryModel, WFAligner};
use crossbeam_channel::{bounded, Receiver, Sender};
use rayon::{
    iter::{ParallelBridge, ParallelIterator},
    ThreadPoolBuilder,
};
use rust_htslib::bam;
use std::{
    cell::RefCell,
    path::Path,
    sync::Arc,
    thread::{self, JoinHandle},
};

#[derive(Debug, Clone)]
struct ThreadContextParams {
    flank_scoring: TrgtScoring,
}

thread_local! {
    static THREAD_CONTEXT: RefCell<Option<ThreadContextParams>> = const { RefCell::new(None) };
}

fn create_thread_local_ga_aligner_with_scoring() -> WFAligner {
    THREAD_CONTEXT.with(|context_cell| {
        let context = context_cell
            .borrow()
            .as_ref()
            .expect("Thread context parameters not initialized for WFA gap affine aligner")
            .clone();
        let scoring = &context.flank_scoring;

        WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(scoring.mism_scr, scoring.gapo_scr, scoring.gape_scr)
            .with_heuristic(Heuristic::None)
            .build()
    })
}

fn create_thread_local_ga_aligner() -> WFAligner {
    WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
        .affine(2, 5, 1)
        .build()
}

fn create_thread_local_ed() -> WFAligner {
    WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryUltraLow)
        .edit()
        .build()
}

thread_local! {
    // Flank locater aligner
    pub static THREAD_WFA_FLANK: RefCell<WFAligner> = RefCell::new(create_thread_local_ga_aligner_with_scoring());
    // Consensus aligner
    pub static THREAD_WFA_CONSENSUS: RefCell<WFAligner> = RefCell::new(create_thread_local_ga_aligner());
    // Edit distance
    pub static THREAD_WFA_ED: RefCell<WFAligner> = RefCell::new(create_thread_local_ed());
}

pub fn trgt(args: GenotypeArgs) -> Result<()> {
    let thread_context = ThreadContextParams {
        flank_scoring: args.aln_scoring,
    };

    let karyotype = Karyotype::new(&args.karyotype)?;

    let bam_header = get_bam_header(&args.reads_path)?;
    if !is_bam_mapped(&bam_header) {
        return Err("Input BAM is not mapped".into());
    }

    let sample_name = args
        .sample_name
        .unwrap_or(get_sample_name(&args.reads_path, &bam_header)?);

    let mut vcf_writer = create_writer(&args.output_prefix, "vcf.gz", |path| {
        VcfWriter::new(path, &sample_name, &bam_header)
    })?;

    let output_flank_len = std::cmp::min(args.flank_len, args.output_flank_len);
    let bam_writer = if !args.disable_bam_output {
        let n_writer_threads = args.bam_compression_threads.unwrap_or({
            match args.num_threads {
                n if n < 12 => 1,
                _ => 2,
            }
        });
        Some(create_writer(
            &args.output_prefix,
            "spanning.bam",
            |path| BamWriter::new(path, bam_header, output_flank_len, n_writer_threads),
        )?)
    } else {
        None
    };

    let workflow_params = Arc::new(workflows::Params {
        search_flank_len: args.flank_len,
        min_read_qual: args.min_hifi_read_qual,
        max_depth: args.max_depth,
        min_flank_id_frac: args.min_flank_id_frac,
    });

    let num_fetcher_threads = args
        .fetcher_threads
        .unwrap_or_else(|| (args.num_threads / 2).clamp(1, 8));

    let group_channel_buffer_size = args
        .group_channel_buffer_size
        .unwrap_or(num_fetcher_threads);

    // Stage 1 -> 2: Locus group
    let (group_sender, group_receiver) = bounded(group_channel_buffer_size);
    // Stage 2 -> 3: Loci with reads
    let (populated_locus_sender, populated_locus_receiver) =
        bounded(args.locus_channel_buffer_size);
    // Stage 3 -> 4: Final results for writer
    let (analysis_result_sender, analysis_result_receiver) =
        bounded(args.result_channel_buffer_size);

    // Stage 1: Grouping loci into LocusGroup's (IO: Read repeat catalog and reference genome)
    let locus_grouper_thread = {
        let genome_path = args.genome_path.clone();
        thread::spawn(move || {
            stream_locus_groups_into_channel(
                &args.repeats_path,
                &genome_path,
                args.flank_len,
                args.genotyper,
                &karyotype,
                group_sender,
                args.max_group_span,
                args.max_group_size,
            )
        })
    };

    // Stage 2: Populate LocusGroup with reads (IO: Read BAM file)
    let fetcher_threads = spawn_fetcher_threads(
        &args.reads_path,
        &args.genome_path,
        num_fetcher_threads,
        args.decompression_threads,
        &workflow_params,
        group_receiver,
        populated_locus_sender,
    );

    // Stage 4: Spawn the writer thread to write BAM/VCF output (IO: Write BAM and VCF)
    let writer_thread = thread::spawn(move || {
        if let Some(mut bam_writer) = bam_writer {
            for (locus, results) in &analysis_result_receiver {
                vcf_writer.write(&locus, &results);
                bam_writer.write(&locus, &results);
            }
        } else {
            for (locus, results) in &analysis_result_receiver {
                vcf_writer.write(&locus, &results);
            }
        }
    });

    // Stage 3: Worker pool for genotyping TRs
    let pool = initialize_thread_pool(args.num_threads, thread_context)?;
    pool.install(|| {
        populated_locus_receiver
            .into_iter()
            .par_bridge()
            .for_each_with(
                analysis_result_sender,
                |s, locus_result| match locus_result {
                    Ok(locus) => process_locus(locus, &workflow_params, s),
                    Err(err) => log::error!("Locus processing: {:#}", err),
                },
            );
    });

    log::debug!("Analysis complete. Shutting down pipeline threads...");

    match writer_thread.join() {
        Ok(_) => log::trace!("Writer thread has finished."),
        Err(_) => log::error!("Writer thread panicked."),
    }

    match locus_grouper_thread.join() {
        Ok(Ok(_)) => log::trace!("Locus grouper thread has finished."),
        Ok(Err(e)) => log::error!("Locus grouping stage failed: {}", e),
        Err(_) => log::error!("Locus grouper thread panicked."),
    }

    for handle in fetcher_threads {
        match handle.join() {
            Ok(_) => {}
            Err(_) => log::error!("Fetcher thread panicked."),
        }
    }
    log::trace!("All fetcher threads have finished.");
    log::debug!("All pipeline threads shut down successfully.");

    Ok(())
}

fn run_fetcher_thread(
    mut bam_reader: bam::IndexedReader,
    group_receiver: Receiver<Result<LocusGroup>>,
    populated_locus_sender: Sender<Result<Locus>>,
    workflow_params: Arc<Params>,
) {
    for group_result in group_receiver {
        let group = match group_result {
            Ok(group) => group,
            Err(e) => {
                let _ = populated_locus_sender.send(Err(e));
                continue;
            }
        };

        match extract_reads_for_group(&mut bam_reader, group, &workflow_params) {
            Ok(populated_loci) => {
                for locus in populated_loci {
                    if populated_locus_sender.send(Ok(locus)).is_err() {
                        return;
                    }
                }
            }
            Err(e) => {
                let _ = populated_locus_sender.send(Err(e));
            }
        }
    }
}

fn spawn_fetcher_threads(
    reads_path: &Path,
    genome_path: &Path,
    num_fetcher_threads: usize,
    decompression_threads: usize,
    workflow_params: &Arc<Params>,
    group_receiver: Receiver<Result<LocusGroup>>,
    populated_locus_sender: Sender<Result<Locus>>,
) -> Vec<JoinHandle<()>> {
    let mut threads = Vec::with_capacity(num_fetcher_threads);
    log::debug!(
        "Initializing read fetcher pool with {} threads...",
        num_fetcher_threads
    );
    for i in 0..num_fetcher_threads {
        let reads_path = reads_path.to_path_buf();
        let genome_path = genome_path.to_path_buf();
        let group_receiver = group_receiver.clone();
        let populated_locus_sender = populated_locus_sender.clone();
        let workflow_params = workflow_params.clone();

        let handle = thread::Builder::new()
            .name(format!("fetcher-{}", i))
            .spawn(move || {
                let bam_reader =
                    match create_bam_reader(&reads_path, &genome_path, decompression_threads) {
                        Ok(reader) => reader,
                        Err(e) => {
                            log::error!("Fetcher thread failed to initialize: {}", e);
                            return;
                        }
                    };
                run_fetcher_thread(
                    bam_reader,
                    group_receiver,
                    populated_locus_sender,
                    workflow_params,
                );
            })
            .unwrap();
        threads.push(handle);
    }
    threads
}

fn process_locus(
    locus: Locus,
    workflow_params: &Arc<Params>,
    sender_result: &Sender<(Locus, LocusResult)>,
) {
    let (locus, result) = analyze_tr(locus, workflow_params);
    match result {
        Ok(results) => {
            if let Err(e) = sender_result.send((locus, results)) {
                log::error!("Failed to send locus result to writer thread: {}", e);
            }
        }
        Err(err) => {
            log::error!("Error analyzing locus {}: {}", locus.id, err);
        }
    }
}

fn initialize_thread_pool(
    num_threads: usize,
    thread_context: ThreadContextParams,
) -> Result<rayon::ThreadPool> {
    log::debug!(
        "Initializing analysis thread pool with {} threads...",
        num_threads
    );
    ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .thread_name(|i| format!("trgt-{}", i))
        .start_handler(move |_thread_index| {
            THREAD_CONTEXT.with(|cell| {
                *cell.borrow_mut() = Some(thread_context.clone());
            });
            log::trace!("Initialized thread {:?}", std::thread::current().id());
        })
        .exit_handler(|_thread_index| {
            THREAD_CONTEXT.with(|cell| {
                *cell.borrow_mut() = None;
            });
        })
        .build()
        .map_err(|e| format!("Failed to initialize thread pool: {}", e))
}
