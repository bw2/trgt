use super::align::SegType;
use super::align_allele::get_allele_align;
use super::params::{base_color, get_meth_colors, Color, ColorMap, PlotParams};
use super::scale::get_multi_scale;
use crate::trvz::align::{Align, AlignOp};
use crate::utils::{locus::InputLocus, read::Betas, read::Read};
use itertools::Itertools;
use pipeplot::{Band, FontConfig, GridLine, Legend, Pipe, PipePlot, Seg, Shape, TextLabel};

/// Map sequence positions to alignment (pipe) positions
/// Returns a vector where result[seq_idx] = pipe_x_position
fn build_seq_to_pipe_map(align: &Align) -> Vec<u32> {
    let mut map = Vec::new();
    let mut pipe_x: u32 = 0;

    for seg in align {
        match seg.op {
            AlignOp::Match | AlignOp::Subst | AlignOp::Del => {
                // For consensus align, Del == HMM Ins: consumes sequence and pipe
                for _ in 0..seg.width {
                    map.push(pipe_x);
                    pipe_x += 1;
                }
            }
            AlignOp::Ins => {
                // Insertions in alignment consume sequence but no pipe (width is 0 today)
            }
        }
    }

    map
}

/// Generate text labels for the consensus sequence (all bases)
/// Uses alignment to map sequence positions to pipe positions
fn get_consensus_labels(allele_seq: &[u8], align: &Align) -> Vec<TextLabel> {
    let seq_to_pipe = build_seq_to_pipe_map(align);

    allele_seq
        .iter()
        .enumerate()
        .filter_map(|(seq_idx, &base)| {
            // Only emit label if we have a valid pipe position
            seq_to_pipe.get(seq_idx).map(|&pipe_pos| TextLabel {
                pos: pipe_pos,
                text: base as char,
                color: base_color(base).to_string(),
            })
        })
        .collect()
}

pub fn plot_alleles(
    locus: &InputLocus,
    what_to_show: &str,
    allele_seqs: &[Vec<u8>],
    reads: &[Read],
    params: PlotParams,
) -> PipePlot {
    let aligns_by_allele: Vec<_> = allele_seqs
        .iter()
        .enumerate()
        .map(|(index, allele_seq)| {
            let allele_reads = reads
                .iter()
                .filter(|r| r.allele == index as i32)
                .collect_vec();
            (allele_seq, get_allele_align(locus, allele_seq, &allele_reads))
        })
        .collect();
    let allele_height = 4;
    let xpos = 0;
    let mut ypos = 0;
    let mut pipes = Vec::new();
    let mut grid_lines = Vec::new();
    let num_motifs = locus.motifs.len();

    // Get the primary motif length for grid spacing
    let motif_len = locus.motifs.first().map(|m| m.len()).unwrap_or(3);

    for (allele_index, (allele_seq, allele)) in aligns_by_allele.iter().enumerate() {
        let (scale_pipes, scale_height) = get_multi_scale(xpos, ypos, allele_height, &allele.seq, locus);
        pipes.extend(scale_pipes);

        ypos += scale_height;  // After scale row(s)
        ypos += allele_height;  // Padding between scale and consensus

        // Track y_start for grid lines (at consensus row)
        let y_start = ypos;

        // Grid labels positioned just above consensus (bottom of label touches top of consensus)
        let label_y = y_start;

        let consensus_labels = get_consensus_labels(allele_seq, &allele.seq);
        let pipe = get_pipe(
            xpos,
            ypos,
            allele_height,
            &allele.seq,
            &Vec::new(),
            &params.colors,
            true,
            consensus_labels,
            num_motifs,
        );
        pipes.push(pipe);
        ypos += allele_height + params.pipe_pad;

        // Add extra padding if pipes are bookended
        if params.pipe_pad == 0 {
            ypos += 1;
        }

        // TODO: Confirm that this allele / index correspondence is always correct
        for (align, betas, labels) in &allele.reads {
            let (colors, betas, labels) = if what_to_show == "meth" {
                (get_meth_colors(&locus.motifs), betas.clone(), Vec::new())
            } else {
                (params.colors.clone(), Vec::new(), labels.clone())
            };

            let pipe = get_pipe(
                xpos,
                ypos,
                params.pipe_height,
                align,
                &betas,
                &colors,
                false,
                labels,
                num_motifs,
            );
            pipes.push(pipe);
            ypos += params.pipe_height + params.pipe_pad;
        }

        // Calculate y_end for grid lines (after all reads)
        let y_end = ypos;

        // Calculate grid line positions for individual motif runs
        // Collect individual motif runs: (start_x, length, motif_idx)
        let mut motif_runs: Vec<(u32, u32, usize)> = Vec::new();
        let mut current_x: u32 = 0;
        let mut current_motif: Option<usize> = None;
        let mut run_start_x: u32 = 0;
        let mut run_length: u32 = 0;

        for seg in allele.seq.iter() {
            match seg.seg_type {
                SegType::LeftFlank | SegType::RightFlank => {
                    // Flush any pending motif run
                    if let Some(motif_idx) = current_motif {
                        if run_length > 0 {
                            motif_runs.push((run_start_x, run_length, motif_idx));
                        }
                        current_motif = None;
                        run_length = 0;
                    }
                    current_x += seg.width as u32;
                }
                SegType::Tr(motif_idx) => {
                    if motif_idx < num_motifs {
                        // Recognized motif
                        if current_motif == Some(motif_idx) {
                            run_length += seg.width as u32;
                        } else {
                            // Flush previous run if any
                            if let Some(prev_idx) = current_motif {
                                if run_length > 0 {
                                    motif_runs.push((run_start_x, run_length, prev_idx));
                                }
                            }
                            run_start_x = current_x;
                            run_length = seg.width as u32;
                            current_motif = Some(motif_idx);
                        }
                    } else {
                        // Unsegmented region
                        if let Some(prev_idx) = current_motif {
                            if run_length > 0 {
                                motif_runs.push((run_start_x, run_length, prev_idx));
                            }
                            current_motif = None;
                            run_length = 0;
                        }
                    }
                    current_x += seg.width as u32;
                }
            }
        }
        // Flush final run
        if let Some(motif_idx) = current_motif {
            if run_length > 0 {
                motif_runs.push((run_start_x, run_length, motif_idx));
            }
        }

        // Add grid lines for each run >= 15 repeats
        for (run_start, run_len, motif_idx) in &motif_runs {
            let this_motif_len = locus.motifs.get(*motif_idx).map(|m| m.len()).unwrap_or(motif_len);
            let run_repeats = *run_len as usize / this_motif_len;

            // Only add labels for runs >= 15 repeats
            if run_repeats < 15 {
                continue;
            }

            // Calculate step size to avoid label overlap
            // Estimate: label "XXXx" is ~4 chars * 8 pixels = ~32 pixels + 15px gap = ~47px per label
            // At x_scale = 2500/longest_pipe, each bp is roughly 2500/longest_pipe pixels
            // Use allele width as estimate for longest_pipe
            let allele_width: u32 = allele.seq.iter().map(|seg| seg.width as u32).sum();
            let pixels_per_bp = if allele_width > 0 { 2500.0 / allele_width as f64 } else { 1.0 };
            let run_width_pixels = *run_len as f64 * pixels_per_bp;
            let min_pixels_per_label = 50.0; // Minimum spacing to avoid overlap

            // Calculate how many labels can fit
            let max_labels = (run_width_pixels / min_pixels_per_label).floor() as usize;

            // Choose step size: start at 5, double until labels fit
            let mut step = 5usize;
            while step < run_repeats {
                let num_labels = run_repeats / step;
                if num_labels <= max_labels.max(1) {
                    break;
                }
                step *= 2;
            }

            // Add grid lines for this run, starting from 0
            let bp_per_step = (step * this_motif_len) as u32;
            let mut repeat_count = step;
            let mut grid_x = *run_start + bp_per_step;
            while grid_x < *run_start + *run_len {
                grid_lines.push(GridLine {
                    xpos: grid_x,
                    y_start,
                    y_end,
                    label_y,
                    label: Some(format!("{}x", repeat_count)),
                });
                repeat_count += step;
                grid_x += bp_per_step;
            }
        }

        if allele_index + 1 != aligns_by_allele.len() {
            ypos += 7;
        }
    }

    let mut labels = Vec::new();
    if what_to_show == "meth" {
        // Only show methylation legend for meth plots
        labels.push(("Methylated".to_string(), Color::Grad(1.0).to_string()));
        labels.push(("Unmethylated".to_string(), Color::Grad(0.0).to_string()));
    } else {
        // Show motif legend for non-meth plots
        for (index, motif) in locus.motifs.iter().enumerate() {
            let color = params.colors.get(&SegType::Tr(index)).unwrap().to_string();
            labels.push((std::str::from_utf8(motif).unwrap().to_string(), color));
        }
    }

    ypos += 1;

    let legend = Legend {
        xpos,
        ypos,
        height: allele_height,
        labels,
    };

    PipePlot {
        pipes,
        legend,
        font: FontConfig::default(),
        grid_lines,
    }
}

fn get_pipe(
    xpos: u32,
    ypos: u32,
    height: u32,
    align: &Align,
    betas: &Betas,
    colors: &ColorMap,
    outline: bool,
    labels: Vec<TextLabel>,
    num_motifs: usize,
) -> Pipe {
    let mut segs = Vec::new();
    for align_seg in align {
        let shape = match align_seg.op {
            AlignOp::Del => Shape::HLine,
            AlignOp::Ins => Shape::VLine(align_seg.insertion_size as u32),
            AlignOp::Match | AlignOp::Subst => Shape::Rect,
        };
        let color = if align_seg.op == AlignOp::Match {
            colors.get(&align_seg.seg_type).unwrap()
        } else if align_seg.op == AlignOp::Subst {
            &Color::Gray
        } else {
            &Color::Black
        };
        let dashed = match align_seg.seg_type {
            SegType::Tr(idx) if idx == num_motifs => true,
            _ => false,
        };
        segs.push(Seg {
            width: align_seg.width as u32,
            color: color.to_string(),
            shape,
            dashed,
        });
    }

    let mut bands = Vec::new();

    for beta in betas {
        // Band width is 2 to cover the entire CpG
        let color = Color::Grad(beta.value);
        bands.push(Band {
            pos: beta.pos as u32,
            width: 2,
            color: color.to_string(),
        });
    }

    Pipe {
        xpos,
        ypos,
        height,
        segs,
        bands,
        outline,
        labels,
    }
}
