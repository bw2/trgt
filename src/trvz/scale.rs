use super::align::{Align, SegType};
use super::params::Color;
use crate::utils::locus::InputLocus;
use pipeplot::{Pipe, Seg, Shape};
use std::collections::HashMap;

pub fn get_scale(mut xpos: u32, ypos: u32, height: u32, align: &Align) -> Pipe {
    let lf_len = align
        .iter()
        .filter(|op| op.seg_type == SegType::LeftFlank)
        .map(|op| op.width)
        .sum::<usize>();
    xpos += lf_len as u32;

    let allele_len = align
        .iter()
        .filter(|op| op.seg_type != SegType::LeftFlank && op.seg_type != SegType::RightFlank)
        .map(|op| op.width)
        .sum::<usize>();

    let mut label = allele_len.to_string();
    label += "bp";
    let segs = vec![Seg {
        width: allele_len as u32,
        color: Color::Black.to_string(),
        shape: Shape::DoubleArrow(Some(label)),
        dashed: false,
    }];

    Pipe {
        xpos,
        ypos,
        height,
        segs,
        bands: Vec::new(),
        outline: false,
        labels: Vec::new(),
    }
}

/// Generate scale bars for each contiguous motif region
/// Iterates through segments in order, advancing position for ALL segment types
/// When a motif has multiple runs, adds an overarching bp label above the individual labels
/// Returns (pipes, height_used) where height_used is height or 2*height if overarching label exists
pub fn get_multi_scale(
    xpos: u32,
    ypos: u32,
    height: u32,
    align: &Align,
    locus: &InputLocus,
) -> (Vec<Pipe>, u32) {
    let mut current_x = xpos;

    // Collect runs by motif index: motif_idx -> [(start_x, length), ...]
    let mut motif_runs: HashMap<usize, Vec<(u32, usize)>> = HashMap::new();

    // Track contiguous runs of recognized motifs only (not unsegmented)
    let mut current_motif: Option<usize> = None;
    let mut run_start_x = current_x;
    let mut run_length: usize = 0;

    for seg in align.iter() {
        match seg.seg_type {
            SegType::LeftFlank | SegType::RightFlank => {
                // Flush any pending motif run
                if let Some(motif_idx) = current_motif {
                    if run_length > 0 {
                        motif_runs
                            .entry(motif_idx)
                            .or_default()
                            .push((run_start_x, run_length));
                    }
                    current_motif = None;
                    run_length = 0;
                }
                // Advance position for flank
                current_x += seg.width as u32;
            }
            SegType::Tr(motif_idx) => {
                if motif_idx < locus.motifs.len() {
                    // Recognized motif
                    if current_motif == Some(motif_idx) {
                        // Continue current run
                        run_length += seg.width;
                    } else {
                        // Flush previous run if any
                        if let Some(prev_idx) = current_motif {
                            if run_length > 0 {
                                motif_runs
                                    .entry(prev_idx)
                                    .or_default()
                                    .push((run_start_x, run_length));
                            }
                        }
                        // Start new run
                        run_start_x = current_x;
                        run_length = seg.width;
                        current_motif = Some(motif_idx);
                    }
                } else {
                    // Unsegmented region (motif_idx == motifs.len())
                    // Flush any pending motif run
                    if let Some(prev_idx) = current_motif {
                        if run_length > 0 {
                            motif_runs
                                .entry(prev_idx)
                                .or_default()
                                .push((run_start_x, run_length));
                        }
                        current_motif = None;
                        run_length = 0;
                    }
                    // No scale bar for unsegmented, but advance position
                }
                // Always advance position for TR segments
                current_x += seg.width as u32;
            }
        }
    }

    // Flush final run if any
    if let Some(motif_idx) = current_motif {
        if run_length > 0 {
            motif_runs
                .entry(motif_idx)
                .or_default()
                .push((run_start_x, run_length));
        }
    }

    // Check if any motif has multiple runs (needs overarching label)
    let has_overarching = motif_runs.values().any(|runs| runs.len() > 1);

    // If overarching label needed: place it at ypos, individual labels at ypos + height
    // Otherwise: place individual labels at ypos
    let individual_ypos = if has_overarching { ypos + height } else { ypos };
    let height_used = if has_overarching { 2 * height } else { height };

    // Calculate pixel scale for label width check
    // Use conservative estimate (cap at 1.5) since actual scale depends on longest pipe in plot
    let total_width: u32 = align.iter().map(|seg| seg.width as u32).sum();
    let pixels_per_bp = if total_width > 0 {
        (2500.0 / total_width as f64).min(1.5)
    } else {
        1.0
    };

    // Create pipes from collected runs
    let mut pipes = Vec::new();

    for (motif_idx, runs) in &motif_runs {
        // Add individual motif annotations (only if label fits within arrows)
        for (start_x, length) in runs {
            if let Some(pipe) = create_scale_pipe_if_fits(
                *start_x,
                individual_ypos,
                height,
                *length,
                &locus.motifs[*motif_idx],
                pixels_per_bp,
            ) {
                pipes.push(pipe);
            }
        }

        // If more than one run, add overarching bp label at the top row
        if runs.len() > 1 {
            let min_x = runs.iter().map(|(x, _)| *x).min().unwrap();
            let max_end = runs.iter().map(|(x, len)| *x + *len as u32).max().unwrap();
            let total_bp: usize = runs.iter().map(|(_, len)| *len).sum();
            let span_width = max_end - min_x;

            pipes.push(create_overarching_bp_pipe(
                min_x,
                ypos,
                height,
                span_width,
                total_bp,
            ));
        }
    }

    (pipes, height_used)
}

/// Create a scale pipe only if the label fits within the arrow span.
/// Returns None if the label would be wider than the arrows (completely covering them).
fn create_scale_pipe_if_fits(
    xpos: u32,
    ypos: u32,
    height: u32,
    region_len: usize,
    motif: &[u8],
    pixels_per_bp: f64,
) -> Option<Pipe> {
    let motif_str = std::str::from_utf8(motif).unwrap_or("?");
    let repeat_count = region_len as f64 / motif.len() as f64;

    // Format: "31 x CAG" or "12.3 x CCG"
    let label = if (repeat_count - repeat_count.round()).abs() < 0.05 {
        format!("{} x {}", repeat_count.round() as u32, motif_str)
    } else {
        format!("{:.1} x {}", repeat_count, motif_str)
    };

    // Calculate label width in pixels (matching svg.rs: ~10px per char + 8px padding)
    let label_width_px = label.len() as f64 * 10.0 + 8.0;

    // Calculate arrow span width in pixels
    let arrow_width_px = region_len as f64 * pixels_per_bp;

    // Only show if arrows extend past the label (label fits within arrows)
    if label_width_px >= arrow_width_px {
        return None;
    }

    let seg = Seg {
        width: region_len as u32,
        color: Color::Black.to_string(),
        shape: Shape::DoubleArrow(Some(label)),
        dashed: false,
    };

    Some(Pipe {
        xpos,
        ypos,
        height,
        segs: vec![seg],
        bands: Vec::new(),
        labels: Vec::new(),
        outline: false,
    })
}

/// Create an overarching bp label that spans multiple motif runs
fn create_overarching_bp_pipe(
    xpos: u32,
    ypos: u32,
    height: u32,
    span_width: u32,
    total_bp: usize,
) -> Pipe {
    let label = format!("{}bp", total_bp);

    let seg = Seg {
        width: span_width,
        color: Color::Black.to_string(),
        shape: Shape::DoubleArrow(Some(label)),
        dashed: false,
    };

    Pipe {
        xpos,
        ypos,
        height,
        segs: vec![seg],
        bands: Vec::new(),
        labels: Vec::new(),
        outline: false,
    }
}
