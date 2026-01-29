use super::align::{Align, SegType};
use super::params::Color;
use crate::utils::locus::InputLocus;
use pipeplot::{Pipe, Seg, Shape};

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
pub fn get_multi_scale(
    xpos: u32,
    ypos: u32,
    height: u32,
    align: &Align,
    locus: &InputLocus,
) -> Vec<Pipe> {
    let mut pipes = Vec::new();
    let mut current_x = xpos;

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
                        pipes.push(create_scale_pipe(
                            run_start_x,
                            ypos,
                            height,
                            run_length,
                            &locus.motifs[motif_idx],
                        ));
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
                                pipes.push(create_scale_pipe(
                                    run_start_x,
                                    ypos,
                                    height,
                                    run_length,
                                    &locus.motifs[prev_idx],
                                ));
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
                            pipes.push(create_scale_pipe(
                                run_start_x,
                                ypos,
                                height,
                                run_length,
                                &locus.motifs[prev_idx],
                            ));
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
            pipes.push(create_scale_pipe(
                run_start_x,
                ypos,
                height,
                run_length,
                &locus.motifs[motif_idx],
            ));
        }
    }

    pipes
}

fn create_scale_pipe(xpos: u32, ypos: u32, height: u32, region_len: usize, motif: &[u8]) -> Pipe {
    let motif_str = std::str::from_utf8(motif).unwrap_or("?");
    let repeat_count = region_len as f64 / motif.len() as f64;

    // Format: "31 x CAG" or "12.3 x CCG"
    let label = if (repeat_count - repeat_count.round()).abs() < 0.05 {
        format!("{} x {}", repeat_count.round() as u32, motif_str)
    } else {
        format!("{:.1} x {}", repeat_count, motif_str)
    };

    let seg = Seg {
        width: region_len as u32,
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
