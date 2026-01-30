use crate::pipeplot::{Color, FontConfig, GridLine, Legend, Pipe, PipePlot, Shape};
use std::{fs, path::Path};

const DEFAULT_X_SCALE: f64 = 2500.0;
const DEFAULT_Y_SCALE: f64 = 4.5;
const DEFAULT_PADDING: f64 = 12.0;

pub fn generate_string(plot: &PipePlot) -> String {
    let longest_pipe = get_longest_pipe(plot);
    let x_scale = DEFAULT_X_SCALE / longest_pipe as f64;
    let scale = (x_scale, DEFAULT_Y_SCALE);
    let mut generator = Generator::new(scale, DEFAULT_PADDING);
    generator.generate(plot);
    generator.buffer
}

pub fn render_from_string(svg_content: &str, path: &Path) -> Result<(), String> {
    fs::write(path, svg_content).map_err(|e| e.to_string())
}

fn get_longest_pipe(plot: &PipePlot) -> u32 {
    plot.pipes
        .iter()
        .map(|pipe| pipe.segs.iter().map(|seg| seg.width).sum())
        .max()
        .unwrap_or(0)
}

struct Generator {
    scale: (f64, f64),
    pad: f64,
    buffer: String,
    force_rectangle_labels: bool,
}

impl Generator {
    fn new(scale: (f64, f64), pad: f64) -> Self {
        Self {
            scale,
            pad,
            buffer: String::with_capacity(10_000),
            force_rectangle_labels: false,
        }
    }

    /// Check if letters would overlap in any consensus pipe (outline: true)
    fn would_labels_overlap(&self, pipe_plot: &PipePlot) -> bool {
        // Find any consensus pipe (identified by outline: true)
        for pipe in &pipe_plot.pipes {
            if pipe.outline && !pipe.labels.is_empty() {
                // Calculate the width per base position in pixels
                let width_per_base = self.to_x(1);

                // Calculate the font size that would be used for this pipe
                let pipe_height_pixels = self.to_y(pipe.height);
                let font_size = (pipe_height_pixels * 0.6).max(6.0).min(14.0);

                // For monospace fonts, character width is roughly 60% of font height
                let char_width = font_size * 0.6;

                // If width per base is less than character width, letters would overlap
                if width_per_base < char_width {
                    return true;
                }
            }
        }
        false
    }

    pub fn generate(&mut self, pipe_plot: &PipePlot) {
        // Detect if letters would overlap in any consensus pipe (outline: true)
        // If so, switch all pipes to rectangle mode for labels
        self.force_rectangle_labels = self.would_labels_overlap(pipe_plot);

        let (width, height) = self.get_dimensions(pipe_plot);
        self.start_svg(width, height);
        self.add_background();

        // Draw grid lines first (behind pipes)
        for grid_line in &pipe_plot.grid_lines {
            self.plot_grid_line(grid_line, &pipe_plot.font);
        }

        for pipe in &pipe_plot.pipes {
            self.plot_pipe(pipe, &pipe_plot.font);
            if pipe.outline {
                self.plot_outline(pipe);
            }
        }
        self.plot_legend(&pipe_plot.legend, &pipe_plot.font);
        self.end_svg();
    }

    fn add_line(&mut self, line: &str) {
        self.buffer.reserve(line.len() + 1);
        self.buffer.push_str(line);
        self.buffer.push('\n');
    }

    fn plot_legend(&mut self, legend: &Legend, font: &FontConfig) {
        let base_x = self.to_x(legend.xpos) + self.pad;
        let base_y = self.to_y(legend.ypos) + self.pad;
        let height = self.to_y(legend.height);
        let mut x = base_x;
        for (label, color) in &legend.labels {
            self.add_rect((x, base_y), (height, height), color, false, false);
            x += height + 2.0;
            let point = format!("x=\"{}\" y=\"{}\"", x, base_y + height - 1.0);
            let font_style = format!(
                r#"font-family="{}" font-weight="{}" font-size="{}""#,
                font.family, font.weight, font.size
            );
            let line = format!("<text {} {} >{}</text>", point, font_style, label);
            self.add_line(&line);
            x += 5.0 * (2 * label.len() as u32 + 1) as f64;
        }
    }

    fn plot_grid_line(&mut self, grid_line: &GridLine, font: &FontConfig) {
        let x = self.to_x(grid_line.xpos) + self.pad;
        let y1 = self.to_y(grid_line.y_start) + self.pad;
        let y2 = self.to_y(grid_line.y_end) + self.pad;
        let label_y = self.to_y(grid_line.label_y) + self.pad;

        // Draw dashed vertical line
        let line_elem = format!(
            r##"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#999999" stroke-width="1" stroke-dasharray="3,3" />"##,
            x, y1, x, y2
        );
        self.add_line(&line_elem);

        // Draw label if present (just above consensus, with 2px gap)
        if let Some(label) = &grid_line.label {
            let text_elem = format!(
                r##"<text x="{}" y="{}" font-family="{}" font-weight="{}" font-size="11px" text-anchor="middle" fill="#666666">{}</text>"##,
                x, label_y - 5.0, font.family, font.weight, label
            );
            self.add_line(&text_elem);
        }
    }

    fn plot_pipe(&mut self, pipe: &Pipe, font: &FontConfig) {
        let x = self.to_x(pipe.xpos) + self.pad;
        let y = self.to_y(pipe.ypos) + self.pad;
        let add_highlight = pipe.height > 1;
        let pipe_height = self.to_y(pipe.height);
        let stroke = if pipe.height > 1 { 1.5 } else { 1.0 };
        // Plot the main pipe
        let mut x_cur = x;

        for seg in &pipe.segs {
            let dims = (self.to_x(seg.width), pipe_height);
            match &seg.shape {
                Shape::Rect => self.add_rect((x_cur, y), dims, &seg.color, add_highlight, seg.dashed),
                Shape::HLine => self.add_hline((x_cur, y), dims, &seg.color, stroke),
                Shape::Tick(label) => self.add_tick((x_cur, y), dims, &seg.color, *label, font),
                Shape::None | Shape::VLine(_) => {}
                Shape::DoubleArrow(label) => {
                    self.add_double_arrow((x_cur, y), dims, &seg.color, stroke, label, font)
                }
            }

            x_cur += self.to_x(seg.width);
        }

        // Plot insertions
        let mut x_cur = x;
        for seg in &pipe.segs {
            let dims = (self.to_x(seg.width), pipe_height);

            if let Shape::VLine(insertion_size) = seg.shape {
                self.add_vline((x_cur, y), dims, &seg.color, insertion_size, font);
            }

            x_cur += self.to_x(seg.width);
        }

        // Plot bands
        for band in &pipe.bands {
            let beta_x = x + self.to_x(band.pos);
            let dims = (self.to_x(1), pipe_height);
            self.add_rect((beta_x, y), dims, &band.color, false, false);
        }

        // Plot labels (rendered last so they appear on top)
        self.plot_labels(pipe, font);
    }

    fn plot_labels(&mut self, pipe: &Pipe, font: &FontConfig) {
        if pipe.labels.is_empty() {
            return;
        }

        let pipe_height_pixels = self.to_y(pipe.height);

        // Use rectangles if:
        // 1. Squished mode (pipe.height <= 1), or
        // 2. Consensus pipe (outline: true) and letters would overlap
        let use_rectangles = pipe.height <= 1 || (pipe.outline && self.force_rectangle_labels);
        if use_rectangles {
            for label in &pipe.labels {
                let x = self.to_x(pipe.xpos + label.pos) + self.pad;
                let y = self.to_y(pipe.ypos) + self.pad;
                let width = self.to_x(1);
                // Draw a small colored rectangle at the label position
                let rect = format!(
                    r#"<rect x="{}" y="{}" width="{}" height="{}" fill="{}" />"#,
                    x, y, width, pipe_height_pixels, label.color
                );
                self.add_line(&rect);
            }
            return;
        }

        // Normal mode: render text labels
        // Scale font size relative to pipe height (~60% of height, clamped 6-14px)
        let font_size = (pipe_height_pixels * 0.6).max(6.0).min(14.0);

        for label in &pipe.labels {
            let x = self.to_x(pipe.xpos + label.pos) + self.pad + self.to_x(1) / 2.0;
            let y = self.to_y(pipe.ypos) + self.pad + pipe_height_pixels / 2.0;

            // Use fill (not stroke) for text color, use font config
            let text_elem = format!(
                r#"<text x="{}" y="{}" dy="0.25em" text-anchor="middle" font-family="{}" font-weight="{}" font-size="{}px" fill="{}">{}</text>"#,
                x, y, font.family, font.weight, font_size, label.color, label.text
            );
            self.add_line(&text_elem);
        }
    }

    fn plot_outline(&mut self, pipe: &Pipe) {
        let height = self.to_y(pipe.height);
        let width = self.to_x(pipe.segs.iter().map(|seg| seg.width).sum());

        let x = self.to_x(pipe.xpos) + self.pad;
        let y = self.to_y(pipe.ypos) + self.pad;

        let dimensions = format!("width=\"{width}\" height=\"{height}\"");
        let pos = format!("x=\"{x}\" y=\"{y}\"");
        let style = r##"stroke="#000000" stroke-width="1.5" fill="transparent""##;
        let line = format!("<rect {} {} {} />", dimensions, pos, style);
        self.add_line(&line);
    }

    fn add_rect(&mut self, pos: (f64, f64), dims: (f64, f64), color: &Color, add_highlight: bool, dashed: bool) {
        let (x, y) = pos;
        let (w, h) = dims;

        let pos = format!("x=\"{}\" y=\"{}\"", x, y);
        let dim = format!("height=\"{}\" width=\"{}\"", h, w);

        let stroke_style = if dashed {
            r##"stroke="#000000" stroke-width="1" stroke-dasharray="2,2""##.to_string()
        } else {
            format!(r#"stroke="{}" stroke-width="0""#, color)
        };

        let style = format!("fill=\"{}\" {}", color, stroke_style);

        let rect = format!("<rect {} {} {} opacity=\"0.9\" />", pos, dim, style);
        self.add_line(&rect);

        if add_highlight {
            let pos = format!("x=\"{}\" y=\"{}\"", x, y + h * 0.18);
            let dim = format!("height=\"{}\" width=\"{}\"", h / 3.0, w);
            let highlight = format!("<rect {} {} fill=\"#F4EDF2\" opacity=\"0.25\" />", pos, dim);
            self.add_line(&highlight);
        }
    }

    fn add_hline(&mut self, pos: (f64, f64), dims: (f64, f64), color: &Color, stroke: f64) {
        let x1 = pos.0;
        let x2 = pos.0 + dims.0;
        let y1 = pos.1 + dims.1 / 2.0;
        let y2 = y1;

        let x1y1 = format!("x1=\"{}\" y1=\"{}\"", x1, y1);
        let x2y2 = format!("x2=\"{}\" y2=\"{}\"", x2, y2);

        let style = format!("stroke=\"{}\" stroke-width=\"{}\"", color, stroke);

        let line = format!("<line {} {} {} />", x1y1, x2y2, style);
        self.add_line(&line);
    }

    fn add_vline(&mut self, pos: (f64, f64), dims: (f64, f64), color: &Color, insertion_size: u32, font: &FontConfig) {
        let x1 = pos.0;
        let x2 = pos.0;
        let y1 = pos.1;
        let y2 = pos.1 + dims.1;

        let x1y1 = format!("x1=\"{}\" y1=\"{}\"", x1, y1);
        let x2y2 = format!("x2=\"{}\" y2=\"{}\"", x2, y2);

        let stroke_width = 2.0_f64.min(self.to_x(1));
        let style = format!("stroke=\"{}\" stroke-width=\"{}\"", color, stroke_width);

        let line = format!("<line {} {} {} />", x1y1, x2y2, style);
        self.add_line(&line);

        // Add label showing insertion size (to the right of the line, centered vertically)
        if insertion_size > 0 {
            let label = format!("{}", insertion_size);
            let label_x = x1 + 2.0;  // Just to the right of the line
            let label_y = (y1 + y2) / 2.0;  // Center vertically on the marker
            let font_size = (dims.1 * 0.5).max(6.0).min(10.0);
            let text_elem = format!(
                r#"<text x="{}" y="{}" font-family="{}" font-weight="{}" font-size="{}px" fill="{}" dominant-baseline="middle">{}</text>"#,
                label_x, label_y, font.family, font.weight, font_size, color, label
            );
            self.add_line(&text_elem);
        }
    }

    fn add_double_arrow(
        &mut self,
        pos: (f64, f64),
        dims: (f64, f64),
        _color: &Color,
        stroke: f64,
        label: &Option<String>,
        font: &FontConfig,
    ) {
        let x1 = pos.0;
        let x2 = pos.0 + dims.0;
        let y_center = pos.1 + dims.1 / 2.0;

        // Use dark gray color for lower contrast
        let arrow_color = "#888888";

        let x1y1 = format!("x1=\"{}\" y1=\"{}\"", x1, y_center);
        let x2y2 = format!("x2=\"{}\" y2=\"{}\"", x2, y_center);

        // Dashed line with long dashes (10px dash, 5px gap)
        let style = format!(
            r##"stroke="{}" stroke-width="{}" stroke-dasharray="10,5""##,
            arrow_color, stroke
        );

        // Draw the arrow line first (behind the label)
        let line = format!("<line {} {} {} />", x1y1, x2y2, style);
        self.add_line(&line);

        let arrow_pt1 = format!("{} {}", x1, y_center);
        let arrow_pt2 = format!("{} {}", x1 + 5.0, y_center + 5.0);
        let arrow_pt3 = format!("{} {}", x1 + 5.0, y_center - 5.0);
        let left_arrow = format!(
            r##"<polygon points="{arrow_pt1}, {arrow_pt2}, {arrow_pt3}" fill="{arrow_color}"/>"##
        );
        self.add_line(&left_arrow);

        let arrow_pt1 = format!("{} {}", x2, y_center);
        let arrow_pt2 = format!("{} {}", x2 - 5.0, y_center - 5.0);
        let arrow_pt3 = format!("{} {}", x2 - 5.0, y_center + 5.0);
        let right_arrow = format!(
            r##"<polygon points="{arrow_pt1}, {arrow_pt2}, {arrow_pt3}" fill="{arrow_color}"/>"##
        );
        self.add_line(&right_arrow);

        if let Some(label) = label {
            let label_x = (x1 + x2) / 2.0;
            // Estimate label width based on character count (rough approximation)
            let char_width = 10.0;
            let label_width = label.len() as f64 * char_width + 8.0;
            let label_height = 16.0;

            // Draw white background rectangle behind the label
            let bg_rect = format!(
                r##"<rect x="{}" y="{}" width="{}" height="{}" fill="white" />"##,
                label_x - label_width / 2.0,
                y_center - label_height / 2.0,
                label_width,
                label_height
            );
            self.add_line(&bg_rect);

            // Draw the label centered on the arrow line
            let font_style = format!(
                r#"font-family="{}" font-weight="{}" font-size="{}" text-anchor="middle" dominant-baseline="middle""#,
                font.family, font.weight, font.size
            );
            let text_elem = format!(
                r#"<text x="{}" y="{}" {} >{}</text>"#,
                label_x, y_center, font_style, label
            );
            self.add_line(&text_elem);
        }
    }

    fn add_tick(
        &mut self,
        pos: (f64, f64),
        dims: (f64, f64),
        color: &Color,
        label: Option<u32>,
        font: &FontConfig,
    ) {
        let x1 = pos.0;
        let x2 = pos.0;
        let y1 = pos.1;
        let y2 = pos.1 + dims.1;

        let x1y1 = format!("x1=\"{}\" y1=\"{}\"", x1, y1);
        let x2y2 = format!("x2=\"{}\" y2=\"{}\"", x2, y2);

        let stroke_width = 1.5;
        let style = format!("stroke=\"{}\" stroke-width=\"{}\"", color, stroke_width);

        let line = format!("<line {} {} {} />", x1y1, x2y2, style);
        self.add_line(&line);

        if let Some(label) = label {
            let point = format!("x=\"{}\" y=\"{}\"", x1, y1 - 2.0); // 2.0 is padding
            let font_style = format!(
                r#"font-family="{}" font-weight="{}" font-size="{}" text-anchor="middle""#,
                font.family, font.weight, font.size
            );
            let line = format!("<text {} {} >{}</text>", point, font_style, label);
            self.add_line(&line);
        }
    }

    fn get_dimensions(&self, pipe_plot: &PipePlot) -> (f64, f64) {
        let width = pipe_plot
            .pipes
            .iter()
            .map(|p| p.xpos + p.segs.iter().map(|s| s.width).sum::<u32>())
            .max()
            .unwrap_or(0);
        let height = pipe_plot.legend.ypos + pipe_plot.legend.height;

        let xdim = self.to_x(width) + 2.0 * self.pad;
        let ydim = self.to_y(height) + 2.0 * self.pad;
        (xdim, ydim)
    }

    fn start_svg(&mut self, width: f64, height: f64) {
        self.add_line(r#"<?xml version="1.0"?>"#);
        let line = format!(
            r#"<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="{}" height="{}">"#,
            width, height
        );
        self.add_line(&line);
    }

    fn end_svg(&mut self) {
        self.add_line("</svg>");
    }

    fn add_background(&mut self) {
        self.add_line(r#"<rect width="100%" height="100%" fill="white"/>"#);
    }

    fn to_x(&self, x: u32) -> f64 {
        x as f64 * self.scale.0
    }

    fn to_y(&self, y: u32) -> f64 {
        y as f64 * self.scale.1
    }
}
