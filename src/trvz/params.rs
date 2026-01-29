use super::align::SegType;
use std::{collections::HashMap, fmt};

pub type ColorMap = HashMap<SegType, Color>;

pub struct PlotParams {
    pub colors: ColorMap,
    pub pipe_height: u32,
    pub pipe_pad: u32,
}

pub fn pick_params(motifs: &[Vec<u8>], is_squished: bool) -> PlotParams {
    let colors = pick_colors(motifs);

    if !is_squished {
        PlotParams {
            colors,
            pipe_height: 4,
            pipe_pad: 1,
        }
    } else {
        PlotParams {
            colors,
            pipe_height: 1,
            pipe_pad: 0,
        }
    }
}

fn pick_colors(motifs: &[Vec<u8>]) -> ColorMap {
    let tr_colors = [
        Color::Blue,
        Color::Purple,
        Color::Orange,
        Color::Pink,
        Color::Yellow,
        Color::Green,
        Color::Red,
        Color::Khaki,
        Color::PaleRed,
        Color::PaleBlue,
    ];

    let mut colors = HashMap::new();
    colors.insert(SegType::LeftFlank, Color::Teal);
    colors.insert(SegType::RightFlank, Color::Teal);

    for index in 0..motifs.len() {
        let color = tr_colors[index % tr_colors.len()].clone();
        colors.insert(SegType::Tr(index), color);
    }
    colors.insert(SegType::Tr(motifs.len()), Color::LightGray);

    colors
}

pub fn get_meth_colors(motifs: &[Vec<u8>]) -> ColorMap {
    let mut colors = HashMap::new();
    colors.insert(SegType::LeftFlank, Color::Teal);
    colors.insert(SegType::RightFlank, Color::Teal);

    for index in 0..motifs.len() + 1 {
        colors.insert(SegType::Tr(index), Color::LightGray);
    }

    colors
}

#[derive(Debug, PartialEq, Clone)]
pub enum Color {
    Purple,
    Blue,
    Orange,
    Teal,
    Gray,
    LightGray,
    Black,
    Green,
    Pink,
    Yellow,
    Red,
    Khaki,
    PaleRed,
    PaleBlue,
    BaseA, // #FF6347 (tomato red)
    BaseT, // #FCA100 (orange/yellow)
    BaseG, // #2F8734 (green)
    BaseC, // #393939 (dark gray/charcoal)
    BaseN, // #000000 (black)
    Grad(f64),
}

impl fmt::Display for Color {
    fn fmt(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Color::Purple => write!(formatter, "#814ED1"),
            Color::Blue => write!(formatter, "#bccbe8"),
            Color::Orange => write!(formatter, "#E16A2C"),
            Color::Teal => write!(formatter, "#e0e0e0"),
            Color::Gray => write!(formatter, "#BFBFBF"),
            Color::LightGray => write!(formatter, "#D1D1D1"),
            Color::Black => write!(formatter, "#000000"),
            Color::Pink => write!(formatter, "#ED3981"),
            Color::Yellow => write!(formatter, "#EFCD17"),
            Color::Green => write!(formatter, "#009D4E"),
            Color::Red => write!(formatter, "#E3371E"),
            Color::Khaki => write!(formatter, "#F0E68C"),
            Color::PaleRed => write!(formatter, "#FF4858"),
            Color::PaleBlue => write!(formatter, "#46B2E8"),
            Color::BaseA => write!(formatter, "#22AA22"),
            Color::BaseT => write!(formatter, "#DD0000"),
            Color::BaseG => write!(formatter, "#CC9900"),
            Color::BaseC => write!(formatter, "#2222DD"),
            Color::BaseN => write!(formatter, "#000000"),
            Color::Grad(value) => write!(formatter, "{}", get_gradient(*value)),
        }
    }
}

fn get_gradient(value: f64) -> String {
    let blue: (u8, u8, u8) = (0, 73, 255);
    let red: (u8, u8, u8) = (255, 0, 0);
    let mix_red = (blue.0 as f64 * (1.0 - value) + red.0 as f64 * value).round() as u8;
    let mix_green = (blue.1 as f64 * (1.0 - value) + red.1 as f64 * value).round() as u8;
    let mix_blue = (blue.2 as f64 * (1.0 - value) + red.2 as f64 * value).round() as u8;

    format!("#{:02X}{:02X}{:02X}", mix_red, mix_green, mix_blue)
}

/// Returns the display color for a nucleotide base
pub fn base_color(base: u8) -> Color {
    match base {
        b'A' | b'a' => Color::BaseA,
        b'T' | b't' => Color::BaseT,
        b'G' | b'g' => Color::BaseG,
        b'C' | b'c' => Color::BaseC,
        b'N' | b'n' => Color::BaseN,
        _ => Color::Black,
    }
}
