#[derive(Debug, PartialEq)]
pub enum Shape {
    Rect,
    HLine,
    VLine(u32),  // Vertical line for insertions, with insertion size in bp
    None,
    Tick(Option<u32>),
    DoubleArrow(Option<String>),
}

pub type Color = String;

/// A text label to render on a pipe
#[derive(Debug, Clone)]
pub struct TextLabel {
    pub pos: u32,       // Position within pipe (in bp units from pipe start)
    pub text: char,     // Single character (A/T/G/C/N)
    pub color: Color,   // Color for this base
}

#[derive(Debug, PartialEq)]
pub struct Seg {
    pub width: u32,
    pub color: Color,
    pub shape: Shape,
    pub dashed: bool,   // if true, draw with dashed stroke
}

#[derive(Debug)]
pub struct Band {
    pub pos: u32, // Position relative to pipe's start
    pub width: u32,
    pub color: Color,
}

#[derive(Debug)]
pub struct Pipe {
    pub xpos: u32,
    pub ypos: u32,
    pub height: u32,
    pub segs: Vec<Seg>,
    pub bands: Vec<Band>,
    pub labels: Vec<TextLabel>,
    pub outline: bool,
}

#[derive(Debug)]
pub struct Legend {
    pub xpos: u32,
    pub ypos: u32,
    pub height: u32,
    pub labels: Vec<(String, String)>,
}

/// A vertical grid line spanning a region
#[derive(Debug)]
pub struct GridLine {
    pub xpos: u32,
    pub y_start: u32,
    pub y_end: u32,
    pub label_y: u32,  // Y position for the label (top of reads section)
    pub label: Option<String>,
}

#[derive(Debug)]
pub struct FontConfig {
    pub family: String,
    pub weight: String,
    pub size: String,
}

impl Default for FontConfig {
    fn default() -> Self {
        Self {
            family: "Arial Black, Helvetica Bold, sans-serif".to_string(),
            weight: "900".to_string(),
            size: "14px".to_string(),
        }
    }
}

#[derive(Debug)]
pub struct PipePlot {
    pub pipes: Vec<Pipe>,
    pub legend: Legend,
    pub font: FontConfig,
    pub grid_lines: Vec<GridLine>,
}

impl PipePlot {
    pub fn set_font_family(&mut self, font_family: &str) {
        self.font.family = font_family.to_owned();
    }
}
