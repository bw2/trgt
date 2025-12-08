#[derive(Debug, Clone, PartialEq)]
pub struct Span {
    pub motif_index: usize,
    pub start: usize,
    pub end: usize,
}

impl Span {
    pub fn len(&self) -> usize {
        self.end - self.start
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

pub type Spans = Vec<Span>;

#[derive(Debug, Clone)]
pub struct Annotation {
    pub labels: Option<Spans>,
    pub motif_counts: Vec<usize>,
    pub purity: f64,
}

impl Annotation {
    #[cfg(test)]
    pub fn base() -> Self {
        Self {
            labels: None,
            motif_counts: vec![0],
            purity: f64::NAN,
        }
    }
}
