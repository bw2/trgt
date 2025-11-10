use super::spans::Span;

pub fn count_motifs(motifs: &[Vec<u8>], labels: &Vec<Span>) -> Vec<usize> {
    let mut motif_counts = vec![0; motifs.len()];
    for span in labels {
        motif_counts[span.motif_index] += 1;
    }
    motif_counts
}

pub fn collapse_labels(spans: Vec<Span>) -> Vec<Span> {
    let mut collapsed = Vec::new();
    for span in spans {
        if collapsed.is_empty() {
            collapsed.push(span);
            continue;
        }

        let last_span = collapsed.last_mut().unwrap();
        if last_span.motif_index == span.motif_index && last_span.end == span.start {
            last_span.end = span.end;
        } else {
            collapsed.push(span);
        }
    }
    collapsed
}

#[inline]
pub fn replace_invalid_bases_inplace(seq: &mut [u8], allowed: &[u8]) {
    let mut is_allowed = [false; 256];
    for &b in allowed {
        is_allowed[b as usize] = true;
    }
    let n = allowed.len();
    for (i, b) in seq.iter_mut().enumerate() {
        if !is_allowed[*b as usize] {
            *b = allowed[i % n];
        }
    }
}

#[inline]
pub fn replace_invalid_bases(seq: &[u8], allowed: &[u8]) -> Vec<u8> {
    let mut v = seq.to_vec();
    replace_invalid_bases_inplace(&mut v, allowed);
    v
}
