use crate::wfaligner::CigarOp;
use arrayvec::ArrayVec;
use std::mem;

pub fn repair_consensus(reference: &[u8], seqs: &[&[u8]], aligns: &[Vec<CigarOp>]) -> Vec<u8> {
    //                                        A  T  C  G  -
    let mut ref_counts: Vec<[u32; 5]> = vec![[0, 0, 0, 0, 0]; reference.len()];
    let mut ref_inserts: Vec<Vec<&[u8]>> = vec![Vec::new(); reference.len() + 1];
    for (seq, ops) in seqs.iter().zip(aligns) {
        let mut x = 0;
        let mut y = 0;
        for &(len, op) in ops {
            match op {
                '=' | 'M' | 'X' => {
                    let piece = &seq[x..x + len];
                    for (offset, base) in piece.iter().enumerate() {
                        let index = base_to_index(*base);
                        ref_counts[y + offset][index] += 1;
                    }
                    x += len;
                    y += len;
                }
                'D' => {
                    for off in 0..len {
                        ref_counts[y + off][4] += 1;
                    }
                    y += len;
                }
                'I' => {
                    let piece = &seq[x..x + len];
                    ref_inserts[y].push(piece);
                    x += len;
                }
                'N' | 'S' | 'H' | 'P' => {
                    panic!("Unexpected CIGAR operation: {}", op);
                }
                _ => panic!("Unknown CIGAR operation: {}", op),
            };
        }
    }

    let consensus_indexes: Vec<usize> = ref_counts
        .into_iter()
        .map(|rec| {
            rec.iter()
                .enumerate()
                .max_by_key(|(_, val)| *val)
                .unwrap()
                .0
        })
        .collect();

    let mut consensus = Vec::with_capacity(reference.len() + reference.len() / 8);
    for ref_pos in 0..=reference.len() {
        {
            let mut bucket = mem::take(&mut ref_inserts[ref_pos]);
            append_ins_consensus(&mut consensus, &mut bucket, seqs.len());
            ref_inserts[ref_pos] = bucket;
        }
        if ref_pos < reference.len() {
            let base_index = consensus_indexes[ref_pos];
            if base_index != 4 {
                consensus.push(index_to_base(base_index));
            }
        }
    }
    consensus
}

#[inline]
fn base_to_index(base: u8) -> usize {
    match base {
        b'A' => 0,
        b'T' => 1,
        b'C' => 2,
        b'G' => 3,
        _ => panic!("Encountered unexpected base"),
    }
}
#[inline]
fn index_to_base(index: usize) -> u8 {
    match index {
        0 => b'A',
        1 => b'T',
        2 => b'C',
        3 => b'G',
        _ => panic!("Encountered unexpected base index"),
    }
}
fn append_ins_consensus(consensus: &mut Vec<u8>, bucket: &mut Vec<&[u8]>, n_reads: usize) {
    debug_assert!(n_reads >= bucket.len());

    if bucket.is_empty() {
        return;
    }
    let reads_without = n_reads - bucket.len();
    if bucket.len() <= reads_without {
        return;
    }

    bucket.sort();

    let mut best_key = bucket[0];
    let mut best_cnt = 1usize;
    let mut cur_key = bucket[0];
    let mut cur_cnt = 1usize;

    for &x in &bucket[1..] {
        if x == cur_key {
            cur_cnt += 1;
        } else {
            if cur_cnt > best_cnt {
                best_cnt = cur_cnt;
                best_key = cur_key;
            }
            cur_key = x;
            cur_cnt = 1;
        }
    }
    if cur_cnt > best_cnt {
        best_cnt = cur_cnt;
        best_key = cur_key;
    }

    if best_cnt > reads_without {
        consensus.extend_from_slice(best_key);
    }
}

pub fn get_consensus(sizes: &ArrayVec<usize, 2>, seqs: &[&[u8]], counts: &[usize]) -> Vec<Vec<u8>> {
    let mut consensuses = Vec::new();

    let allele = get_closest_size(seqs, sizes[0]).unwrap();
    let consensus = get_most_frequent_seq(seqs, counts, allele);
    consensuses.push(consensus);

    if sizes.len() != 1 && sizes[0] != sizes[1] {
        let allele = get_closest_size(seqs, sizes[1]).unwrap();
        let consensus = get_most_frequent_seq(seqs, counts, allele);
        consensuses.push(consensus);
    }

    consensuses
}

fn get_closest_size(seqs: &[&[u8]], allele: usize) -> Option<usize> {
    let mut closest_size = None;
    for seq in seqs {
        let read_len = seq.len();

        if closest_size.is_none() {
            closest_size = Some(read_len);
            continue;
        }

        if closest_size.unwrap().abs_diff(allele) > read_len.abs_diff(allele) {
            closest_size = Some(read_len);
        }
    }

    closest_size
}

fn get_most_frequent_seq(seqs: &[&[u8]], counts: &[usize], length: usize) -> Vec<u8> {
    seqs.iter()
        .zip(counts)
        .filter(|rec| rec.0.len() == length)
        .max_by_key(|rec| rec.1)
        .unwrap()
        .0
        .to_vec()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_repair(reference: &str, reads: &[&str], cigars: &[Vec<CigarOp>]) -> String {
        let store: Vec<Vec<u8>> = reads.iter().map(|s| s.as_bytes().to_vec()).collect();
        let seqs: Vec<&[u8]> = store.iter().map(|v| v.as_slice()).collect();
        let consensus = repair_consensus(reference.as_bytes(), &seqs, cigars);
        String::from_utf8(consensus).unwrap()
    }

    #[test]
    fn test_trailing_insertion_is_emitted() {
        let cig = vec![vec![(2, 'M'), (10, 'I')], vec![(2, 'M'), (10, 'I')]];
        let consensus = run_repair("AC", &["ACGGGGGGGGGG", "ACGGGGGGGGGG"], &cig);
        assert_eq!(consensus, "ACGGGGGGGGGG");
    }

    #[test]
    fn test_emits_trailing_insertion_with_majority() {
        let reference = b"ACG";
        let seqs: Vec<&[u8]> = vec![b"ACGTT", b"ACGTT"];
        let aligns: Vec<Vec<CigarOp>> = vec![vec![(3, '='), (2, 'I')], vec![(3, '='), (2, 'I')]];
        let out = repair_consensus(reference, &seqs, &aligns);
        assert_eq!(out, b"ACGTT");
    }

    #[test]
    fn test_suppresses_trailing_insertion_without_majority() {
        let reference = b"ACG";
        let seqs: Vec<&[u8]> = vec![b"ACGTT", b"ACG"];
        let aligns: Vec<Vec<CigarOp>> = vec![vec![(3, '='), (2, 'I')], vec![(3, '=')]];
        let out = repair_consensus(reference, &seqs, &aligns);
        assert_eq!(out, b"ACG");
    }

    #[test]
    fn test_insertion_between_bases_is_emitted() {
        let cig = vec![
            vec![(1, 'M'), (2, 'I'), (1, 'M')],
            vec![(1, 'M'), (2, 'I'), (1, 'M')],
        ];
        let got = run_repair("AC", &["AGGC", "AGGC"], &cig);
        assert_eq!(got, "AGGC");
    }

    #[test]
    fn test_majority_deletion_wins() {
        let cig = vec![
            vec![(1, 'M'), (1, 'D'), (2, 'M')],
            vec![(1, 'M'), (1, 'D'), (2, 'M')],
            vec![(1, 'M'), (1, 'D'), (2, 'M')],
            vec![(4, 'M')],
        ];
        let got = run_repair("ACGT", &["AGT", "AGT", "AGT", "ACGT"], &cig);
        assert_eq!(got, "AGT");
    }

    #[test]
    fn test_mismatches_in_consensus() {
        let cig = vec![
            vec![(1, 'X')],
            vec![(1, 'X')],
            vec![(1, 'X')],
            vec![(1, 'M')],
        ];
        let got = run_repair("A", &["G", "G", "G", "A"], &cig);
        assert_eq!(got, "G");
    }

    #[test]
    fn test_insertion_majority_threshold() {
        let cig1 = vec![vec![(1, 'I'), (1, 'M')], vec![(1, 'M')]];
        let got1 = run_repair("A", &["GA", "A"], &cig1);
        assert_eq!(got1, "A");

        let cig2 = vec![vec![(1, 'I'), (1, 'M')], vec![(1, 'I'), (1, 'M')]];
        let got2 = run_repair("A", &["GA", "GA"], &cig2);
        assert_eq!(got2, "GA");
    }
}
