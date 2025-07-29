use super::{consensus, Gt, TrSize};
use crate::trgt::reads::{extract_mismatch_offsets, HiFiRead};
use crate::utils::{align, median, GenomicRegion};
use itertools::Itertools;
use std::{cmp::Ordering, collections::BTreeMap};

type Profile = Vec<Option<bool>>;

pub fn genotype(
    reads: &[HiFiRead],
    tr_seqs: &[&str],
    region: &GenomicRegion,
) -> Option<(Gt, Vec<String>, Vec<i32>)> {
    let (trs_by_allele, mut allele_assignment) = get_trs_with_hp(reads, tr_seqs)
        .or_else(|| get_trs_with_clustering(reads, tr_seqs, region))?;
    let mut gt = Gt::new();
    let mut alleles = Vec::new();

    for trs in trs_by_allele {
        let (backbone, frequency) = simple_consensus(&trs)?;
        const MIN_FREQ_TO_ALIGN: f64 = 0.5;
        let allele = if frequency < MIN_FREQ_TO_ALIGN {
            let aligns = align(&backbone, &trs);
            consensus::repair_consensus(&backbone, &trs, &aligns)
        } else {
            backbone.to_string()
        };

        let min_tr_len = trs.iter().map(|tr| tr.len()).min().unwrap();
        let max_tr_len = trs.iter().map(|tr| tr.len()).max().unwrap();

        let size = TrSize::new(allele.len(), (min_tr_len, max_tr_len));
        gt.push(size);
        alleles.push(allele);
    }

    // Smaller allele should always appear first
    if alleles[0].len() > alleles[1].len() {
        gt.swap(0, 1);
        alleles.swap(0, 1);
        allele_assignment = allele_assignment.into_iter().map(|a| (a + 1) % 2).collect();
    }

    Some((gt, alleles, allele_assignment))
}

fn get_trs_with_hp<'a>(
    reads: &[HiFiRead],
    tr_seqs: &[&'a str],
) -> Option<([Vec<&'a str>; 2], Vec<i32>)> {
    let mut allele_assignment = Vec::new();
    let mut trs_by_allele = [Vec::new(), Vec::new()];
    let mut assignment_tie_breaker: usize = 1;
    let mut num_unassigned = 0;
    for (read, tr_seq) in reads.iter().zip(tr_seqs.iter()) {
        match read.hp_tag {
            Some(1) => {
                allele_assignment.push(0_i32);
                trs_by_allele[0].push(*tr_seq);
            }
            Some(2) => {
                allele_assignment.push(1_i32);
                trs_by_allele[1].push(*tr_seq);
            }
            _ => {
                assignment_tie_breaker = (assignment_tie_breaker + 1) % 2;
                allele_assignment.push(assignment_tie_breaker as i32);
                trs_by_allele[assignment_tie_breaker].push(*tr_seq);
                num_unassigned += 1;
            }
        }
    }

    let prop_assigned = (reads.len() - num_unassigned) as f64 / reads.len() as f64;
    if !trs_by_allele[0].is_empty() && !trs_by_allele[1].is_empty() && prop_assigned >= 0.7 {
        Some((trs_by_allele, allele_assignment))
    } else {
        None
    }
}

fn get_trs_with_clustering<'a>(
    reads: &[HiFiRead],
    tr_seqs: &[&'a str],
    region: &GenomicRegion,
) -> Option<([Vec<&'a str>; 2], Vec<i32>)> {
    if tr_seqs.is_empty() {
        return None;
    }

    let analysis_region = get_analysis_region(reads, region);
    let snvs = call_snvs(reads, region, 0.20, analysis_region);
    let profiles = get_profiles(reads, &snvs, region);

    let candidate_gts = get_candidate_gts(&profiles);

    // if there is just one genotype, it must be homozygous
    if candidate_gts.len() <= 1 {
        return None;
    }

    let gt_logliks = candidate_gts
        .iter()
        .map(|gt| get_loglik(gt, &profiles))
        .collect_vec();

    let top_gt = &candidate_gts
        .iter()
        .zip(gt_logliks.iter())
        .max_by(|(_gt1, ll1), (_gt2, ll2)| ll1.partial_cmp(ll2).unwrap())
        .unwrap()
        .0;

    // No SNPs were called since both haplotypes are identical
    if top_gt.0 == top_gt.1 {
        return None;
    }

    let mut allele_assignment = Vec::new();
    let mut assignment_tie_breaker = 1;
    let mut trs_by_allele = [Vec::new(), Vec::new()];
    for (index, profile) in profiles.iter().enumerate() {
        let dist1 = get_dist(profile, &top_gt.0);
        let dist2 = get_dist(profile, &top_gt.1);

        match dist1.cmp(&dist2) {
            Ordering::Less => {
                allele_assignment.push(0);
                trs_by_allele[0].push(tr_seqs[index]);
            }
            Ordering::Greater => {
                allele_assignment.push(1);
                trs_by_allele[1].push(tr_seqs[index]);
            }
            Ordering::Equal => {
                assignment_tie_breaker = (assignment_tie_breaker + 1) % 2;
                allele_assignment.push(assignment_tie_breaker);
                trs_by_allele[0].push(tr_seqs[index]);
                trs_by_allele[1].push(tr_seqs[index]);
            }
        }
    }
    Some((trs_by_allele, allele_assignment))
}

fn get_dist(read: &[Option<bool>], allele: &[bool]) -> usize {
    read.iter()
        .zip(allele.iter())
        .map(|(p, h)| (p.as_ref() == Some(h)) as usize)
        .sum()
}

/// Determine consensus for the input sequences
///
/// Return the most frequent sequence and its relative frequency.
/// If multiple sequences have the same frequency,
/// return the one whose length is closest to the median.
///
fn simple_consensus(seqs: &[&str]) -> Option<(String, f64)> {
    let median_len = median(&seqs.iter().map(|s| s.len() as i32).collect_vec())? as usize;
    let mut seq_to_count = BTreeMap::new();
    for seq in seqs {
        *seq_to_count.entry(*seq).or_insert(0) += 1;
    }
    let top_group_size = *seq_to_count.values().max()?;
    let consensus = seq_to_count
        .into_iter()
        .filter(|(_, c)| *c == top_group_size)
        .map(|(s, _)| (s, s.len().abs_diff(median_len)))
        .min_by_key(|(_, delta)| *delta)
        .map(|(s, _)| s.to_string())?;
    let top_group_frequency = (top_group_size as f64) / (seqs.len() as f64);
    Some((consensus, top_group_frequency))
}

fn get_loglik(gt: &(Vec<bool>, Vec<bool>), profiles: &[Profile]) -> f64 {
    profiles
        .iter()
        .map(|profile| {
            let term1 = eval_profile_given_hap(profile, &gt.0);
            let term2 = eval_profile_given_hap(profile, &gt.1);
            ln_sum_exp(term1, term2) - std::f64::consts::LN_2
        })
        .sum()
}

fn eval_profile_given_hap(profile: &Profile, hap: &[bool]) -> f64 {
    const MATCH_PROB: f64 = 0.9;
    const MISMATCH_PROB: f64 = 1.0 - MATCH_PROB;
    let match_ln = MATCH_PROB.ln();
    let mismatch_ln = MISMATCH_PROB.ln();
    let diff_ln = match_ln - mismatch_ln;
    let (matches, total_compared) = profile
        .iter()
        .zip(hap.iter())
        .filter_map(|(p, h)| p.as_ref().map(|p_val| p_val == h))
        .fold((0, 0), |(matches, total), is_match| {
            (matches + is_match as u32, total + 1)
        });
    (matches as f64 * diff_ln) + (total_compared as f64 * mismatch_ln)
}

fn ln_sum_exp(term1: f64, term2: f64) -> f64 {
    let max_term = term1.max(term2);
    max_term + ((term1 - max_term).exp() + (term2 - max_term).exp()).ln()
}

fn get_candidate_gts(profiles: &[Profile]) -> Vec<(Vec<bool>, Vec<bool>)> {
    let haps = profiles
        .iter()
        .filter(|p| p.iter().all(Option::is_some))
        .collect_vec();

    const PUTATIVE_HAP_FRAC: f64 = 0.40;
    if (haps.len() as f64 / profiles.len() as f64) < PUTATIVE_HAP_FRAC {
        return Vec::new();
    }

    let mut haps: Vec<Vec<bool>> = haps
        .into_iter()
        .map(|p| p.iter().map(|v| v.unwrap()).collect())
        .collect();

    haps.sort();
    haps.dedup();

    haps.into_iter()
        .combinations_with_replacement(2)
        .map(|mut combo| (combo.remove(0), combo.remove(0)))
        .collect()
}

fn get_analysis_region(reads: &[HiFiRead], region: &GenomicRegion) -> (i64, i64) {
    // min fraction of reads that must fully cover the left (right) flank
    const COV_READ_FRAC: f64 = 0.85;
    let skip_count = (reads.len() as f64 * (1.0 - COV_READ_FRAC)).round() as usize;

    let start = reads
        .iter()
        .map(|r| (r.ref_start - region.start as i64))
        .sorted()
        .nth_back(skip_count)
        .unwrap();

    let end = reads
        .iter()
        .map(|r| (r.ref_end - region.end as i64))
        .sorted()
        .nth(skip_count)
        .unwrap();

    (region.start as i64 + start, region.end as i64 + end)
}

fn get_profiles(reads: &[HiFiRead], snvs: &[u32], region: &GenomicRegion) -> Vec<Profile> {
    let region_start = region.start as i64;
    let region_end = region.end as i64;
    let snv_offsets = extract_mismatch_offsets(snvs, region.start, region.end);
    debug_assert_eq!(snvs.len(), snv_offsets.len());

    reads
        .iter()
        .map(|read| {
            let read_mismatches: std::collections::HashSet<u32> = read
                .mismatch_positions
                .as_deref()
                .unwrap_or(&[])
                .iter()
                .cloned()
                .collect();

            let read_start_offset = read.ref_start - region_start;
            let read_end_offset = read.ref_end - region_end;

            snvs.iter()
                .zip(snv_offsets.iter())
                .map(|(&snv_pos, &snv_offset)| {
                    let snv_offset = snv_offset as i64;
                    if snv_offset < read_start_offset || snv_offset > read_end_offset {
                        None
                    } else {
                        Some(read_mismatches.contains(&snv_pos))
                    }
                })
                .collect()
        })
        .collect()
}

fn call_snvs(
    reads: &[HiFiRead],
    region: &GenomicRegion,
    min_freq: f64,
    analysis_region: (i64, i64),
) -> Vec<u32> {
    let (lb, rb) = analysis_region;
    let offset_counts = reads
        .iter()
        .filter_map(|r| r.mismatch_positions.as_ref())
        .flatten()
        .filter(|&&pos| {
            (pos >= lb as u32 && pos < region.start) || (pos >= region.end && pos <= rb as u32)
        })
        .counts();

    let total_reads = reads.len() as f64;
    offset_counts
        .into_iter()
        .filter(|(_, count)| *count as f64 / total_reads >= min_freq)
        .map(|(offset, _)| *offset)
        .sorted()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;

    struct TestRead<'a> {
        left_flank: &'a str,
        tr_seq: &'a str,
        right_flank: &'a str,
        ref_start: i64,
        hp_tag: Option<u8>,
    }

    fn make_read(spec: &TestRead) -> HiFiRead {
        let left_len = spec.left_flank.len();
        let tr_len = spec.tr_seq.len();

        let mut mismatch_positions = Vec::new();
        for (i, c) in spec.left_flank.chars().enumerate() {
            if c == 'X' {
                mismatch_positions.push((spec.ref_start as usize + i) as u32);
            }
        }
        for (i, c) in spec.right_flank.chars().enumerate() {
            if c == 'X' {
                mismatch_positions.push((spec.ref_start as usize + left_len + tr_len + i) as u32);
            }
        }

        let ref_end = spec.ref_start + (left_len + tr_len + spec.right_flank.len()) as i64;

        HiFiRead {
            id: "read".to_string(),
            is_reverse: false,
            bases: spec.tr_seq.as_bytes().to_vec(),
            quals: "(".repeat(tr_len).as_bytes().to_vec(),
            meth: None,
            read_qual: None,
            mismatch_positions: Some(mismatch_positions),
            cigar: None,
            hp_tag: spec.hp_tag,
            mapq: 60,
            ref_start: spec.ref_start,
            ref_end,
        }
    }

    fn make_reads(specs: &[TestRead]) -> Vec<HiFiRead> {
        specs.iter().map(make_read).collect_vec()
    }

    #[test]
    fn if_het_snvs_then_genotype() {
        let reads = make_reads(&[
            TestRead {
                left_flank: "XX====",
                tr_seq: "TATATATA",
                right_flank: "===X===",
                ref_start: 0,
                hp_tag: None,
            },
            TestRead {
                left_flank: "XX=X==",
                tr_seq: "TATATATA",
                right_flank: "===X===",
                ref_start: 0,
                hp_tag: None,
            },
            TestRead {
                left_flank: "XX===",
                tr_seq: "TATATATATA",
                right_flank: "X=X===",
                ref_start: 0,
                hp_tag: None,
            },
            TestRead {
                left_flank: "XX===",
                tr_seq: "TATATATATA",
                right_flank: "X=X===",
                ref_start: 0,
                hp_tag: None,
            },
            TestRead {
                left_flank: "XX===",
                tr_seq: "TATATATATA",
                right_flank: "X==",
                ref_start: 0,
                hp_tag: None,
            },
            TestRead {
                left_flank: "======",
                tr_seq: "TATATATA",
                right_flank: "===X===",
                ref_start: 0,
                hp_tag: None,
            },
        ]);
        let tr_seqs = reads
            .iter()
            .map(|r| std::str::from_utf8(&r.bases).unwrap())
            .collect_vec();
        let region = GenomicRegion::new("chr1", 6, 14).unwrap();
        let result = genotype(&reads, &tr_seqs, &region);

        let gt = Gt::from([TrSize::new(8, (8, 8)), TrSize::new(10, (10, 10))]);

        let alleles = vec!["TATATATA".to_string(), "TATATATATA".to_string()];
        let assignment = vec![0, 0, 1, 1, 1, 0];

        assert_eq!(result, Some((gt, alleles, assignment)));
    }

    #[test]
    fn if_hom_snvs_then_none() {
        let reads = make_reads(&[
            TestRead {
                left_flank: "XX====",
                tr_seq: "TATATATATA",
                right_flank: "=X=X===",
                ref_start: 0,
                hp_tag: None,
            },
            TestRead {
                left_flank: "XX====",
                tr_seq: "TATATATATA",
                right_flank: "=X=X===",
                ref_start: 0,
                hp_tag: None,
            },
            TestRead {
                left_flank: "XX====",
                tr_seq: "TATATATATA",
                right_flank: "=X=X===",
                ref_start: 0,
                hp_tag: None,
            },
            TestRead {
                left_flank: "XX====",
                tr_seq: "TATATATATA",
                right_flank: "=X=X===",
                ref_start: 0,
                hp_tag: None,
            },
        ]);

        let tr_seqs = reads
            .iter()
            .map(|r| std::str::from_utf8(&r.bases).unwrap())
            .collect_vec();
        let region = GenomicRegion::new("chr1", 6, 16).unwrap();
        let result = genotype(&reads, &tr_seqs, &region);
        assert_eq!(result, None);
    }

    #[test]
    fn if_no_snvs_then_none() {
        let reads = make_reads(&[
            TestRead {
                left_flank: "======",
                tr_seq: "TATATATA",
                right_flank: "=======",
                ref_start: 0,
                hp_tag: None,
            },
            TestRead {
                left_flank: "======",
                tr_seq: "TATATATATA",
                right_flank: "=======",
                ref_start: 0,
                hp_tag: None,
            },
        ]);
        let tr_seqs = reads
            .iter()
            .map(|r| std::str::from_utf8(&r.bases).unwrap())
            .collect_vec();
        let region = GenomicRegion::new("chr1", 6, 14).unwrap();
        let result = genotype(&reads, &tr_seqs, &region);
        assert_eq!(result, None);
    }

    #[test]
    fn if_hp_tags_then_genotype() {
        let reads = make_reads(&[
            TestRead {
                left_flank: "======",
                tr_seq: "TATATA",
                right_flank: "======",
                ref_start: 0,
                hp_tag: Some(1),
            },
            TestRead {
                left_flank: "======",
                tr_seq: "TATATA",
                right_flank: "======",
                ref_start: 0,
                hp_tag: Some(1),
            },
            TestRead {
                left_flank: "======",
                tr_seq: "TATATATATA",
                right_flank: "======",
                ref_start: 0,
                hp_tag: Some(2),
            },
            TestRead {
                left_flank: "======",
                tr_seq: "TATATATATA",
                right_flank: "======",
                ref_start: 0,
                hp_tag: Some(2),
            },
            TestRead {
                left_flank: "======",
                tr_seq: "TATATATATA",
                right_flank: "======",
                ref_start: 0,
                hp_tag: Some(2),
            },
        ]);
        let tr_seqs = reads
            .iter()
            .map(|r| std::str::from_utf8(&r.bases).unwrap())
            .collect_vec();
        let region = GenomicRegion::new("chr1", 6, 12).unwrap();
        let result = genotype(&reads, &tr_seqs, &region);

        let gt = Gt::from([TrSize::new(6, (6, 6)), TrSize::new(10, (10, 10))]);
        let alleles = vec!["TATATA".to_string(), "TATATATATA".to_string()];
        let assignment = vec![0, 0, 1, 1, 1];

        assert_eq!(result, Some((gt, alleles, assignment)));
    }

    #[test]
    fn if_hp_tags_insufficient_then_cluster() {
        let reads = make_reads(&[
            TestRead {
                left_flank: "XX====",
                tr_seq: "TATATATA",
                right_flank: "===X===",
                ref_start: 0,
                hp_tag: Some(1),
            },
            TestRead {
                left_flank: "XX=X==",
                tr_seq: "TATATATA",
                right_flank: "===X===",
                ref_start: 0,
                hp_tag: None,
            },
            TestRead {
                left_flank: "XX===",
                tr_seq: "TATATATATA",
                right_flank: "X=X===",
                ref_start: 0,
                hp_tag: Some(2),
            },
            TestRead {
                left_flank: "XX===",
                tr_seq: "TATATATATA",
                right_flank: "X=X===",
                ref_start: 0,
                hp_tag: None,
            },
            TestRead {
                left_flank: "XX===",
                tr_seq: "TATATATATA",
                right_flank: "X==",
                ref_start: 0,
                hp_tag: None,
            },
            TestRead {
                left_flank: "======",
                tr_seq: "TATATATA",
                right_flank: "===X===",
                ref_start: 0,
                hp_tag: None,
            },
        ]);
        let tr_seqs = reads
            .iter()
            .map(|r| std::str::from_utf8(&r.bases).unwrap())
            .collect_vec();
        let region = GenomicRegion::new("chr1", 6, 14).unwrap();
        let result = genotype(&reads, &tr_seqs, &region);

        let gt = Gt::from([TrSize::new(8, (8, 8)), TrSize::new(10, (10, 10))]);

        let alleles = vec!["TATATATA".to_string(), "TATATATATA".to_string()];
        let assignment = vec![0, 0, 1, 1, 1, 0];

        assert_eq!(result, Some((gt, alleles, assignment)));
    }

    #[test]
    fn snv_at_tr_start_is_ignored() {
        let region = GenomicRegion::new("chr1", 100, 110).unwrap();
        let reads = vec![HiFiRead {
            id: "read1".to_string(),
            is_reverse: false,
            bases: "GATTACA".as_bytes().to_vec(),
            quals: "!!!!!!!".as_bytes().to_vec(),
            meth: None,
            read_qual: None,
            mismatch_positions: Some(vec![
                99,  // outside (include)
                100, // boundary (exclude)
                101, // inside (exclude)
            ]),
            cigar: None,
            hp_tag: None,
            mapq: 60,
            ref_start: 90,
            ref_end: 120,
        }];

        let analysis_region = get_analysis_region(&reads, &region);
        let snvs = call_snvs(&reads, &region, 0.20, analysis_region);
        assert_eq!(snvs, vec![99]);
    }

    #[test]
    fn snv_at_tr_end_is_used() {
        let region = GenomicRegion::new("chr1", 100, 110).unwrap();
        let reads = vec![HiFiRead {
            id: "read1".to_string(),
            is_reverse: false,
            bases: "GATTACA".as_bytes().to_vec(),
            quals: "!!!!!!!".as_bytes().to_vec(),
            meth: None,
            read_qual: None,
            //
            mismatch_positions: Some(vec![
                109, // inside (exclude)
                110, // boundary (include)
                111, // outside (include)
            ]),
            cigar: None,
            hp_tag: None,
            mapq: 60,
            ref_start: 90,
            ref_end: 120,
        }];

        let analysis_region = get_analysis_region(&reads, &region);
        let snvs = call_snvs(&reads, &region, 0.20, analysis_region);
        assert_eq!(snvs, vec![110, 111]);
    }

    #[test]
    fn test_get_profiles_basic() {
        let region = GenomicRegion::new("chr1", 100, 110).unwrap();
        let snvs = vec![95, 115];

        let reads = vec![
            HiFiRead {
                id: "read1".to_string(),
                is_reverse: false,
                bases: "GATTACA".as_bytes().to_vec(),
                quals: "!!!!!!!".as_bytes().to_vec(),
                meth: None,
                read_qual: None,
                mismatch_positions: Some(vec![95, 115]),
                cigar: None,
                hp_tag: None,
                mapq: 60,
                ref_start: 90,
                ref_end: 120,
            },
            HiFiRead {
                id: "read2".to_string(),
                is_reverse: false,
                bases: "GATTACA".as_bytes().to_vec(),
                quals: "!!!!!!!".as_bytes().to_vec(),
                meth: None,
                read_qual: None,
                mismatch_positions: Some(vec![95]), // Has only first mismatch
                cigar: None,
                hp_tag: None,
                mapq: 60,
                ref_start: 90,
                ref_end: 120,
            },
        ];

        let profiles = get_profiles(&reads, &snvs, &region);

        assert_eq!(profiles.len(), 2);
        assert_eq!(profiles[0], vec![Some(true), Some(true)]);
        assert_eq!(profiles[1], vec![Some(true), Some(false)]);
    }

    #[test]
    fn test_get_profiles_no_coverage() {
        let region = GenomicRegion::new("chr1", 100, 110).unwrap();
        let snvs = vec![95, 115];

        let reads = vec![HiFiRead {
            id: "read1".to_string(),
            is_reverse: false,
            bases: "GATTACA".as_bytes().to_vec(),
            quals: "!!!!!!!".as_bytes().to_vec(),
            meth: None,
            read_qual: None,
            mismatch_positions: Some(vec![95]),
            cigar: None,
            hp_tag: None,
            mapq: 60,
            ref_start: 96,
            ref_end: 114,
        }];

        let profiles = get_profiles(&reads, &snvs, &region);
        assert_eq!(profiles.len(), 1);
        assert_eq!(profiles[0], vec![None, None]);
    }

    #[test]
    fn test_get_profiles_partial_coverage() {
        let region = GenomicRegion::new("chr1", 100, 110).unwrap();
        let snvs = vec![95, 115];

        let reads = vec![HiFiRead {
            id: "read1".to_string(),
            is_reverse: false,
            bases: "GATTACA".as_bytes().to_vec(),
            quals: "!!!!!!!".as_bytes().to_vec(),
            meth: None,
            read_qual: None,
            mismatch_positions: Some(vec![95]),
            cigar: None,
            hp_tag: None,
            mapq: 60,
            ref_start: 90,
            ref_end: 114, // Does not cover second SNV
        }];

        let profiles = get_profiles(&reads, &snvs, &region);
        assert_eq!(profiles.len(), 1);
        assert_eq!(profiles[0], vec![Some(true), None]);
    }

    #[test]
    fn test_get_profiles_no_mismatches() {
        let region = GenomicRegion::new("chr1", 100, 110).unwrap();
        let snvs = vec![95, 115];

        let reads = vec![HiFiRead {
            id: "read1".to_string(),
            is_reverse: false,
            bases: "GATTACA".as_bytes().to_vec(),
            quals: "!!!!!!!".as_bytes().to_vec(),
            meth: None,
            read_qual: None,
            mismatch_positions: None,
            cigar: None,
            hp_tag: None,
            mapq: 60,
            ref_start: 90,
            ref_end: 120,
        }];

        let profiles = get_profiles(&reads, &snvs, &region);
        assert_eq!(profiles.len(), 1);
        assert_eq!(profiles[0], vec![Some(false), Some(false)]);
    }
}
