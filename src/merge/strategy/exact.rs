use crate::merge::genotype_types::{remap_genotype, AlleleRegistry};
use crate::utils::Result;
use rust_htslib::bcf::record::GenotypeAllele;

#[allow(clippy::type_complexity)]
pub fn merge_exact<'a>(
    vcf_gts: Vec<Vec<Vec<GenotypeAllele>>>,
    sample_alleles: &'a [Vec<&'a [u8]>],
) -> Result<(Vec<Vec<Vec<GenotypeAllele>>>, Vec<&'a [u8]>)> {
    let registry = AlleleRegistry::build(sample_alleles)?;
    let merged_alleles = registry.merged_alleles();

    let out_sample_gts = vcf_gts
        .into_iter()
        .zip(sample_alleles.iter())
        .map(|(file_gts, file_alleles)| {
            let mapping = registry.create_mapping(file_alleles);
            let mut remapped_file_gts = Vec::with_capacity(file_gts.len());
            for mut sample_gt in file_gts {
                remap_genotype(&mut sample_gt, &mapping)?;
                remapped_file_gts.push(sample_gt);
            }
            Ok(remapped_file_gts)
        })
        .collect::<Result<Vec<_>>>()?;

    Ok((out_sample_gts, merged_alleles))
}

#[cfg(test)]
mod tests {
    use super::*;

    const REF_SEQ: &[u8] =
        b"CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA";
    const A: &[u8] = b"CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA";
    const B: &[u8] = b"CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA";
    const C: &[u8] = b"CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA";
    const D: &[u8] = b"CGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGCGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTG";
    const E: &[u8] = b"CGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGCGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTG";

    #[test]
    fn test_merge_exact() {
        let ref_seq = REF_SEQ;
        let a = A;
        let b = B;
        let c = C;
        let d = D;
        let e = E;

        let sample_gts = vec![
            vec![vec![
                GenotypeAllele::Unphased(1),
                GenotypeAllele::Unphased(2),
            ]],
            vec![vec![
                GenotypeAllele::Unphased(1),
                GenotypeAllele::Unphased(2),
            ]],
            vec![vec![
                GenotypeAllele::Unphased(0),
                GenotypeAllele::Unphased(0),
            ]],
            vec![vec![
                GenotypeAllele::UnphasedMissing,
                GenotypeAllele::UnphasedMissing,
            ]],
            vec![vec![
                GenotypeAllele::Unphased(1),
                GenotypeAllele::Unphased(2),
            ]],
        ];

        let sample_alleles = vec![
            vec![ref_seq, a, b],
            vec![ref_seq, a, c],
            vec![ref_seq],
            vec![ref_seq],
            vec![ref_seq, d, e],
        ];

        let (out_gts, sorted_alleles) = merge_exact(sample_gts, &sample_alleles).unwrap();

        assert_eq!(sorted_alleles.len(), 6);
        assert_eq!(sorted_alleles[0], ref_seq);
        assert_eq!(sorted_alleles[1], a);
        assert_eq!(sorted_alleles[2], b);
        assert_eq!(sorted_alleles[3], c);
        assert_eq!(sorted_alleles[4], d);
        assert_eq!(sorted_alleles[5], e);

        assert_eq!(
            out_gts[0][0],
            [GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(2)]
        );
        assert_eq!(
            out_gts[1][0],
            [GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(3)]
        );
        assert_eq!(
            out_gts[2][0],
            [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)]
        );
        assert_eq!(
            out_gts[3][0],
            [
                GenotypeAllele::UnphasedMissing,
                GenotypeAllele::UnphasedMissing
            ]
        );
        assert_eq!(
            out_gts[4][0],
            [GenotypeAllele::Unphased(4), GenotypeAllele::Unphased(5)]
        );
    }

    #[test]
    fn test_merge_exact_phasing() {
        let ref_seq = REF_SEQ;
        let a = A;
        let b = B;
        let c = D;
        let d = E;

        let sample_gts = vec![
            vec![vec![
                GenotypeAllele::Unphased(1),
                GenotypeAllele::Unphased(2),
            ]],
            vec![vec![GenotypeAllele::Phased(1), GenotypeAllele::Phased(2)]],
            vec![vec![GenotypeAllele::Phased(1), GenotypeAllele::Phased(2)]],
        ];

        let sample_alleles = vec![
            vec![ref_seq, a, b],
            vec![ref_seq, b, c],
            vec![ref_seq, c, d],
        ];

        let (out_gts, sorted_alleles) = merge_exact(sample_gts, &sample_alleles).unwrap();

        assert_eq!(sorted_alleles[0], ref_seq);

        assert_eq!(sorted_alleles.len(), 5);

        assert_eq!(
            out_gts[0][0],
            [GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(2)]
        );
        assert_eq!(
            out_gts[1][0],
            [GenotypeAllele::Phased(2), GenotypeAllele::Phased(3)]
        );
        assert_eq!(
            out_gts[2][0],
            [GenotypeAllele::Phased(3), GenotypeAllele::Phased(4)]
        );
    }

    #[test]
    fn test_merge_phase_ordering() {
        let ref_seq = REF_SEQ;
        let a = A;
        let b = B;

        let sample_gts = vec![
            vec![vec![GenotypeAllele::Phased(2), GenotypeAllele::Phased(1)]],
            vec![vec![GenotypeAllele::Phased(1), GenotypeAllele::Phased(0)]],
        ];

        let sample_alleles = vec![vec![ref_seq, a, b], vec![ref_seq, b]];

        let (out_gts, sorted_alleles) = merge_exact(sample_gts, &sample_alleles).unwrap();
        assert_eq!(sorted_alleles.len(), 3);
        assert_eq!(sorted_alleles[0], ref_seq);
        assert_eq!(sorted_alleles[1], a);
        assert_eq!(sorted_alleles[2], b);

        assert_eq!(
            out_gts[0][0],
            [GenotypeAllele::Phased(2), GenotypeAllele::Phased(1)]
        );
        assert_eq!(
            out_gts[1][0],
            [GenotypeAllele::Phased(2), GenotypeAllele::Phased(0)]
        );
    }
}
