use crate::utils::Result;
use rust_htslib::bcf::record::GenotypeAllele;
use std::collections::HashMap;

pub fn remap_genotype(g: &mut [GenotypeAllele], mapping: &[Option<i32>]) -> Result<()> {
    for allele in g {
        match allele {
            GenotypeAllele::PhasedMissing | GenotypeAllele::UnphasedMissing => {}
            GenotypeAllele::Phased(idx) | GenotypeAllele::Unphased(idx) => {
                *idx = lookup_idx(*idx, mapping)?;
            }
        }
    }
    Ok(())
}

#[inline]
fn lookup_idx(idx: i32, mapping: &[Option<i32>]) -> Result<i32> {
    let uidx = idx as usize;
    mapping
        .get(uidx)
        .and_then(|x| *x)
        .ok_or_else(|| format!("Allele index {} not found in mapping", idx))
}

#[derive(Debug)]
pub struct AlleleRegistry<'a> {
    ref_allele: &'a [u8],
    /// Sorted list of alternate alleles (by length, then lexicographically)
    alt_alleles: Vec<&'a [u8]>,
    /// Maps allele sequence to index in merged allele list (0 = ref, 1+ = alt)
    allele_to_index: HashMap<&'a [u8], usize>,
}

impl<'a> AlleleRegistry<'a> {
    pub fn build(allele_lists: &[Vec<&'a [u8]>]) -> Result<Self> {
        let ref_allele = Self::validate_ref_allele(allele_lists)?;

        let mut alt_alleles: Vec<&'a [u8]> = allele_lists
            .iter()
            .flat_map(|als| als.iter().skip(1))
            .copied()
            .collect();

        alt_alleles.sort_unstable_by(|a, b| {
            let len_cmp = a.len().cmp(&b.len());
            if len_cmp == std::cmp::Ordering::Equal {
                a.cmp(b)
            } else {
                len_cmp
            }
        });
        alt_alleles.dedup();

        let mut allele_to_index = HashMap::with_capacity(1 + alt_alleles.len());
        allele_to_index.insert(ref_allele, 0);
        for (i, allele) in alt_alleles.iter().enumerate() {
            allele_to_index.insert(*allele, i + 1);
        }
        Ok(AlleleRegistry {
            ref_allele,
            alt_alleles,
            allele_to_index,
        })
    }

    fn validate_ref_allele(allele_lists: &[Vec<&'a [u8]>]) -> Result<&'a [u8]> {
        let mut ref_allele: Option<&'a [u8]> = None;
        for als in allele_lists {
            if let Some(first) = als.first() {
                if let Some(existing) = ref_allele {
                    if existing != *first {
                        return Err(format!(
                            "Reference alleles do not match: '{}' vs '{}'",
                            String::from_utf8_lossy(existing),
                            String::from_utf8_lossy(first),
                        ));
                    }
                } else {
                    ref_allele = Some(*first);
                }
            }
        }
        ref_allele.ok_or_else(|| "No reference allele found".to_string())
    }

    pub fn merged_alleles(&self) -> Vec<&'a [u8]> {
        let mut result = Vec::with_capacity(1 + self.alt_alleles.len());
        result.push(self.ref_allele);
        result.extend(self.alt_alleles.iter().copied());
        result
    }

    pub fn create_mapping(&self, allele_list: &[&'a [u8]]) -> Vec<Option<i32>> {
        let mut mapping = vec![None; allele_list.len()];
        for (old_idx, allele) in allele_list.iter().enumerate() {
            if let Some(&new_idx) = self.allele_to_index.get(allele) {
                mapping[old_idx] = Some(new_idx as i32);
            }
        }
        mapping
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genotype_allele_remap() {
        // Mapping: 0->0, 1->2, 2->1
        let mapping = vec![Some(0), Some(2), Some(1)];

        let mut genotype = vec![GenotypeAllele::Unphased(1)];
        remap_genotype(&mut genotype, &mapping).unwrap();
        assert_eq!(genotype[0], GenotypeAllele::Unphased(2));

        let mut missing_genotype = vec![GenotypeAllele::UnphasedMissing];
        remap_genotype(&mut missing_genotype, &mapping).unwrap();
        assert_eq!(missing_genotype[0], GenotypeAllele::UnphasedMissing);

        let mut phased_genotype = vec![GenotypeAllele::Phased(2)];
        remap_genotype(&mut phased_genotype, &mapping).unwrap();
        assert_eq!(phased_genotype[0], GenotypeAllele::Phased(1));
    }

    #[test]
    fn test_allele_registry_basic() {
        let allele_lists = vec![
            vec![b"CAGCAG".as_ref(), b"CAG".as_ref(), b"CAGCAGCAG".as_ref()],
            vec![b"CAGCAG".as_ref(), b"CCG".as_ref()],
            vec![
                b"CAGCAG".as_ref(),
                b"CAG".as_ref(),
                b"CAGCAGCAGCAG".as_ref(),
            ],
        ];

        let registry = AlleleRegistry::build(&allele_lists).unwrap();
        let merged = registry.merged_alleles();

        // Ref allele should be first
        assert_eq!(merged[0], b"CAGCAG");
        // Then sorted by length, then lexicographically
        assert_eq!(merged[1], b"CAG");
        assert_eq!(merged[2], b"CCG");
        assert_eq!(merged[3], b"CAGCAGCAG");
        assert_eq!(merged[4], b"CAGCAGCAGCAG");
        assert_eq!(merged.len(), 5);
    }

    #[test]
    fn test_allele_registry_mapping() {
        let allele_lists = vec![
            vec![b"REF".as_ref(), b"ALT1".as_ref(), b"ALT2".as_ref()],
            vec![b"REF".as_ref(), b"ALT2".as_ref(), b"ALT3".as_ref()],
        ];

        let registry = AlleleRegistry::build(&allele_lists).unwrap();

        // First VCF: 0->0 (REF), 1->1 (ALT1), 2->2 (ALT2)
        let mapping1 = registry.create_mapping(&allele_lists[0]);
        assert_eq!(mapping1[0], Some(0));
        assert_eq!(mapping1[1], Some(1));
        assert_eq!(mapping1[2], Some(2));

        // Second VCF: 0->0 (REF), 1->2 (ALT2), 2->3 (ALT3)
        let mapping2 = registry.create_mapping(&allele_lists[1]);
        assert_eq!(mapping2[0], Some(0));
        assert_eq!(mapping2[1], Some(2));
        assert_eq!(mapping2[2], Some(3));
    }

    #[test]
    fn test_allele_registry_ref_mismatch() {
        let allele_lists = vec![
            vec![b"REF1".as_ref(), b"ALT1".as_ref()],
            vec![b"REF2".as_ref(), b"ALT2".as_ref()],
        ];

        let result = AlleleRegistry::build(&allele_lists);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("do not match"));
    }

    #[test]
    fn test_allele_registry_empty_list() {
        let allele_lists = vec![
            vec![b"REF".as_ref(), b"ALT1".as_ref()],
            vec![], // Missing record
        ];

        let registry = AlleleRegistry::build(&allele_lists).unwrap();
        let merged = registry.merged_alleles();
        assert_eq!(merged[0], b"REF");
        assert_eq!(merged.len(), 2); // REF + ALT1
    }
}
