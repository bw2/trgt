use super::field_registry::{
    FieldCategory, FieldDescriptor, FieldNumber, FieldValue, VcfType, MISSING_FLOAT,
    MISSING_INTEGER, VECTOR_END_FLOAT,
};
use crate::utils::Result;
use rust_htslib::bcf::Record;

pub static INFO_FIELDS: &[FieldDescriptor<'static>] = &[
    FieldDescriptor::builtin(
        b"TRID",
        FieldCategory::Info,
        VcfType::String,
        FieldNumber::Fixed(1),
        "Tandem repeat ID",
    )
    .required(),
    FieldDescriptor::builtin(
        b"END",
        FieldCategory::Info,
        VcfType::Integer,
        FieldNumber::Fixed(1),
        "End position of the variant described in this record",
    )
    .required(),
    FieldDescriptor::builtin(
        b"MOTIFS",
        FieldCategory::Info,
        VcfType::String,
        FieldNumber::Variable,
        "Motifs that the tandem repeat is composed of",
    )
    .required(),
    FieldDescriptor::builtin(
        b"STRUC",
        FieldCategory::Info,
        VcfType::String,
        FieldNumber::Fixed(1),
        "Structure of the region",
    ),
];

pub static FORMAT_FIELDS: &[FieldDescriptor<'static>] = &[
    FieldDescriptor::builtin(
        b"GT",
        FieldCategory::Format,
        VcfType::String,
        FieldNumber::Fixed(1),
        "Genotype",
    )
    .required()
    .with_special_handling(),
    FieldDescriptor::builtin(
        b"AL",
        FieldCategory::Format,
        VcfType::Integer,
        FieldNumber::Variable,
        "Length of each allele",
    )
    .required(),
    FieldDescriptor::builtin_with_aliases(
        b"ALLR",
        FieldCategory::Format,
        VcfType::String,
        FieldNumber::Variable,
        "Length range per allele",
        &[b"ALCI"],
    )
    .required(),
    FieldDescriptor::builtin(
        b"SD",
        FieldCategory::Format,
        VcfType::Integer,
        FieldNumber::Variable,
        "Number of spanning reads supporting per allele",
    )
    .required(),
    FieldDescriptor::builtin(
        b"MC",
        FieldCategory::Format,
        VcfType::String,
        FieldNumber::Variable,
        "Motif counts per allele",
    )
    .required(),
    FieldDescriptor::builtin(
        b"MS",
        FieldCategory::Format,
        VcfType::String,
        FieldNumber::Variable,
        "Motif spans per allele",
    )
    .required(),
    FieldDescriptor::builtin(
        b"AP",
        FieldCategory::Format,
        VcfType::Float,
        FieldNumber::Variable,
        "Allele purity per allele",
    )
    .required(),
    FieldDescriptor::builtin(
        b"AM",
        FieldCategory::Format,
        VcfType::Float,
        FieldNumber::Variable,
        "Mean methylation level per allele",
    )
    .required()
    .with_converter(convert_am_field),
    FieldDescriptor::builtin(
        b"PS",
        FieldCategory::Format,
        VcfType::Integer,
        FieldNumber::Fixed(1),
        "Phase set identifier",
    ),
    FieldDescriptor::builtin(
        b"PF",
        FieldCategory::Format,
        VcfType::String,
        FieldNumber::Fixed(1),
        "Phasing flag",
    ),
];

fn convert_am_field(record: &Record, field_id: &[u8]) -> Result<FieldValue> {
    if let Ok(float_field) = record.format(field_id).float() {
        let mut out = Vec::new();
        for sample_values in float_field.iter() {
            out.extend_from_slice(sample_values);
            if sample_values.len() <= 1 {
                out.push(VECTOR_END_FLOAT);
            }
        }
        return Ok(FieldValue::Float(out));
    }

    // Fall back to integer AM from old TRGT
    if let Ok(int_field) = record.format(field_id).integer() {
        let mut out = Vec::new();
        for sample_values in int_field.iter() {
            for &val in sample_values.iter() {
                if val == MISSING_INTEGER {
                    out.push(MISSING_FLOAT);
                } else {
                    // old AM was 0..255
                    out.push(val as f32 / 255.0);
                }
            }
            if sample_values.len() <= 1 {
                out.push(VECTOR_END_FLOAT);
            }
        }
        return Ok(FieldValue::Float(out));
    }

    Err(format!(
        "Could not read AM field as Float or Integer from record at position {}",
        record.pos()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::merge::field_registry::{
        MISSING_FLOAT, MISSING_INTEGER, VECTOR_END_FLOAT, VECTOR_END_INTEGER,
    };
    use crate::utils::test_util::{TestVcfBuilder, TestVcfRecord};
    use rust_htslib::bcf::Read;

    #[test]
    fn test_info_fields_defined() {
        let expected_fields: Vec<&[u8]> = vec![b"TRID", b"END", b"MOTIFS", b"STRUC"];
        let actual_fields: Vec<&[u8]> = INFO_FIELDS.iter().map(|f| f.id.as_ref()).collect();

        assert_eq!(
            actual_fields.len(),
            expected_fields.len(),
            "Expected {} INFO fields, found {}",
            expected_fields.len(),
            actual_fields.len()
        );

        for expected_id in &expected_fields {
            assert!(
                actual_fields.iter().any(|f| *f == *expected_id),
                "Missing expected INFO field: {}",
                String::from_utf8_lossy(expected_id)
            );
        }
    }

    #[test]
    fn test_format_fields_defined() {
        let expected_fields: Vec<&[u8]> = vec![
            b"GT", b"AL", b"ALLR", b"SD", b"MC", b"MS", b"AP", b"AM", b"PS", b"PF",
        ];
        let actual_fields: Vec<&[u8]> = FORMAT_FIELDS.iter().map(|f| f.id.as_ref()).collect();

        assert_eq!(
            actual_fields.len(),
            expected_fields.len(),
            "Expected {} FORMAT fields, found {}",
            expected_fields.len(),
            actual_fields.len()
        );

        for expected_id in &expected_fields {
            assert!(
                actual_fields.iter().any(|f| *f == *expected_id),
                "Missing expected FORMAT field: {}",
                String::from_utf8_lossy(expected_id)
            );
        }
    }

    #[test]
    fn test_allr_has_alci_alias() {
        let allr_field = FORMAT_FIELDS
            .iter()
            .find(|f| f.id.as_ref() == b"ALLR")
            .unwrap();
        let aliases: Vec<&[u8]> = allr_field.aliases.to_vec();
        assert_eq!(aliases.len(), 1);
        assert_eq!(aliases[0], b"ALCI");
    }

    #[test]
    fn test_am_has_converter() {
        let am_field = FORMAT_FIELDS
            .iter()
            .find(|f| f.id.as_ref() == b"AM")
            .unwrap();
        assert!(am_field.converter.is_some());
    }

    #[test]
    fn test_required_fields() {
        let required_info_fields: Vec<&[u8]> =
            vec![b"TRID".as_ref(), b"END".as_ref(), b"MOTIFS".as_ref()];
        for field_id in required_info_fields {
            let field = INFO_FIELDS
                .iter()
                .find(|f| f.id.as_ref() == field_id)
                .unwrap();
            assert!(
                field.required,
                "INFO field {} should be required",
                String::from_utf8_lossy(field.id.as_ref())
            );
        }

        let required_format_fields: Vec<&[u8]> = vec![
            b"GT".as_ref(),
            b"AL".as_ref(),
            b"ALLR".as_ref(),
            b"SD".as_ref(),
            b"MC".as_ref(),
            b"MS".as_ref(),
            b"AP".as_ref(),
            b"AM".as_ref(),
        ];
        for field_id in required_format_fields {
            let field = FORMAT_FIELDS
                .iter()
                .find(|f| f.id.as_ref() == field_id)
                .unwrap();
            assert!(
                field.required,
                "FORMAT field {} should be required",
                String::from_utf8_lossy(field.id.as_ref())
            );
        }
    }

    #[test]
    fn test_gt_has_special_handling() {
        let gt_field = FORMAT_FIELDS
            .iter()
            .find(|f| f.id.as_ref() == b"GT")
            .unwrap();
        assert!(
            gt_field.special_handling,
            "GT field should have special handling"
        );
        assert!(gt_field.required, "GT field should be required");
    }

    #[test]
    fn test_optional_fields() {
        let struc_field = INFO_FIELDS
            .iter()
            .find(|f| f.id.as_ref() == b"STRUC")
            .unwrap();
        assert!(!struc_field.required, "INFO field STRUC should be optional");

        let ps_field = FORMAT_FIELDS
            .iter()
            .find(|f| f.id.as_ref() == b"PS")
            .unwrap();
        assert!(!ps_field.required, "FORMAT field PS should be optional");

        let pf_field = FORMAT_FIELDS
            .iter()
            .find(|f| f.id.as_ref() == b"PF")
            .unwrap();
        assert!(!pf_field.required, "FORMAT field PF should be optional");
    }

    #[test]
    fn test_convert_am_field() {
        let temp_file = TestVcfBuilder::new()
            .contig("chr1", 100)
            .add_format("AM", ".", "Float", "Mean methylation level per allele")
            .sample("S1")
            .record(
                TestVcfRecord::new()
                    .alleles(&["A", "ATTT"])
                    .format_float("AM", &[0.5]),
            )
            .build();

        let mut reader = rust_htslib::bcf::Reader::from_path(temp_file.path()).unwrap();
        let mut record = reader.empty_record();
        let _ = reader.read(&mut record).unwrap();

        let field = convert_am_field(&record, b"AM").unwrap();
        let FieldValue::Float(values) = field else {
            panic!("Expected AM field to convert into Float values");
        };
        assert_eq!(values.len(), 2);
        assert!((values[0] - 0.5).abs() < f32::EPSILON);
        assert_eq!(values[1].to_bits(), VECTOR_END_FLOAT.to_bits());
    }

    #[test]
    fn test_convert_am_field_multisample() {
        let temp_file = TestVcfBuilder::new()
            .contig("chr1", 100)
            .add_format("AM", ".", "Float", "Mean methylation level per allele")
            .sample("S1")
            .sample("S2")
            .record(
                TestVcfRecord::new()
                    .alleles(&["A", "ATTT"])
                    .format_float("AM", &[0.2, 0.4, 0.9, VECTOR_END_FLOAT]),
            )
            .build();

        let mut reader = rust_htslib::bcf::Reader::from_path(temp_file.path()).unwrap();
        let mut record = reader.empty_record();
        let _ = reader.read(&mut record).unwrap();

        let field = convert_am_field(&record, b"AM").unwrap();
        let FieldValue::Float(values) = field else {
            panic!("Expected AM field to convert into Float values");
        };
        // S1 -> [0.2, 0.4], S2 -> [0.9, VECTOR_END]
        assert_eq!(values.len(), 4);
        assert!((values[0] - 0.2).abs() < f32::EPSILON);
        assert!((values[1] - 0.4).abs() < f32::EPSILON);
        assert!((values[2] - 0.9).abs() < f32::EPSILON);
        assert_eq!(values[3].to_bits(), VECTOR_END_FLOAT.to_bits());
    }

    #[test]
    fn convert_am_field_converts_integer() {
        let temp_file = TestVcfBuilder::new()
            .contig("chr1", 100)
            .add_format(
                "AM",
                ".",
                "Integer",
                "Mean methylation level per allele (0-255)",
            )
            .sample("S1")
            .sample("S2")
            .record(
                TestVcfRecord::new()
                    .alleles(&["A", "ATTT"])
                    .format_integer("AM", &[0, 255, MISSING_INTEGER, VECTOR_END_INTEGER]),
            )
            .build();

        let mut reader = rust_htslib::bcf::Reader::from_path(temp_file.path()).unwrap();
        let mut record = reader.empty_record();
        let _ = reader.read(&mut record).unwrap();

        let field = convert_am_field(&record, b"AM").unwrap();
        let FieldValue::Float(values) = field else {
            panic!("Expected AM field to convert into Float values");
        };
        // S1 -> [0, 255] -> [0.0, 1.0], S2 -> [.] -> [MISSING, VECTOR_END]
        assert_eq!(values.len(), 4);
        assert!((values[0] - 0.0).abs() < f32::EPSILON);
        assert!((values[1] - 1.0).abs() < f32::EPSILON);
        assert_eq!(values[2].to_bits(), MISSING_FLOAT.to_bits());
        assert_eq!(values[3].to_bits(), VECTOR_END_FLOAT.to_bits());
    }
}
