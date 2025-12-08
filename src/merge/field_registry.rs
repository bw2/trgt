use super::vcf_writer::VcfWriter;
use crate::utils::Result;
use rust_htslib::{bcf::Record, htslib};
use std::{
    borrow::Cow,
    collections::{HashMap, HashSet},
};

pub const MISSING_INTEGER: i32 = i32::MIN;
pub const VECTOR_END_INTEGER: i32 = i32::MIN + 1;
pub const MISSING_FLOAT: f32 = f32::from_bits(0x7F80_0001);
pub const VECTOR_END_FLOAT: f32 = f32::from_bits(0x7F80_0002);

const EMPTY_ALIASES: &[&[u8]] = &[];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FieldCategory {
    Info,
    Format,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VcfType {
    Integer,
    Float,
    String,
    Character,
    Flag,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FieldNumber {
    Fixed(usize),  // Fixed number (e.g., 1, 2)
    AlleleCount,   // 'A'
    AlleleAndRef,  // 'R'
    GenotypeCount, // 'G'
    Variable,      // '.'
}

/// Cache for VCF field name to numeric ID mappings
#[derive(Debug, Default)]
pub struct FieldIdCache {
    cache: HashMap<Vec<u8>, Option<i32>>,
}

impl FieldIdCache {
    pub fn new() -> Self {
        Self {
            cache: HashMap::new(),
        }
    }

    /// Get the numeric field ID for a field name, using the cache if available
    #[inline(always)]
    pub fn get_or_lookup(
        &mut self,
        header: &rust_htslib::bcf::header::HeaderView,
        field_id: &[u8],
    ) -> Option<i32> {
        if let Some(&cached) = self.cache.get(field_id) {
            return cached;
        }

        let id = header.name_to_id(field_id).ok().map(|id| id.0 as i32);
        self.cache.insert(field_id.to_vec(), id);
        id
    }
}

#[derive(Debug, Clone)]
pub enum FieldValue {
    Integer(Vec<i32>),
    Float(Vec<f32>),
    String(Vec<Vec<u8>>),
    // Character is a single ASCII character. If Number=1 this will have len 1. If Number=. it can have more.
    Character(Vec<u8>),
    // Flag is INFO-only and presence/absence is the value (0 or 1).
    Flag(bool),
}

impl FieldValue {
    pub fn vcf_type(&self) -> VcfType {
        match self {
            FieldValue::Integer(_) => VcfType::Integer,
            FieldValue::Float(_) => VcfType::Float,
            FieldValue::String(_) => VcfType::String,
            FieldValue::Character(_) => VcfType::Character,
            FieldValue::Flag(_) => VcfType::Flag,
        }
    }

    pub fn clear(&mut self) {
        match self {
            FieldValue::Integer(v) => v.clear(),
            FieldValue::Float(v) => v.clear(),
            FieldValue::String(v) => v.clear(),
            FieldValue::Character(v) => v.clear(),
            FieldValue::Flag(v) => *v = false,
        }
    }

    pub fn push_missing_and_end(&mut self) {
        match self {
            FieldValue::Integer(v) => {
                v.push(MISSING_INTEGER);
                v.push(VECTOR_END_INTEGER);
            }
            FieldValue::Float(v) => {
                v.push(MISSING_FLOAT);
                v.push(VECTOR_END_FLOAT);
            }
            FieldValue::String(v) => {
                v.push(vec![b'.']);
            }
            FieldValue::Character(v) => {
                v.push(b'.');
            }
            FieldValue::Flag(_) => {
                // Flags are INFO-only, do nothing
            }
        }
    }

    pub fn push_missing_for_n(&mut self, n: usize) {
        for _ in 0..n {
            self.push_missing_and_end();
        }
    }
}

pub type FieldConverter = fn(&Record, &[u8]) -> Result<FieldValue>;
pub type FieldAliases<'a> = &'a [&'a [u8]];

// Two constructors for different use cases:
// - builtin(): Compile-time constants (zero allocations)
// - new(): Runtime-allocated (for custom VCF fields)
// This distinction is useful since built-in descriptors can be tagged with custom behavior (e.g., special handling of different versions),
// at the same time we need to be able to handle custom fields discovered from VCF headers, which should be the uncommon case
#[derive(Clone)]
pub struct FieldDescriptor<'a> {
    pub id: Cow<'a, [u8]>,
    pub category: FieldCategory,
    pub vcf_type: VcfType,
    pub number: FieldNumber,
    pub required: bool,
    pub description: Cow<'a, str>,
    // Alternative field IDs needed for older TRGT versions
    pub aliases: FieldAliases<'a>,
    // Custom converters needed for older TRGT versions
    pub converter: Option<FieldConverter>,
    // Tag for special handling (e.g., GT is written separately from the field data registry)
    pub special_handling: bool,
}

// Compile-time constants: zero allocations, for built-in fields defined in static arrays (see field_descriptors.rs)
impl<'a> FieldDescriptor<'a> {
    pub const fn builtin(
        id: &'a [u8],
        category: FieldCategory,
        vcf_type: VcfType,
        number: FieldNumber,
        description: &'a str,
    ) -> Self {
        FieldDescriptor {
            id: Cow::Borrowed(id),
            category,
            vcf_type,
            number,
            required: false,
            description: Cow::Borrowed(description),
            aliases: EMPTY_ALIASES,
            converter: None,
            special_handling: false,
        }
    }

    pub const fn builtin_with_aliases(
        id: &'a [u8],
        category: FieldCategory,
        vcf_type: VcfType,
        number: FieldNumber,
        description: &'a str,
        aliases: FieldAliases<'a>,
    ) -> Self {
        FieldDescriptor {
            id: Cow::Borrowed(id),
            category,
            vcf_type,
            number,
            required: false,
            description: Cow::Borrowed(description),
            aliases,
            converter: None,
            special_handling: false,
        }
    }

    pub const fn required(mut self) -> Self {
        self.required = true;
        self
    }

    pub const fn with_converter(mut self, converter: FieldConverter) -> Self {
        self.converter = Some(converter);
        self
    }

    pub const fn with_special_handling(mut self) -> Self {
        self.special_handling = true;
        self
    }

    pub fn create_empty_value(&self) -> FieldValue {
        match self.vcf_type {
            VcfType::Integer => FieldValue::Integer(Vec::new()),
            VcfType::Float => FieldValue::Float(Vec::new()),
            VcfType::String => FieldValue::String(Vec::new()),
            VcfType::Character => FieldValue::Character(Vec::new()),
            VcfType::Flag => FieldValue::Flag(false),
        }
    }

    pub fn create_missing_values(&self, n_samples: usize) -> FieldValue {
        let mut value = self.create_empty_value();
        for _ in 0..n_samples {
            value.push_missing_and_end();
        }
        value
    }

    /// Try to find the field ID in the record using a cache to avoid repeated FFI lookups.
    ///
    /// Returns the primary field ID if present, otherwise checks aliases. Returns None if
    /// neither the primary ID nor any alias is present in the record.
    pub fn try_read_field_id(&self, record: &Record, cache: &mut FieldIdCache) -> Option<&[u8]> {
        let primary = self.id.as_ref();
        if self.field_exists(record, primary, cache) {
            return Some(primary);
        }

        self.aliases
            .iter()
            .find(|&&alias| self.field_exists(record, alias, cache))
            .copied()
    }

    #[inline(always)]
    fn field_exists(&self, record: &Record, field_id: &[u8], cache: &mut FieldIdCache) -> bool {
        let Some(field_idx) = cache.get_or_lookup(record.header(), field_id) else {
            return false;
        };

        // SAFETY: `record.inner` is a htslib record, and we got `field_idx` from the cache
        // (which was originally obtained from the header)
        unsafe {
            match self.category {
                FieldCategory::Info => !htslib::bcf_get_info_id(record.inner, field_idx).is_null(),
                FieldCategory::Format => !htslib::bcf_get_fmt_id(record.inner, field_idx).is_null(),
            }
        }
    }
}

// Runtime-allocated, use for custom fields discovered from VCF headers
impl FieldDescriptor<'static> {
    pub fn new(
        id: &[u8],
        category: FieldCategory,
        vcf_type: VcfType,
        number: FieldNumber,
        description: &str,
    ) -> Self {
        FieldDescriptor {
            id: Cow::Owned(id.to_vec()),
            category,
            vcf_type,
            number,
            required: false,
            description: Cow::Owned(description.to_string()),
            aliases: EMPTY_ALIASES,
            converter: None,
            special_handling: false,
        }
    }

    pub fn borrowed_id(&self) -> Option<&'static [u8]> {
        match &self.id {
            Cow::Borrowed(slice) => Some(slice),
            Cow::Owned(_) => None,
        }
    }

    pub fn from_header_values(
        id: &str,
        type_str: &str,
        number_str: &str,
        description: &str,
        category: FieldCategory,
    ) -> Option<Self> {
        let vcf_type = match type_str {
            "String" => VcfType::String,
            "Integer" => VcfType::Integer,
            "Float" => VcfType::Float,
            "Character" => VcfType::Character,
            "Flag" => VcfType::Flag,
            _ => return None,
        };

        let number = match number_str {
            "." => FieldNumber::Variable,
            "A" => FieldNumber::AlleleCount,
            "R" => FieldNumber::AlleleAndRef,
            "G" => FieldNumber::GenotypeCount,
            n => {
                if let Ok(num) = n.parse::<usize>() {
                    FieldNumber::Fixed(num)
                } else {
                    return None;
                }
            }
        };

        if vcf_type == VcfType::Flag {
            if category != FieldCategory::Info {
                return None;
            }
            if !matches!(number, FieldNumber::Fixed(0)) {
                return None;
            }
        }

        Some(FieldDescriptor::new(
            id.as_bytes(),
            category,
            vcf_type,
            number,
            description,
        ))
    }
}

pub fn build_predefined_field_ids(
    descriptors: &'static [FieldDescriptor<'static>],
) -> HashSet<&'static [u8]> {
    let mut predefined_ids = HashSet::new();
    for descriptor in descriptors {
        let id = descriptor
            .borrowed_id()
            .expect("Expected borrowed identifier for predefined descriptor");
        predefined_ids.insert(id);
        for &alias in descriptor.aliases {
            predefined_ids.insert(alias);
        }
    }
    predefined_ids
}

pub struct FieldDataRegistry {
    id_to_index: HashMap<Vec<u8>, usize>,
    ordered_ids: Vec<Vec<u8>>,
    values: Vec<FieldValue>,
}

impl FieldDataRegistry {
    pub fn from_descriptors(descriptors: &[FieldDescriptor<'static>]) -> Self {
        let mut id_to_index = HashMap::with_capacity(descriptors.len());
        let mut ordered_ids = Vec::with_capacity(descriptors.len());
        let mut values = Vec::with_capacity(descriptors.len());

        for (idx, descriptor) in descriptors.iter().enumerate() {
            let id_vec = descriptor.id.as_ref().to_vec();
            id_to_index.insert(id_vec.clone(), idx);
            ordered_ids.push(id_vec);
            values.push(descriptor.create_empty_value());
        }

        FieldDataRegistry {
            id_to_index,
            ordered_ids,
            values,
        }
    }

    pub fn index_of(&self, field_id: &[u8]) -> Option<usize> {
        self.id_to_index.get(field_id).copied()
    }

    pub fn get(&self, field_id: &[u8]) -> Option<&FieldValue> {
        self.index_of(field_id).map(|idx| &self.values[idx])
    }

    pub fn get_mut(&mut self, field_id: &[u8]) -> Option<&mut FieldValue> {
        if let Some(idx) = self.index_of(field_id) {
            Some(&mut self.values[idx])
        } else {
            None
        }
    }

    pub fn get_by_index(&self, index: usize) -> &FieldValue {
        &self.values[index]
    }

    pub fn get_by_index_mut(&mut self, index: usize) -> &mut FieldValue {
        &mut self.values[index]
    }

    pub fn iter(&self) -> impl Iterator<Item = (&[u8], &FieldValue)> {
        self.ordered_ids
            .iter()
            .map(|id| id.as_slice())
            .zip(self.values.iter())
    }

    pub fn iter_mut_with_id(&mut self) -> impl Iterator<Item = (&[u8], &mut FieldValue)> {
        self.ordered_ids
            .iter()
            .map(|id| id.as_slice())
            .zip(self.values.iter_mut())
    }

    pub fn append_record(
        &mut self,
        record: &Record,
        descriptors: &[FieldDescriptor<'static>],
        cache: &mut FieldIdCache,
    ) -> Result<()> {
        let n = record.sample_count() as usize;
        for (idx, descriptor) in descriptors.iter().enumerate() {
            if let Err(err) = self.append_field(record, descriptor, idx, cache) {
                if descriptor.required {
                    let field_name = String::from_utf8_lossy(descriptor.id.as_ref());
                    return Err(format!(
                        "Failed to read required FORMAT field {}: {}",
                        field_name, err
                    ));
                }
                self.add_missing_values_for_field_idx(idx, n);
            }
        }
        Ok(())
    }

    pub fn push_missing_for_all_fields(&mut self, n: usize) {
        for value in &mut self.values {
            value.push_missing_for_n(n);
        }
    }

    pub fn write_format_fields(&self, writer: &mut VcfWriter) -> Result<()> {
        for (field_id, field_value) in self.iter() {
            Self::write_field(writer, field_id, field_value)?;
        }
        Ok(())
    }

    pub fn clear(&mut self) {
        for value in &mut self.values {
            value.clear();
        }
    }

    pub fn len(&self) -> usize {
        self.values.len()
    }

    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    fn append_field(
        &mut self,
        record: &Record,
        descriptor: &FieldDescriptor<'static>,
        descriptor_idx: usize,
        cache: &mut FieldIdCache,
    ) -> Result<()> {
        if let Some(converter) = descriptor.converter {
            let new_values = converter(record, descriptor.id.as_ref())?;
            let existing = self.get_by_index_mut(descriptor_idx);
            match (existing, new_values) {
                (FieldValue::Integer(e), FieldValue::Integer(n)) => e.extend(n),
                (FieldValue::Float(e), FieldValue::Float(n)) => e.extend(n),
                (FieldValue::String(e), FieldValue::String(n)) => e.extend(n),
                _ => return Err("Type mismatch between existing and converted field values".into()),
            }
            return Ok(());
        }

        let field_id = descriptor.try_read_field_id(record, cache).ok_or_else(|| {
            format!(
                "{} field not found",
                String::from_utf8_lossy(descriptor.id.as_ref())
            )
        })?;

        let target = self.get_by_index_mut(descriptor_idx);
        match (descriptor.vcf_type, target) {
            (VcfType::Integer, FieldValue::Integer(values)) => {
                extend_integer_field(record, field_id, values)?
            }
            (VcfType::Float, FieldValue::Float(values)) => {
                extend_float_field(record, field_id, values)?
            }
            (VcfType::String, FieldValue::String(values)) => {
                extend_string_field(record, field_id, values)?
            }
            (VcfType::Character, FieldValue::Character(values)) => {
                extend_character_field(record, field_id, values)?
            }
            (VcfType::Flag, FieldValue::Flag(flag)) => {
                *flag = record.info(descriptor.id.as_ref()).flag().unwrap_or(false);
            }
            _ => unreachable!("descriptor and registry are created from the same source"),
        }

        Ok(())
    }

    fn add_missing_values_for_field_idx(&mut self, descriptor_idx: usize, n: usize) {
        self.get_by_index_mut(descriptor_idx).push_missing_for_n(n);
    }

    fn write_field(
        writer: &mut VcfWriter,
        field_id: &[u8],
        field_value: &FieldValue,
    ) -> Result<()> {
        match field_value {
            FieldValue::Integer(values) => {
                writer
                    .dummy_record_mut()
                    .push_format_integer(field_id, values)
                    .map_err(|e| e.to_string())?;
            }
            FieldValue::Float(values) => {
                writer
                    .dummy_record_mut()
                    .push_format_float(field_id, values)
                    .map_err(|e| e.to_string())?;
            }
            FieldValue::String(values) => {
                writer
                    .dummy_record_mut()
                    .push_format_string(field_id, values)
                    .map_err(|e| e.to_string())?;
            }
            FieldValue::Character(values) => {
                let tmp: Vec<Vec<u8>> = values.iter().map(|c| vec![*c]).collect();
                writer
                    .dummy_record_mut()
                    .push_format_string(field_id, &tmp)
                    .map_err(|e| e.to_string())?;
            }
            FieldValue::Flag(flag) => {
                if *flag {
                    return Err("FORMAT Flag fields are not supported in VCF".to_string());
                }
            }
        }
        Ok(())
    }
}

fn extend_integer_field(record: &Record, field_id: &[u8], out: &mut Vec<i32>) -> Result<()> {
    let field = record
        .format(field_id)
        .integer()
        .map_err(|e| e.to_string())?;
    for sample_values in field.iter() {
        out.extend_from_slice(sample_values);
        if sample_values.len() <= 1 {
            out.push(VECTOR_END_INTEGER);
        }
    }
    Ok(())
}

fn extend_float_field(record: &Record, field_id: &[u8], out: &mut Vec<f32>) -> Result<()> {
    let field = record.format(field_id).float().map_err(|e| e.to_string())?;
    for sample_values in field.iter() {
        out.extend_from_slice(sample_values);
        if sample_values.len() <= 1 {
            out.push(VECTOR_END_FLOAT);
        }
    }
    Ok(())
}

fn extend_string_field(record: &Record, field_id: &[u8], out: &mut Vec<Vec<u8>>) -> Result<()> {
    let field = record
        .format(field_id)
        .string()
        .map_err(|e| e.to_string())?;
    for sample_value in field.iter() {
        out.push(sample_value.to_vec());
    }
    Ok(())
}

fn extend_character_field(record: &Record, field_id: &[u8], out: &mut Vec<u8>) -> Result<()> {
    let field = record
        .format(field_id)
        .string()
        .map_err(|e| e.to_string())?;
    for sample_value in field.iter() {
        out.push(sample_value.first().copied().unwrap_or(b'.'));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::test_util::{TestVcfBuilder, TestVcfRecord};
    use rust_htslib::bcf::{self, record::GenotypeAllele, Read, Reader};

    fn get_record_from_builder(builder: TestVcfBuilder) -> bcf::Record {
        let temp_vcf = builder.build();
        let reader = Reader::from_path(temp_vcf.path()).expect("Failed to open temp VCF");
        reader.empty_record()
    }

    #[test]
    fn test_field_value_all_types() {
        let mut int_value = FieldValue::Integer(vec![1, 2, 3]);
        assert_eq!(int_value.vcf_type(), VcfType::Integer);
        int_value.push_missing_and_end();
        if let FieldValue::Integer(v) = &int_value {
            assert_eq!(v.len(), 5);
            assert_eq!(v[3], MISSING_INTEGER);
            assert_eq!(v[4], VECTOR_END_INTEGER);
        } else {
            panic!("Expected Integer variant");
        }
        int_value.clear();
        if let FieldValue::Integer(v) = &int_value {
            assert_eq!(v.len(), 0);
        }

        let mut float_value = FieldValue::Float(vec![1.5]);
        assert_eq!(float_value.vcf_type(), VcfType::Float);
        float_value.push_missing_and_end();
        if let FieldValue::Float(v) = &float_value {
            assert_eq!(v.len(), 3);
            assert_eq!(v[1].to_bits(), MISSING_FLOAT.to_bits());
            assert_eq!(v[2].to_bits(), VECTOR_END_FLOAT.to_bits());
        } else {
            panic!("Expected Float variant");
        }
        float_value.clear();
        if let FieldValue::Float(v) = &float_value {
            assert!(v.is_empty());
        }

        let mut string_value = FieldValue::String(vec![b"abc".to_vec()]);
        assert_eq!(string_value.vcf_type(), VcfType::String);
        string_value.push_missing_and_end();
        if let FieldValue::String(v) = &string_value {
            assert_eq!(v.len(), 2);
            assert_eq!(v[1], vec![b'.']);
        } else {
            panic!("Expected String variant");
        }
        string_value.clear();
        if let FieldValue::String(v) = &string_value {
            assert!(v.is_empty());
        }
    }

    #[test]
    fn test_field_descriptor() {
        let descriptor = FieldDescriptor::new(
            b"AL",
            FieldCategory::Format,
            VcfType::Integer,
            FieldNumber::Variable,
            "Allele length",
        )
        .required();

        assert_eq!(descriptor.id.as_ref(), b"AL");
        assert!(descriptor.required);
        assert_eq!(descriptor.vcf_type, VcfType::Integer);

        let empty_value = descriptor.create_empty_value();
        assert!(matches!(empty_value, FieldValue::Integer(_)));

        let missing_values = descriptor.create_missing_values(2);
        if let FieldValue::Integer(v) = missing_values {
            assert_eq!(v.len(), 4); // 2 samples * (missing + vector_end)
        }
    }

    #[test]
    fn test_field_data_registry() {
        let descriptor = FieldDescriptor::new(
            b"AL",
            FieldCategory::Format,
            VcfType::Integer,
            FieldNumber::Variable,
            "Allele length",
        );

        let descriptors = [descriptor];
        let mut registry = FieldDataRegistry::from_descriptors(&descriptors);

        assert!(registry.get(b"AL").is_some());
        assert!(registry.get(b"MISSING").is_none());

        let idx = registry.index_of(b"AL").unwrap();
        if let FieldValue::Integer(values) = registry.get_by_index(idx) {
            assert!(values.is_empty());
        } else {
            panic!("Expected integer values");
        }

        if let FieldValue::Integer(values) = registry.get_by_index_mut(idx) {
            values.extend([1, 2, 3]);
        }

        if let Some(FieldValue::Integer(values)) = registry.get(b"AL") {
            assert_eq!(values.as_slice(), &[1, 2, 3]);
        } else {
            panic!("Expected AL field after mutation");
        }

        registry.clear();
        if let Some(FieldValue::Integer(v)) = registry.get(b"AL") {
            assert!(v.is_empty());
        } else {
            panic!("Expected AL field after clear");
        }
    }

    #[test]
    fn test_field_data_registry_iter_mut() {
        let descriptors = [
            FieldDescriptor::new(
                b"AL",
                FieldCategory::Format,
                VcfType::Integer,
                FieldNumber::Variable,
                "Allele length",
            ),
            FieldDescriptor::new(
                b"AP",
                FieldCategory::Format,
                VcfType::Float,
                FieldNumber::Variable,
                "Allele purity",
            ),
        ];

        let mut registry = FieldDataRegistry::from_descriptors(&descriptors);

        for (field_id, value) in registry.iter_mut_with_id() {
            if field_id == b"AL" {
                if let FieldValue::Integer(values) = value {
                    values.extend([5, 6]);
                } else {
                    panic!("Expected integer values for AL");
                }
            } else if field_id == b"AP" {
                if let FieldValue::Float(values) = value {
                    values.extend([0.1, 0.2]);
                } else {
                    panic!("Expected float values for AP");
                }
            } else {
                panic!("Unexpected field id");
            }
        }

        if let Some(FieldValue::Integer(values)) = registry.get(b"AL") {
            assert_eq!(values.as_slice(), &[5, 6]);
        } else {
            panic!("Expected AL field after iteration");
        }

        if let Some(FieldValue::Float(values)) = registry.get(b"AP") {
            assert_eq!(values.as_slice(), &[0.1, 0.2]);
        } else {
            panic!("Expected AP field after iteration");
        }
    }

    #[test]
    fn test_append_record_missing_optional_field() {
        let temp_file = TestVcfBuilder::new()
            .contig("chr1", 1)
            .add_format("REQUIRED", ".", "Integer", "Required")
            .add_format("OPTIONAL", ".", "String", "Optional")
            .sample("Sample1")
            .sample("Sample2")
            .record(
                TestVcfRecord::new()
                    .pos(999)
                    .alleles(&["A", "AT"])
                    .genotype(&[
                        GenotypeAllele::Unphased(0),
                        GenotypeAllele::Unphased(1),
                        GenotypeAllele::Unphased(0),
                        GenotypeAllele::Unphased(1),
                    ])
                    .format_integer("REQUIRED", &[10, 15, 20, 25]),
            )
            .build();

        let mut reader = Reader::from_path(temp_file.path()).unwrap();
        let mut record = reader.empty_record();
        let _ = reader.read(&mut record).unwrap();

        let descriptors = vec![
            FieldDescriptor::builtin(
                b"REQUIRED",
                FieldCategory::Format,
                VcfType::Integer,
                FieldNumber::Variable,
                "Required",
            )
            .required(),
            FieldDescriptor::builtin(
                b"OPTIONAL",
                FieldCategory::Format,
                VcfType::String,
                FieldNumber::Variable,
                "Optional",
            ),
        ];

        let mut registry = FieldDataRegistry::from_descriptors(&descriptors);
        let mut cache = FieldIdCache::new();
        registry
            .append_record(&record, &descriptors, &mut cache)
            .expect("no error");

        if let Some(FieldValue::Integer(values)) = registry.get(b"REQUIRED") {
            assert_eq!(values, &[10, 15, 20, 25]);
        } else {
            panic!("REQUIRED missing");
        }

        if let Some(FieldValue::String(values)) = registry.get(b"OPTIONAL") {
            assert_eq!(values.len(), 2);
            assert!(values.iter().all(|v| v.as_slice() == b"."));
        } else {
            panic!("OPTIONAL missing");
        }
    }

    #[test]
    fn test_try_read_field_id_missing() {
        // Field is defined in header but not set in record
        let record = get_record_from_builder(
            TestVcfBuilder::new().add_info("PRIMARY", "1", "Integer", "Primary"),
        );

        let descriptor = FieldDescriptor::new(
            b"PRIMARY",
            FieldCategory::Info,
            VcfType::Integer,
            FieldNumber::Fixed(1),
            "Primary",
        );
        let mut cache = FieldIdCache::new();
        assert!(descriptor.try_read_field_id(&record, &mut cache).is_none());

        // Field is not defined in header at all
        let record2 = get_record_from_builder(TestVcfBuilder::new());
        let descriptor2 = FieldDescriptor::new(
            b"PRIMARY",
            FieldCategory::Info,
            VcfType::Integer,
            FieldNumber::Fixed(1),
            "Primary",
        );
        let mut cache2 = FieldIdCache::new();
        assert!(descriptor2
            .try_read_field_id(&record2, &mut cache2)
            .is_none());
    }

    #[test]
    fn test_field_descriptor_builder_helpers() {
        fn dummy_converter(_record: &Record, _field_id: &[u8]) -> Result<FieldValue> {
            Ok(FieldValue::Integer(vec![42]))
        }

        let descriptor = FieldDescriptor::new(
            b"DP",
            FieldCategory::Info,
            VcfType::Integer,
            FieldNumber::Fixed(1),
            "Depth",
        )
        .with_converter(dummy_converter)
        .with_special_handling();

        assert!(descriptor.aliases.is_empty());
        assert!(descriptor.converter.is_some());
        assert!(descriptor.special_handling);
    }

    #[test]
    fn test_uses_alias_when_primary_missing() {
        static TEST_ALIASES: &[&[u8]] = &[b"ALCI"];
        let descriptor = FieldDescriptor::builtin_with_aliases(
            b"ALLR",
            FieldCategory::Format,
            VcfType::String,
            FieldNumber::Variable,
            "Allele length range",
            TEST_ALIASES,
        );

        let mut record = get_record_from_builder(
            TestVcfBuilder::new()
                .add_format("ALCI", ".", "String", "Alias") // Matches alias
                .sample("S1"),
        );

        let format_data = [&b"10-20"[..]];
        record.push_format_string(b"ALCI", &format_data).unwrap();

        let mut cache = FieldIdCache::new();
        assert_eq!(
            descriptor.try_read_field_id(&record, &mut cache),
            Some(b"ALCI".as_slice())
        );
    }

    #[test]
    fn test_field_id_cache() {
        let mut cache = FieldIdCache::new();
        assert!(cache.cache.is_empty());

        let temp_vcf = TestVcfBuilder::new()
            .add_info("EXISTS", "1", "Integer", "Existing field")
            .build();
        let reader = Reader::from_path(temp_vcf.path()).unwrap();
        let header = reader.header();

        // Missing field returns None and caches it
        let missing_result = cache.get_or_lookup(header, b"MISSING");
        assert!(missing_result.is_none());
        assert!(cache.cache.get(b"MISSING".as_slice()).unwrap().is_none());

        // Existing field returns Some and caches it
        let exists_result = cache.get_or_lookup(header, b"EXISTS");
        assert!(exists_result.is_some());
        let first_id = exists_result.unwrap();

        // Second lookup returns same cached value
        let second_result = cache.get_or_lookup(header, b"EXISTS");
        assert_eq!(first_id, second_result.unwrap());
    }

    #[test]
    fn test_try_read_field_id_info_and_format() {
        // Test INFO field
        let mut info_record = get_record_from_builder(
            TestVcfBuilder::new().add_info("PRIMARY", "1", "Integer", "Primary"),
        );
        info_record.push_info_integer(b"PRIMARY", &[5]).unwrap();

        let info_descriptor = FieldDescriptor::new(
            b"PRIMARY",
            FieldCategory::Info,
            VcfType::Integer,
            FieldNumber::Fixed(1),
            "Primary",
        );

        let mut cache = FieldIdCache::new();
        assert_eq!(
            info_descriptor.try_read_field_id(&info_record, &mut cache),
            Some(b"PRIMARY".as_slice())
        );

        // Test FORMAT field
        let mut format_record = get_record_from_builder(
            TestVcfBuilder::new()
                .add_format("FMT", "1", "String", "Format")
                .sample("S1"),
        );
        format_record
            .push_format_string(b"FMT", &[&b"A"[..]])
            .unwrap();

        let format_descriptor = FieldDescriptor::new(
            b"FMT",
            FieldCategory::Format,
            VcfType::String,
            FieldNumber::Fixed(1),
            "Format",
        );

        let mut cache2 = FieldIdCache::new();
        assert_eq!(
            format_descriptor.try_read_field_id(&format_record, &mut cache2),
            Some(b"FMT".as_slice())
        );
    }
}
