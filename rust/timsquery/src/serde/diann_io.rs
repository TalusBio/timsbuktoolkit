use crate::TimsElutionGroup;
use crate::ion::{
    IonAnnot,
    IonParsingError,
};
use serde::Deserialize;
use std::path::Path;
use tinyvec::tiny_vec;
use tracing::{
    debug,
    error,
    info,
    warn,
};
use arrow::array::{Float32Array, Int64Array, LargeStringArray, StringArray};

#[derive(Debug)]
pub enum DiannReadingError {
    IoError,
    CsvError,
    DiannPrecursorParsingError,
    ParquetError(String),
    ArrowError(String),
}

impl std::fmt::Display for DiannReadingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DiannReadingError::IoError => write!(f, "IO error"),
            DiannReadingError::CsvError => write!(f, "CSV parsing error"),
            DiannReadingError::DiannPrecursorParsingError => write!(f, "DIA-NN precursor parsing error"),
            DiannReadingError::ParquetError(msg) => write!(f, "Parquet error: {}", msg),
            DiannReadingError::ArrowError(msg) => write!(f, "Arrow error: {}", msg),
        }
    }
}

impl std::error::Error for DiannReadingError {}

#[derive(Debug)]
pub enum DiannPrecursorParsingError {
    IonParsingError,
    IonOverCapacity,
    EmptyIonString,
    Other,
}

impl From<IonParsingError> for DiannPrecursorParsingError {
    fn from(err: IonParsingError) -> Self {
        error!("Ion parsing error: {:?}", err);
        DiannPrecursorParsingError::IonParsingError
    }
}

impl From<DiannPrecursorParsingError> for DiannReadingError {
    fn from(_err: DiannPrecursorParsingError) -> Self {
        DiannReadingError::DiannPrecursorParsingError
    }
}

impl From<csv::Error> for DiannReadingError {
    fn from(err: csv::Error) -> Self {
        error!("CSV reading error: {:?}", err);
        DiannReadingError::CsvError
    }
}

impl From<std::io::Error> for DiannReadingError {
    fn from(err: std::io::Error) -> Self {
        error!("IO error: {:?}", err);
        DiannReadingError::IoError
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct DiannPrecursorExtras {
    pub modified_peptide: String,
    pub stripped_peptide: String,
    pub protein_id: String,
    pub is_decoy: bool,
    pub relative_intensities: Vec<(IonAnnot, f32)>,
}

/// Represents a single row from a DIA-NN library TSV file
#[derive(Debug, Clone, Deserialize)]
struct DiannLibraryRow {
    #[serde(rename = "ModifiedPeptide")]
    modified_peptide: String,
    #[serde(rename = "StrippedPeptide")]
    #[serde(alias = "PeptideSequence")]
    stripped_peptide: String,
    #[serde(rename = "PrecursorMz")]
    precursor_mz: f64,
    #[serde(rename = "PrecursorCharge")]
    precursor_charge: i32,
    #[serde(rename = "Tr_recalibrated")]
    tr_recalibrated: f64,
    #[serde(rename = "IonMobility")]
    ion_mobility: f64,
    #[serde(rename = "ProteinID")]
    #[serde(alias = "ProteinGroup")]
    protein_id: String,
    #[serde(rename = "Decoy")]
    #[serde(alias = "decoy")]
    decoy: i32,
    #[serde(rename = "FragmentMz")]
    #[serde(alias = "ProductMz")]
    fragment_mz: f64,
    #[serde(rename = "FragmentType")]
    fragment_type: String,
    #[serde(rename = "FragmentNumber")]
    #[serde(alias = "FragmentSeriesNumber")]
    fragment_number: i32,
    #[serde(rename = "FragmentCharge")]
    fragment_charge: i32,
    #[serde(rename = "FragmentLossType")]
    fragment_loss_type: String,
    #[serde(rename = "RelativeIntensity")]
    #[serde(alias = "LibraryIntensity")]
    relative_intensity: f32,
}

impl DiannLibraryRow {
    /// Check if this row belongs to the same precursor as another row
    /// This is much faster than generating the string key
    fn is_same_precursor(&self, other: &DiannLibraryRow) -> bool {
        self.modified_peptide == other.modified_peptide
            && self.precursor_mz == other.precursor_mz
            && self.precursor_charge == other.precursor_charge
            && self.tr_recalibrated == other.tr_recalibrated
            && self.ion_mobility == other.ion_mobility
            && self.protein_id == other.protein_id
            && self.decoy == other.decoy
    }
}

/// Check if a file is a DIA-NN parquet library
pub fn sniff_diann_parquet_library_file<T: AsRef<Path>>(file: T) -> bool {
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    let file_handle = match std::fs::File::open(file.as_ref()) {
        Ok(f) => f,
        Err(err) => {
            debug!(
                "Failed to open file {} for DIA-NN parquet sniffing: {:?}",
                file.as_ref().display(),
                err
            );
            return false;
        }
    };

    let builder = match ParquetRecordBatchReaderBuilder::try_new(file_handle) {
        Ok(b) => b,
        Err(_) => return false,
    };

    let schema = builder.schema();
    let field_names: Vec<String> = schema.fields().iter().map(|f| f.name().to_string()).collect();

    // Required columns for DiaNN 2.2 parquet format
    let required_columns = vec![
        "Modified.Sequence",
        "Stripped.Sequence",
        "Precursor.Charge",
        "RT",
        "IM",
        "Precursor.Mz",
        "Product.Mz",
        "Relative.Intensity",
        "Fragment.Type",
        "Fragment.Charge",
        "Fragment.Series.Number",
        "Fragment.Loss.Type",
        "Protein.Group",
        "Decoy",
    ];

    required_columns.iter().all(|col| field_names.contains(&col.to_string()))
}

pub fn sniff_diann_library_file<T: AsRef<Path>>(file: T) -> bool {
    let file_result = std::fs::File::open(file.as_ref());
    let file = match file_result {
        Ok(f) => f,
        Err(err) => {
            debug!(
                "Failed to open file {} for DIA-NN sniffing: {:?}",
                file.as_ref().display(),
                err
            );
            return false;
        }
    };

    let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').from_reader(file);

    // Get headers
    let headers = match rdr.headers() {
        Ok(h) => h,
        Err(_) => return false,
    };

    let columns: Vec<String> = headers.iter().map(|s| s.to_string()).collect();

    // Define required columns with their aliases
    let required_with_aliases = vec![
        vec!["ModifiedPeptide"],
        vec!["StrippedPeptide", "PeptideSequence"],
        vec!["PrecursorMz"],
        vec!["PrecursorCharge"],
        vec!["Tr_recalibrated"],
        vec!["IonMobility"],
        vec!["ProteinID", "ProteinGroup"],
        vec!["Decoy", "decoy"],
        vec!["FragmentMz", "ProductMz"],
        vec!["FragmentType"],
        vec!["FragmentNumber", "FragmentSeriesNumber"],
        vec!["FragmentCharge"],
        vec!["FragmentLossType"],
        vec!["RelativeIntensity", "LibraryIntensity"],
    ];

    // Check if all required columns (or their aliases) are present
    required_with_aliases
        .iter()
        .all(|aliases| aliases.iter().any(|col| columns.contains(&col.to_string())))
}

struct ParsingBuffers {
    fragment_labels: Vec<IonAnnot>,
}

/// Holds all the extracted column data from a parquet file
struct ParquetColumnData<'a> {
    modified_sequences: &'a [String],
    stripped_sequences: &'a [String],
    precursor_charges: &'a [i64],
    rts: &'a [f32],
    ims: &'a [f32],
    precursor_mzs: &'a [f32],
    product_mzs: &'a [f32],
    relative_intensities: &'a [f32],
    fragment_types: &'a [String],
    fragment_charges: &'a [i64],
    fragment_numbers: &'a [i64],
    fragment_loss_types: &'a [String],
    protein_groups: &'a [String],
    decoys: &'a [i64],
}

pub fn read_library_file<T: AsRef<Path>>(
    file: T,
) -> Result<Vec<(TimsElutionGroup<IonAnnot>, DiannPrecursorExtras)>, DiannReadingError> {
    let file_handle = std::fs::File::open(file.as_ref())?;

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file_handle);

    info!("Reading file content from {}", file.as_ref().display());

    // Optimization: Stream rows and group adjacent ones instead of loading everything into memory
    let mut elution_groups = Vec::new();
    let mut current_group: Vec<DiannLibraryRow> = Vec::with_capacity(20);
    // Reuse buffer for fragment labels to avoid allocation per group
    let mut buffers = ParsingBuffers {
        fragment_labels: Vec::with_capacity(20),
    };

    let mut group_id = 0;

    for result in rdr.deserialize() {
        let row: DiannLibraryRow = result?;

        if let Some(last_row) = current_group.first()
            && !row.is_same_precursor(last_row)
        {
            // New group found, parse the collected group
            let eg = parse_precursor_group(group_id, &current_group, &mut buffers)?;
            elution_groups.push(eg);

            // Start new group
            group_id += 1;
            current_group.clear();
        }

        current_group.push(row);
    }

    // Process the last group
    if !current_group.is_empty() {
        let eg = parse_precursor_group(group_id, &current_group, &mut buffers)?;
        elution_groups.push(eg);
    }

    info!("Parsed {} elution groups", elution_groups.len());
    Ok(elution_groups)
}

fn parse_precursor_group(
    id: u64,
    rows: &[DiannLibraryRow],
    buffers: &mut ParsingBuffers,
) -> Result<(TimsElutionGroup<IonAnnot>, DiannPrecursorExtras), DiannPrecursorParsingError> {
    if rows.is_empty() {
        error!("Empty precursor group encountered on {id}");
        return Err(DiannPrecursorParsingError::Other);
    }

    // All rows in this group share the same precursor info, so take first row
    let first_row = &rows[0];
    let mobility = first_row.ion_mobility as f32;

    // In diann retention times are usually in either biognosys-iRT times (which is loosely minutes)
    // or actual minutes - either way we convert to seconds here
    let rt_seconds = first_row.tr_recalibrated as f32 * 60.0;
    let precursor_mz = first_row.precursor_mz;
    let precursor_charge: u8 =
        first_row
            .precursor_charge
            .try_into()
            .map_err(|e: std::num::TryFromIntError| {
                error!("Failed to convert PrecursorCharge to u8: {:?}", e);
                DiannPrecursorParsingError::IonOverCapacity
            })?;

    // Extract fragment information from all rows
    let mut fragment_mzs = Vec::with_capacity(rows.len());
    buffers.fragment_labels.clear();
    let mut relative_intensities = Vec::with_capacity(rows.len());
    let mut num_unknown_losses = 0;

    for (i, row) in rows.iter().enumerate() {
        let fragment_mz = row.fragment_mz;
        let frag_charge: u8 =
            row.fragment_charge
                .try_into()
                .map_err(|e: std::num::TryFromIntError| {
                    error!("Failed to convert FragmentCharge to u8: {:?}", e);
                    DiannPrecursorParsingError::IonOverCapacity
                })?;
        let rel_intensity = row.relative_intensity;

        // Check for loss type - warn if not "noloss"
        if row.fragment_loss_type != "noloss" {
            warn!(
                "Unsupported fragment loss type '{}' at row {}; falling back to calling it an unknown ion",
                row.fragment_loss_type, i
            );

            num_unknown_losses += 1;
            let ion_annot = IonAnnot::try_new('?', Some(num_unknown_losses), frag_charge as i8, 0)?;
            buffers.fragment_labels.push(ion_annot);
            fragment_mzs.push(fragment_mz);
            relative_intensities.push((ion_annot, rel_intensity));
            continue;
        }

        let frag_char = row.fragment_type.chars().next().ok_or_else(|| {
            error!(
                "Empty FragmentType at row {}; cannot parse ion annotation",
                i
            );
            DiannPrecursorParsingError::EmptyIonString
        })?;

        let frag_num = row.fragment_number.try_into().map_err(|_| {
            error!(
                "Invalid fragment number (I expect all of then < 255): {}",
                row.fragment_number
            );
            DiannPrecursorParsingError::IonOverCapacity
        })?;

        let ion_annot = IonAnnot::try_new(frag_char, Some(frag_num), frag_charge as i8, 0)?;

        buffers.fragment_labels.push(ion_annot);
        fragment_mzs.push(fragment_mz);
        relative_intensities.push((ion_annot, rel_intensity));
    }

    let is_decoy = match first_row.decoy {
        0 => false,
        1 => true,
        other => {
            error!("Unexpected Decoy value: {}", other);
            return Err(DiannPrecursorParsingError::Other);
        }
    };

    let precursor_extras = DiannPrecursorExtras {
        modified_peptide: first_row.modified_peptide.clone(),
        stripped_peptide: first_row.stripped_peptide.clone(),
        protein_id: first_row.protein_id.clone(),
        is_decoy,
        relative_intensities,
    };

    let eg = TimsElutionGroup::builder()
        .id(id)
        .mobility_ook0(mobility)
        .rt_seconds(rt_seconds)
        .fragment_labels(buffers.fragment_labels.as_slice().into())
        .fragment_mzs(fragment_mzs)
        .precursor_labels(tiny_vec![0]) // Single monoisotopic precursor
        .precursor(precursor_mz, precursor_charge)
        .try_build()
        .unwrap();

    Ok((eg, precursor_extras))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_sniff_diann_library_file() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let file_path = PathBuf::from(manifest_dir)
            .join("tests")
            .join("diann_io_files")
            .join("sample_lib.txt");

        let result = sniff_diann_library_file(file_path);
        assert!(result, "File should be detected as DIA-NN library");

        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let file_path = PathBuf::from(manifest_dir).join("Cargo.toml");
        let result = sniff_diann_library_file(file_path);
        assert!(
            !result,
            "Cargo.toml should not be detected as DIA-NN library"
        );
    }

    #[test]
    fn test_read_library_file() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let file_path = PathBuf::from(manifest_dir)
            .join("tests")
            .join("diann_io_files")
            .join("sample_lib.txt");

        let result = read_library_file(file_path);
        assert!(
            result.is_ok(),
            "Failed to read library file: {:?}",
            result.err()
        );

        let mut elution_groups = result.unwrap();

        // Sample file has 2 unique precursors: MGRYSGK and HGDTGRR
        assert_eq!(elution_groups.len(), 2, "Expected 2 elution groups");

        // Find the MGRYSGK group (should be first if sorted, but let's be safe)
        // Since we don't have access to peptide sequence in TimsElutionGroup,
        // we identify by known properties from the sample data

        // MGRYSGK: PrecursorMz=399.699980472937, IonMobility=0.7825, Tr=3.78, 5 fragments
        // HGDTGRR: PrecursorMz=399.701900228177, IonMobility=0.7709, Tr=3.51, 4 fragments

        // Order is not stable, so lets sort by RT to have a predictable order
        elution_groups.sort_by(|a, b| a.0.rt_seconds().partial_cmp(&b.0.rt_seconds()).unwrap());

        // More specific assertions if we can identify groups by ID or other means:
        let mgrysgk = &elution_groups[1].0;
        assert_eq!(mgrysgk.fragment_count(), 5);
        assert!((mgrysgk.rt_seconds() - 3.78 * 60.).abs() < 0.01);
        assert!((mgrysgk.mobility_ook0() - 0.7825).abs() < 0.001);

        let hgdtgrr = &elution_groups[0].0;
        assert_eq!(hgdtgrr.fragment_count(), 4);
        assert!((hgdtgrr.rt_seconds() - 3.51 * 60.).abs() < 0.01);
        assert!((hgdtgrr.mobility_ook0() - 0.7709).abs() < 0.001);
    }

    #[test]
    fn test_fragment_annotations_parsed_correctly() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let file_path = PathBuf::from(manifest_dir)
            .join("tests")
            .join("diann_io_files")
            .join("sample_lib.txt");

        let mut elution_groups = read_library_file(file_path).expect("Failed to read library");
        elution_groups.sort_by(|a, b| a.0.rt_seconds().partial_cmp(&b.0.rt_seconds()).unwrap());

        // More specific assertions if we can identify groups by ID or other means:
        let mgrysgk = &elution_groups[1].0;
        let hgdtgrr = &elution_groups[0].0;

        let mut mgrysgk_expected_labels = vec!["y6", "b6", "b3", "b5", "y4"]
            .into_iter()
            .map(|s| IonAnnot::try_from(s).unwrap())
            .collect::<Vec<_>>();
        let mut actual_labels: Vec<IonAnnot> = mgrysgk
            .iter_fragments()
            .map(|(label, _mz)| label.clone())
            .collect();
        mgrysgk_expected_labels.sort();
        actual_labels.sort();
        assert_eq!(actual_labels, mgrysgk_expected_labels);

        let mut hgdtgrr_expected_labels = vec!["y6", "b3", "y4", "y5"]
            .into_iter()
            .map(|s| IonAnnot::try_from(s).unwrap())
            .collect::<Vec<_>>();
        let mut actual_labels: Vec<IonAnnot> = hgdtgrr
            .iter_fragments()
            .map(|(label, _mz)| label.clone())
            .collect();
        hgdtgrr_expected_labels.sort();
        actual_labels.sort();
        assert_eq!(actual_labels, hgdtgrr_expected_labels);
    }

    #[test]
    fn test_precursor_properties() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let file_path = PathBuf::from(manifest_dir)
            .join("tests")
            .join("diann_io_files")
            .join("sample_lib.txt");

        let _elution_groups = read_library_file(file_path).expect("Failed to read library");
        // TODO ... implement the actual assertions
    }

    #[test]
    fn test_speclib2_tsv_parsing() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let file_path = PathBuf::from(manifest_dir)
            .join("tests")
            .join("diann_io_files")
            .join("sample_lib.tsv");

        // Check that sniff detects it correctly
        let is_diann = sniff_diann_library_file(&file_path);
        assert!(
            is_diann,
            "speclib2 TSV file should be detected as DIA-NN library"
        );

        // This test mainly checks that the name aliases wotk correctly
        let _elution_groups = read_library_file(file_path).expect("Failed to read library");
    }

    #[test]
    fn test_sniff_parquet_library_file() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let file_path = PathBuf::from(manifest_dir)
            .join("tests")
            .join("diann_io_files")
            .join("sample_pq_speclib.parquet");

        let result = sniff_diann_parquet_library_file(file_path);
        assert!(result, "File should be detected as DIA-NN parquet library");
    }

    #[test]
    fn test_read_parquet_library_file() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let file_path = PathBuf::from(manifest_dir)
            .join("tests")
            .join("diann_io_files")
            .join("sample_pq_speclib.parquet");

        let result = read_parquet_library_file(file_path);
        assert!(
            result.is_ok(),
            "Failed to read parquet library file: {:?}",
            result.err()
        );

        let elution_groups = result.unwrap();

        // Sample file has 3 unique precursors
        assert_eq!(elution_groups.len(), 3, "Expected 3 elution groups");

        // Verify basic properties of first precursor
        let (first_eg, first_extras) = &elution_groups[0];
        assert_eq!(first_extras.modified_peptide, "GREEWESAALQNANTK");
        assert_eq!(first_extras.stripped_peptide, "GREEWESAALQNANTK");
        assert!(!first_extras.is_decoy);
        assert_eq!(first_eg.fragment_count(), 4);
    }
}

/// Read a DIA-NN spectral library from a parquet file (DiaNN 2.2+ format)
pub fn read_parquet_library_file<T: AsRef<Path>>(
    file: T,
) -> Result<Vec<(TimsElutionGroup<IonAnnot>, DiannPrecursorExtras)>, DiannReadingError> {
    use arrow::record_batch::RecordBatch;
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    let file_handle = std::fs::File::open(file.as_ref())?;

    let builder = ParquetRecordBatchReaderBuilder::try_new(file_handle)
        .map_err(|e| DiannReadingError::ParquetError(format!("Failed to create parquet reader: {}", e)))?;

    let reader = builder.build()
        .map_err(|e| DiannReadingError::ParquetError(format!("Failed to build parquet reader: {}", e)))?;

    info!("Reading parquet file content from {}", file.as_ref().display());

    // Read all batches
    let mut all_batches: Vec<RecordBatch> = Vec::new();
    for batch_result in reader {
        let batch = batch_result
            .map_err(|e| DiannReadingError::ParquetError(format!("Failed to read batch: {}", e)))?;
        all_batches.push(batch);
    }

    if all_batches.is_empty() {
        return Ok(Vec::new());
    }

    // Helper function to extract columns from batches
    fn get_string_column(batches: &[RecordBatch], column_name: &str) -> Result<Vec<String>, DiannReadingError> {
        let mut result = Vec::new();
        for batch in batches {
            let column = batch
                .column_by_name(column_name)
                .ok_or_else(|| DiannReadingError::ArrowError(format!("Column {} not found", column_name)))?;

            // Try LargeStringArray first (common in parquet files)
            if let Some(large_string_array) = column.as_any().downcast_ref::<LargeStringArray>() {
                result.extend(large_string_array.iter().map(|s| s.unwrap_or("").to_string()));
            } else if let Some(string_array) = column.as_any().downcast_ref::<StringArray>() {
                result.extend(string_array.iter().map(|s| s.unwrap_or("").to_string()));
            } else {
                return Err(DiannReadingError::ArrowError(format!("Column {} is not a string array", column_name)));
            }
        }
        Ok(result)
    }

    fn get_float_column(batches: &[RecordBatch], column_name: &str) -> Result<Vec<f32>, DiannReadingError> {
        let mut result = Vec::new();
        for batch in batches {
            let column = batch
                .column_by_name(column_name)
                .ok_or_else(|| DiannReadingError::ArrowError(format!("Column {} not found", column_name)))?;
            let float_array = column
                .as_any()
                .downcast_ref::<Float32Array>()
                .ok_or_else(|| DiannReadingError::ArrowError(format!("Column {} is not a float array", column_name)))?;
            result.extend(float_array.iter().map(|v| v.unwrap_or(0.0)));
        }
        Ok(result)
    }

    fn get_int_column(batches: &[RecordBatch], column_name: &str) -> Result<Vec<i64>, DiannReadingError> {
        let mut result = Vec::new();
        for batch in batches {
            let column = batch
                .column_by_name(column_name)
                .ok_or_else(|| DiannReadingError::ArrowError(format!("Column {} not found", column_name)))?;
            let int_array = column
                .as_any()
                .downcast_ref::<Int64Array>()
                .ok_or_else(|| DiannReadingError::ArrowError(format!("Column {} is not an int array", column_name)))?;
            result.extend(int_array.iter().map(|v| v.unwrap_or(0)));
        }
        Ok(result)
    }

    // Extract all columns
    let modified_sequences = get_string_column(&all_batches, "Modified.Sequence")?;
    let stripped_sequences = get_string_column(&all_batches, "Stripped.Sequence")?;
    let precursor_charges = get_int_column(&all_batches, "Precursor.Charge")?;
    let rts = get_float_column(&all_batches, "RT")?;
    let ims = get_float_column(&all_batches, "IM")?;
    let precursor_mzs = get_float_column(&all_batches, "Precursor.Mz")?;
    let product_mzs = get_float_column(&all_batches, "Product.Mz")?;
    let relative_intensities = get_float_column(&all_batches, "Relative.Intensity")?;
    let fragment_types = get_string_column(&all_batches, "Fragment.Type")?;
    let fragment_charges = get_int_column(&all_batches, "Fragment.Charge")?;
    let fragment_numbers = get_int_column(&all_batches, "Fragment.Series.Number")?;
    let fragment_loss_types = get_string_column(&all_batches, "Fragment.Loss.Type")?;
    let protein_groups = get_string_column(&all_batches, "Protein.Group")?;
    let decoys = get_int_column(&all_batches, "Decoy")?;

    let num_rows = modified_sequences.len();

    // Create column data struct
    let columns = ParquetColumnData {
        modified_sequences: &modified_sequences,
        stripped_sequences: &stripped_sequences,
        precursor_charges: &precursor_charges,
        rts: &rts,
        ims: &ims,
        precursor_mzs: &precursor_mzs,
        product_mzs: &product_mzs,
        relative_intensities: &relative_intensities,
        fragment_types: &fragment_types,
        fragment_charges: &fragment_charges,
        fragment_numbers: &fragment_numbers,
        fragment_loss_types: &fragment_loss_types,
        protein_groups: &protein_groups,
        decoys: &decoys,
    };

    // Group fragments by precursor (same logic as TSV reader)
    let mut elution_groups = Vec::new();
    let mut current_group_indices = Vec::new();
    let mut buffers = ParsingBuffers {
        fragment_labels: Vec::with_capacity(20),
    };
    let mut group_id = 0;

    for i in 0..num_rows {
        let is_new_precursor = if let Some(&last_idx) = current_group_indices.last() {
            columns.modified_sequences[i] != columns.modified_sequences[last_idx]
                || columns.precursor_mzs[i] != columns.precursor_mzs[last_idx]
                || columns.precursor_charges[i] != columns.precursor_charges[last_idx]
                || columns.rts[i] != columns.rts[last_idx]
                || columns.ims[i] != columns.ims[last_idx]
                || columns.protein_groups[i] != columns.protein_groups[last_idx]
                || columns.decoys[i] != columns.decoys[last_idx]
        } else {
            false
        };

        if is_new_precursor {
            // Parse the collected group
            let eg = parse_precursor_group_from_parquet(
                group_id,
                &current_group_indices,
                &columns,
                &mut buffers,
            )?;
            elution_groups.push(eg);

            group_id += 1;
            current_group_indices.clear();
        }

        current_group_indices.push(i);
    }

    // Process the last group
    if !current_group_indices.is_empty() {
        let eg = parse_precursor_group_from_parquet(
            group_id,
            &current_group_indices,
            &columns,
            &mut buffers,
        )?;
        elution_groups.push(eg);
    }

    info!("Parsed {} elution groups from parquet file", elution_groups.len());
    Ok(elution_groups)
}

fn parse_precursor_group_from_parquet(
    id: u64,
    indices: &[usize],
    columns: &ParquetColumnData,
    buffers: &mut ParsingBuffers,
) -> Result<(TimsElutionGroup<IonAnnot>, DiannPrecursorExtras), DiannPrecursorParsingError> {
    if indices.is_empty() {
        error!("Empty precursor group encountered on {id}");
        return Err(DiannPrecursorParsingError::Other);
    }

    let first_idx = indices[0];
    let mobility = columns.ims[first_idx];

    // RT is in minutes in parquet format, convert to seconds
    let rt_seconds = columns.rts[first_idx] * 60.0;
    let precursor_mz = columns.precursor_mzs[first_idx] as f64;
    let precursor_charge: u8 =
        columns.precursor_charges[first_idx]
            .try_into()
            .map_err(|e: std::num::TryFromIntError| {
                error!("Failed to convert PrecursorCharge to u8: {:?}", e);
                DiannPrecursorParsingError::IonOverCapacity
            })?;

    // Extract fragment information
    let mut fragment_mzs = Vec::with_capacity(indices.len());
    buffers.fragment_labels.clear();
    let mut rel_intensities = Vec::with_capacity(indices.len());
    let mut num_unknown_losses = 0;

    for (i, &idx) in indices.iter().enumerate() {
        let fragment_mz = columns.product_mzs[idx] as f64;
        let frag_charge: u8 =
            columns.fragment_charges[idx]
                .try_into()
                .map_err(|e: std::num::TryFromIntError| {
                    error!("Failed to convert FragmentCharge to u8: {:?}", e);
                    DiannPrecursorParsingError::IonOverCapacity
                })?;
        let rel_intensity = columns.relative_intensities[idx];

        // Check for loss type - warn if not "noloss" or "unknown"
        if columns.fragment_loss_types[idx] != "noloss" && columns.fragment_loss_types[idx] != "unknown" {
            warn!(
                "Unsupported fragment loss type '{}' at row {}; falling back to calling it an unknown ion",
                columns.fragment_loss_types[idx], i
            );

            num_unknown_losses += 1;
            let ion_annot = IonAnnot::try_new('?', Some(num_unknown_losses), frag_charge as i8, 0)?;
            buffers.fragment_labels.push(ion_annot);
            fragment_mzs.push(fragment_mz);
            rel_intensities.push((ion_annot, rel_intensity));
            continue;
        }

        let frag_char = columns.fragment_types[idx].chars().next().ok_or_else(|| {
            error!(
                "Empty FragmentType at row {}; cannot parse ion annotation",
                i
            );
            DiannPrecursorParsingError::EmptyIonString
        })?;

        let frag_num = columns.fragment_numbers[idx].try_into().map_err(|_| {
            error!(
                "Invalid fragment number (I expect all of then < 255): {}",
                columns.fragment_numbers[idx]
            );
            DiannPrecursorParsingError::IonOverCapacity
        })?;

        let ion_annot = IonAnnot::try_new(frag_char, Some(frag_num), frag_charge as i8, 0)?;

        buffers.fragment_labels.push(ion_annot);
        fragment_mzs.push(fragment_mz);
        rel_intensities.push((ion_annot, rel_intensity));
    }

    let is_decoy = match columns.decoys[first_idx] {
        0 => false,
        1 => true,
        other => {
            error!("Unexpected Decoy value: {}", other);
            return Err(DiannPrecursorParsingError::Other);
        }
    };

    let precursor_extras = DiannPrecursorExtras {
        modified_peptide: columns.modified_sequences[first_idx].clone(),
        stripped_peptide: columns.stripped_sequences[first_idx].clone(),
        protein_id: columns.protein_groups[first_idx].clone(),
        is_decoy,
        relative_intensities: rel_intensities,
    };

    let eg = TimsElutionGroup::builder()
        .id(id)
        .mobility_ook0(mobility)
        .rt_seconds(rt_seconds)
        .fragment_labels(buffers.fragment_labels.as_slice().into())
        .fragment_mzs(fragment_mzs)
        .precursor_labels(tiny_vec![0])
        .precursor(precursor_mz, precursor_charge)
        .try_build()
        .unwrap();

    Ok((eg, precursor_extras))
}
