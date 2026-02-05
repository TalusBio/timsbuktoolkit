use crate::TimsElutionGroup;
use crate::ion::{
    IonAnnot,
    IonParsingError,
};
use serde::Deserialize;
use std::path::Path;
use tinyvec::tiny_vec;
use tracing::{
    error,
    info,
    warn,
};

#[derive(Debug)]
pub enum SpectronautReadingError {
    IoError,
    CsvError,
    SpectronautPrecursorParsingError,
}

/// Error type for library format detection (sniffing)
#[derive(Debug)]
pub enum LibrarySniffError {
    /// Failed to open or read the file
    IoError(String),
    /// Failed to parse CSV/TSV headers
    InvalidFormat(String),
    /// File is valid but missing required columns
    MissingColumns(Vec<String>),
}

impl std::fmt::Display for SpectronautReadingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SpectronautReadingError::IoError => write!(f, "IO error"),
            SpectronautReadingError::CsvError => write!(f, "CSV parsing error"),
            SpectronautReadingError::SpectronautPrecursorParsingError => {
                write!(f, "Spectronaut precursor parsing error")
            }
        }
    }
}

impl std::error::Error for SpectronautReadingError {}

#[derive(Debug)]
pub enum SpectronautPrecursorParsingError {
    IonParsingError,
    IonOverCapacity,
    EmptyIonString,
    Other,
}

impl From<IonParsingError> for SpectronautPrecursorParsingError {
    fn from(err: IonParsingError) -> Self {
        error!("Ion parsing error: {:?}", err);
        SpectronautPrecursorParsingError::IonParsingError
    }
}

impl From<SpectronautPrecursorParsingError> for SpectronautReadingError {
    fn from(_err: SpectronautPrecursorParsingError) -> Self {
        SpectronautReadingError::SpectronautPrecursorParsingError
    }
}

impl From<csv::Error> for SpectronautReadingError {
    fn from(err: csv::Error) -> Self {
        error!("CSV reading error: {:?}", err);
        SpectronautReadingError::CsvError
    }
}

impl From<std::io::Error> for SpectronautReadingError {
    fn from(err: std::io::Error) -> Self {
        error!("IO error: {:?}", err);
        SpectronautReadingError::IoError
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct SpectronautPrecursorExtras {
    pub modified_peptide: String,
    pub stripped_peptide: String,
    pub protein_id: String,
    pub is_decoy: bool,
    pub relative_intensities: Vec<(IonAnnot, f32)>,
}

/// Represents a single row from a Spectronaut library TSV file
#[derive(Debug, Clone, Deserialize)]
struct SpectronautLibraryRow {
    #[serde(rename = "ModifiedPeptide")]
    modified_peptide: String,
    #[serde(rename = "StrippedPeptide")]
    stripped_peptide: String,
    #[serde(rename = "PrecursorMz")]
    precursor_mz: f64,
    #[serde(rename = "PrecursorCharge")]
    precursor_charge: i32,
    #[serde(rename = "iRT")]
    irt: f64,
    #[serde(rename = "IonMobility")]
    ion_mobility: f64,
    #[serde(rename = "ProteinGroups")]
    protein_groups: String,
    #[serde(rename = "FragmentMz")]
    fragment_mz: f64,
    #[serde(rename = "FragmentType")]
    fragment_type: String,
    #[serde(rename = "FragmentNumber")]
    fragment_number: i32,
    #[serde(rename = "FragmentCharge")]
    fragment_charge: i32,
    #[serde(rename = "FragmentLossType")]
    fragment_loss_type: String,
    #[serde(rename = "RelativeIntensity")]
    relative_intensity: f32,
    #[serde(rename = "ExcludeFromAssay")]
    exclude_from_assay: String,
}

impl SpectronautLibraryRow {
    /// Check if this row belongs to the same precursor as another row
    fn is_same_precursor(&self, other: &SpectronautLibraryRow) -> bool {
        self.modified_peptide == other.modified_peptide
            && self.precursor_mz == other.precursor_mz
            && self.precursor_charge == other.precursor_charge
            && self.irt == other.irt
            && self.ion_mobility == other.ion_mobility
            && self.protein_groups == other.protein_groups
    }

    /// Check if this fragment should be excluded from the assay
    fn is_excluded(&self) -> bool {
        self.exclude_from_assay == "True"
    }
}

/// Check if a file is a Spectronaut library TSV by looking for Spectronaut-specific columns.
///
/// Returns `Ok(())` if the file matches Spectronaut format, or `Err` with details about why not.
pub fn sniff_spectronaut_library_file<T: AsRef<Path>>(file: T) -> Result<(), LibrarySniffError> {
    let file_handle = std::fs::File::open(file.as_ref()).map_err(|e| {
        LibrarySniffError::IoError(format!("Failed to open {}: {}", file.as_ref().display(), e))
    })?;

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file_handle);

    let headers = rdr.headers().map_err(|e| {
        LibrarySniffError::InvalidFormat(format!("Failed to parse TSV headers: {}", e))
    })?;

    let columns: Vec<String> = headers.iter().map(|s| s.to_string()).collect();

    // Required columns for Spectronaut libraries (includes Spectronaut-specific markers)
    let required_columns = [
        "ModifiedPeptide",
        "StrippedPeptide",
        "PrecursorMz",
        "PrecursorCharge",
        "iRT", // Spectronaut-specific (DIA-NN uses Tr_recalibrated)
        "IonMobility",
        "ProteinGroups",
        "FragmentMz",
        "FragmentType",
        "FragmentNumber",
        "FragmentCharge",
        "FragmentLossType",
        "RelativeIntensity",
        "ExcludeFromAssay", // Spectronaut-specific
    ];

    // Find missing columns
    let missing: Vec<String> = required_columns
        .iter()
        .filter(|col| !columns.contains(&col.to_string()))
        .map(|s| s.to_string())
        .collect();

    if missing.is_empty() {
        Ok(())
    } else {
        Err(LibrarySniffError::MissingColumns(missing))
    }
}

struct ParsingBuffers {
    fragment_labels: Vec<IonAnnot>,
}

pub fn read_library_file<T: AsRef<Path>>(
    file: T,
) -> Result<Vec<(TimsElutionGroup<IonAnnot>, SpectronautPrecursorExtras)>, SpectronautReadingError>
{
    let file_handle = std::fs::File::open(file.as_ref())?;

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file_handle);

    info!("Reading file content from {}", file.as_ref().display());

    let mut elution_groups = Vec::new();
    let mut current_group: Vec<SpectronautLibraryRow> = Vec::with_capacity(20);
    let mut buffers = ParsingBuffers {
        fragment_labels: Vec::with_capacity(20),
    };

    let mut group_id = 0;

    for result in rdr.deserialize() {
        let row: SpectronautLibraryRow = result?;

        if let Some(last_row) = current_group.first()
            && !row.is_same_precursor(last_row)
        {
            // New group found, parse the collected group
            if let Some(eg) = parse_precursor_group(group_id, &current_group, &mut buffers)? {
                elution_groups.push(eg);
                group_id += 1;
            }

            current_group.clear();
        }

        current_group.push(row);
    }

    // Process the last group
    if !current_group.is_empty() {
        if let Some(eg) = parse_precursor_group(group_id, &current_group, &mut buffers)? {
            elution_groups.push(eg);
        }
    }

    info!("Parsed {} elution groups", elution_groups.len());
    Ok(elution_groups)
}

fn parse_precursor_group(
    id: u64,
    rows: &[SpectronautLibraryRow],
    buffers: &mut ParsingBuffers,
) -> Result<
    Option<(TimsElutionGroup<IonAnnot>, SpectronautPrecursorExtras)>,
    SpectronautPrecursorParsingError,
> {
    if rows.is_empty() {
        error!("Empty precursor group encountered on {id}");
        return Err(SpectronautPrecursorParsingError::Other);
    }

    // Filter out excluded fragments
    let included_rows: Vec<&SpectronautLibraryRow> =
        rows.iter().filter(|r| !r.is_excluded()).collect();

    // If all fragments are excluded, skip this precursor
    if included_rows.is_empty() {
        warn!(
            "All fragments excluded for precursor group {}, skipping",
            id
        );
        return Ok(None);
    }

    let first_row = &rows[0];
    let mobility = first_row.ion_mobility as f32;

    // iRT values are in a similar scale to DIA-NN's Tr_recalibrated (loosely minutes)
    // Convert to seconds
    let rt_seconds = first_row.irt as f32 * 60.0;
    let precursor_mz = first_row.precursor_mz;
    let precursor_charge: u8 =
        first_row
            .precursor_charge
            .try_into()
            .map_err(|e: std::num::TryFromIntError| {
                error!("Failed to convert PrecursorCharge to u8: {:?}", e);
                SpectronautPrecursorParsingError::IonOverCapacity
            })?;

    let mut fragment_mzs = Vec::with_capacity(included_rows.len());
    buffers.fragment_labels.clear();
    let mut relative_intensities = Vec::with_capacity(included_rows.len());
    let mut num_unknown_losses = 0;

    for (i, row) in included_rows.iter().enumerate() {
        let fragment_mz = row.fragment_mz;
        let frag_charge: u8 =
            row.fragment_charge
                .try_into()
                .map_err(|e: std::num::TryFromIntError| {
                    error!("Failed to convert FragmentCharge to u8: {:?}", e);
                    SpectronautPrecursorParsingError::IonOverCapacity
                })?;
        let rel_intensity = row.relative_intensity;

        // Check for loss type - treat non-"noloss" as unknown ions (same as DIA-NN)
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
            SpectronautPrecursorParsingError::EmptyIonString
        })?;

        let frag_num = row.fragment_number.try_into().map_err(|_| {
            error!(
                "Invalid fragment number (I expect all of them < 255): {}",
                row.fragment_number
            );
            SpectronautPrecursorParsingError::IonOverCapacity
        })?;

        let ion_annot = IonAnnot::try_new(frag_char, Some(frag_num), frag_charge as i8, 0)?;

        buffers.fragment_labels.push(ion_annot);
        fragment_mzs.push(fragment_mz);
        relative_intensities.push((ion_annot, rel_intensity));
    }

    // Spectronaut libraries are target-only, so is_decoy is always false
    let is_decoy = false;

    let precursor_extras = SpectronautPrecursorExtras {
        modified_peptide: first_row.modified_peptide.clone(),
        stripped_peptide: first_row.stripped_peptide.clone(),
        protein_id: first_row.protein_groups.clone(),
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

    Ok(Some((eg, precursor_extras)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_sniff_spectronaut_library_file() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let file_path = PathBuf::from(manifest_dir)
            .join("tests")
            .join("spectronaut_io_files")
            .join("sample_lib.tsv");

        let result = sniff_spectronaut_library_file(&file_path);
        assert!(
            result.is_ok(),
            "File should be detected as Spectronaut library: {:?}",
            result.err()
        );

        // Cargo.toml should not be detected as Spectronaut library
        let file_path = PathBuf::from(manifest_dir).join("Cargo.toml");
        let result = sniff_spectronaut_library_file(file_path);
        assert!(
            result.is_err(),
            "Cargo.toml should not be detected as Spectronaut library"
        );
    }

    #[test]
    fn test_sniff_diann_not_spectronaut() {
        // DIA-NN files should not be detected as Spectronaut
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let file_path = PathBuf::from(manifest_dir)
            .join("tests")
            .join("diann_io_files")
            .join("sample_lib.txt");

        let result = sniff_spectronaut_library_file(file_path);
        assert!(
            result.is_err(),
            "DIA-NN library should not be detected as Spectronaut library"
        );
        // Verify we get helpful error context
        if let Err(LibrarySniffError::MissingColumns(cols)) = result {
            assert!(
                cols.contains(&"iRT".to_string()),
                "Should report missing iRT column"
            );
        }
    }

    #[test]
    fn test_read_library_file() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let file_path = PathBuf::from(manifest_dir)
            .join("tests")
            .join("spectronaut_io_files")
            .join("sample_lib.tsv");

        let result = read_library_file(&file_path);
        assert!(
            result.is_ok(),
            "Failed to read library file: {:?}",
            result.err()
        );

        let elution_groups = result.unwrap();

        // Sample file has 2 unique precursors:
        // KTVTAMDVVYALKR (charge 3) and MRECISIHVGQAGVQIGNACWELYCLEHGIQPDGQMPSDK (charge 5)
        assert_eq!(elution_groups.len(), 2, "Expected 2 elution groups");
    }

    #[test]
    fn test_fragment_annotations_parsed_correctly() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let file_path = PathBuf::from(manifest_dir)
            .join("tests")
            .join("spectronaut_io_files")
            .join("sample_lib.tsv");

        let mut elution_groups = read_library_file(file_path).expect("Failed to read library");
        elution_groups.sort_by(|a, b| a.0.rt_seconds().partial_cmp(&b.0.rt_seconds()).unwrap());

        // First precursor (KTVTAMDVVYALKR) has iRT=44.467922 -> rt_seconds ~ 2668
        // Second precursor (MRECISIHVGQAGVQIGNACWELYCLEHGIQPDGQMPSDK) has iRT=83.00864 -> rt_seconds ~ 4980

        let first_eg = &elution_groups[0].0;
        let second_eg = &elution_groups[1].0;

        // Check that the first elution group has the expected RT (approximately 44.467922 * 60)
        assert!(
            (first_eg.rt_seconds() - 44.467922 * 60.0).abs() < 0.01,
            "First elution group RT mismatch"
        );

        // Check that fragment count matches expected (excluded fragments should be filtered)
        // For KTVTAMDVVYALKR: 6 fragments with ExcludeFromAssay=False out of many total rows
        assert!(
            first_eg.fragment_count() >= 5,
            "Expected at least 5 fragments for first precursor, got {}",
            first_eg.fragment_count()
        );

        // Check second precursor RT
        assert!(
            (second_eg.rt_seconds() - 83.00864 * 60.0).abs() < 0.01,
            "Second elution group RT mismatch"
        );
    }

    #[test]
    fn test_exclude_from_assay_filtering() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let file_path = PathBuf::from(manifest_dir)
            .join("tests")
            .join("spectronaut_io_files")
            .join("sample_lib.tsv");

        let elution_groups = read_library_file(file_path).expect("Failed to read library");

        // Find the KTVTAMDVVYALKR precursor (first one)
        let (eg, extras) = &elution_groups[0];

        // Verify that the extras are populated correctly
        assert_eq!(extras.stripped_peptide, "KTVTAMDVVYALKR");
        assert!(
            !extras.is_decoy,
            "Spectronaut libraries should be target-only"
        );

        // The sample file has 36 rows for this precursor, but many have ExcludeFromAssay=True
        // Only 6 rows have ExcludeFromAssay=False (looking at the sample: rows 2,4,5,6,7,8)
        // So we should have 6 fragments
        assert_eq!(
            eg.fragment_count(),
            6,
            "Expected 6 non-excluded fragments for KTVTAMDVVYALKR"
        );

        // Verify that the relative intensities in extras match the fragment count
        assert_eq!(
            extras.relative_intensities.len(),
            eg.fragment_count(),
            "Relative intensities count should match fragment count"
        );
    }

    #[test]
    fn test_precursor_extras_populated() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let file_path = PathBuf::from(manifest_dir)
            .join("tests")
            .join("spectronaut_io_files")
            .join("sample_lib.tsv");

        let elution_groups = read_library_file(file_path).expect("Failed to read library");

        for (_, extras) in &elution_groups {
            // All Spectronaut library entries should be targets (not decoys)
            assert!(!extras.is_decoy);

            // Should have non-empty peptide sequences
            assert!(!extras.modified_peptide.is_empty());
            assert!(!extras.stripped_peptide.is_empty());

            // Should have protein ID
            assert!(!extras.protein_id.is_empty());

            // Should have relative intensities
            assert!(!extras.relative_intensities.is_empty());
        }
    }
}
