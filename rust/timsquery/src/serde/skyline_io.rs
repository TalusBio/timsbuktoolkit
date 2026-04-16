use crate::TimsElutionGroup;
use crate::ion::{
    IonAnnot,
    IonParsingError,
};
use serde::{
    Deserialize,
    Deserializer,
};
use std::path::Path;
use tinyvec::tiny_vec;
use tracing::{
    debug,
    error,
    info,
    warn,
};

#[derive(Debug)]
pub enum SkylineReadingError {
    IoError,
    CsvError,
    SkylinePrecursorParsingError,
}

#[derive(Debug)]
#[allow(dead_code)] // fields surfaced only through Debug for tracing
pub enum SkylineSniffError {
    IoError(String),
    InvalidFormat(String),
    MissingColumns(Vec<String>),
}

impl std::fmt::Display for SkylineReadingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SkylineReadingError::IoError => write!(f, "IO error"),
            SkylineReadingError::CsvError => write!(f, "CSV parsing error"),
            SkylineReadingError::SkylinePrecursorParsingError => {
                write!(f, "Skyline precursor parsing error")
            }
        }
    }
}

impl std::error::Error for SkylineReadingError {}

#[derive(Debug)]
pub enum SkylinePrecursorParsingError {
    IonParsingError,
    IonOverCapacity,
    Other,
}

impl From<IonParsingError> for SkylinePrecursorParsingError {
    fn from(err: IonParsingError) -> Self {
        error!("Ion parsing error: {:?}", err);
        SkylinePrecursorParsingError::IonParsingError
    }
}

impl From<SkylinePrecursorParsingError> for SkylineReadingError {
    fn from(_err: SkylinePrecursorParsingError) -> Self {
        SkylineReadingError::SkylinePrecursorParsingError
    }
}

impl From<csv::Error> for SkylineReadingError {
    fn from(err: csv::Error) -> Self {
        error!("CSV reading error: {:?}", err);
        SkylineReadingError::CsvError
    }
}

impl From<std::io::Error> for SkylineReadingError {
    fn from(err: std::io::Error) -> Self {
        error!("IO error: {:?}", err);
        SkylineReadingError::IoError
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct SkylinePrecursorExtras {
    pub modified_peptide: String,
    pub stripped_peptide: String,
    pub protein_id: String,
    pub is_decoy: bool,
    pub relative_intensities: Vec<(IonAnnot, f32)>,
}

/// Deserializer that treats Skyline's `#N/A` sentinel as `None`.
fn na_to_none<'de, D>(deserializer: D) -> Result<Option<f32>, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    let trimmed = s.trim();
    if trimmed.is_empty() || trimmed == "#N/A" {
        return Ok(None);
    }
    trimmed.parse::<f32>().map(Some).map_err(serde::de::Error::custom)
}

/// Represents a single row from a Skyline Peptide Transition List CSV.
#[derive(Debug, Clone, Deserialize)]
struct SkylineLibraryRow {
    #[serde(rename = "Protein Name")]
    protein_name: String,
    #[serde(rename = "Peptide Modified Sequence")]
    peptide_modified_sequence: String,
    #[serde(rename = "Precursor Mz")]
    precursor_mz: f64,
    #[serde(rename = "Precursor Charge")]
    precursor_charge: i32,
    #[serde(rename = "Product Mz")]
    product_mz: f64,
    #[serde(rename = "Product Charge")]
    product_charge: i32,
    #[serde(rename = "Fragment Ion Type")]
    fragment_ion_type: String,
    #[serde(rename = "Fragment Ion Ordinal")]
    fragment_ion_ordinal: i32,
    #[serde(rename = "Library Intensity", default, deserialize_with = "na_to_none")]
    library_intensity: Option<f32>,
    #[serde(rename = "Transition Is Decoy")]
    transition_is_decoy: String,
}

impl SkylineLibraryRow {
    fn is_same_precursor(&self, other: &SkylineLibraryRow) -> bool {
        self.peptide_modified_sequence == other.peptide_modified_sequence
            && self.precursor_charge == other.precursor_charge
            && self.precursor_mz == other.precursor_mz
            && self.protein_name == other.protein_name
            && self.transition_is_decoy == other.transition_is_decoy
    }

    fn is_precursor_row(&self) -> bool {
        self.fragment_ion_type.eq_ignore_ascii_case("precursor")
    }

    fn is_decoy(&self) -> Result<bool, SkylinePrecursorParsingError> {
        match self.transition_is_decoy.trim() {
            "True" | "true" | "TRUE" | "1" => Ok(true),
            "False" | "false" | "FALSE" | "0" => Ok(false),
            other => {
                error!("Unexpected `Transition Is Decoy` value: {}", other);
                Err(SkylinePrecursorParsingError::Other)
            }
        }
    }
}

/// Remove bracketed modification annotations, e.g. `C[+57.02]AM` -> `CAM`.
fn strip_modifications(modified_seq: &str) -> String {
    let mut out = String::with_capacity(modified_seq.len());
    let mut depth: i32 = 0;
    for ch in modified_seq.chars() {
        match ch {
            '[' | '(' | '{' => depth += 1,
            ']' | ')' | '}' => {
                if depth > 0 {
                    depth -= 1;
                }
            }
            _ if depth == 0 => out.push(ch),
            _ => {}
        }
    }
    out
}

/// Check if a file is a Skyline Peptide Transition List CSV.
pub fn sniff_skyline_library_file<T: AsRef<Path>>(file: T) -> Result<(), SkylineSniffError> {
    let file_handle = std::fs::File::open(file.as_ref()).map_err(|e| {
        SkylineSniffError::IoError(format!("Failed to open {}: {}", file.as_ref().display(), e))
    })?;

    let mut rdr = csv::ReaderBuilder::new().delimiter(b',').from_reader(file_handle);

    let headers = rdr.headers().map_err(|e| {
        SkylineSniffError::InvalidFormat(format!("Failed to parse CSV headers: {}", e))
    })?;

    let columns: Vec<String> = headers.iter().map(|s| s.to_string()).collect();

    // Required columns for a Skyline transition list. The combination of
    // "Peptide Modified Sequence", "Fragment Ion Type/Ordinal" and
    // "Transition Is Decoy" is distinctive to Skyline exports.
    let required_columns = [
        "Protein Name",
        "Peptide Modified Sequence",
        "Precursor Mz",
        "Precursor Charge",
        "Product Mz",
        "Product Charge",
        "Fragment Ion Type",
        "Fragment Ion Ordinal",
        "Transition Is Decoy",
    ];

    let missing: Vec<String> = required_columns
        .iter()
        .filter(|col| !columns.contains(&col.to_string()))
        .map(|s| s.to_string())
        .collect();

    if missing.is_empty() {
        Ok(())
    } else {
        Err(SkylineSniffError::MissingColumns(missing))
    }
}

struct ParsingBuffers {
    fragment_labels: Vec<IonAnnot>,
}

pub fn read_library_file<T: AsRef<Path>>(
    file: T,
) -> Result<Vec<(TimsElutionGroup<IonAnnot>, SkylinePrecursorExtras)>, SkylineReadingError> {
    let file_handle = std::fs::File::open(file.as_ref())?;

    let mut rdr = csv::ReaderBuilder::new().delimiter(b',').from_reader(file_handle);

    info!("Reading Skyline transition list from {}", file.as_ref().display());
    warn!(
        "Skyline transition lists do not carry retention time or ion mobility; \
         falling back to 0.0. Use an RT-unrestricted search or calibrate separately."
    );

    let mut elution_groups = Vec::new();
    let mut current_group: Vec<SkylineLibraryRow> = Vec::with_capacity(20);
    let mut buffers = ParsingBuffers {
        fragment_labels: Vec::with_capacity(20),
    };
    let mut group_id = 0u64;

    for result in rdr.deserialize() {
        let row: SkylineLibraryRow = result?;

        if let Some(last_row) = current_group.first()
            && !row.is_same_precursor(last_row)
        {
            if let Some(eg) = parse_precursor_group(group_id, &current_group, &mut buffers)? {
                elution_groups.push(eg);
                group_id += 1;
            }
            current_group.clear();
        }

        current_group.push(row);
    }

    if !current_group.is_empty()
        && let Some(eg) = parse_precursor_group(group_id, &current_group, &mut buffers)?
    {
        elution_groups.push(eg);
    }

    info!("Parsed {} elution groups", elution_groups.len());
    Ok(elution_groups)
}

fn parse_precursor_group(
    id: u64,
    rows: &[SkylineLibraryRow],
    buffers: &mut ParsingBuffers,
) -> Result<
    Option<(TimsElutionGroup<IonAnnot>, SkylinePrecursorExtras)>,
    SkylinePrecursorParsingError,
> {
    if rows.is_empty() {
        error!("Empty precursor group encountered on {id}");
        return Err(SkylinePrecursorParsingError::Other);
    }

    // Fragments are every non-precursor row. Precursor (isotope) rows are
    // ignored — downstream computes the isotope envelope from the peptide
    // sequence, matching the DIA-NN path.
    let fragment_rows: Vec<&SkylineLibraryRow> =
        rows.iter().filter(|r| !r.is_precursor_row()).collect();

    if fragment_rows.is_empty() {
        warn!(
            "Skyline precursor group {} has no product-ion rows; skipping",
            id
        );
        return Ok(None);
    }

    let first_row = &rows[0];
    let precursor_mz = first_row.precursor_mz;
    let precursor_charge: u8 =
        first_row
            .precursor_charge
            .try_into()
            .map_err(|e: std::num::TryFromIntError| {
                error!("Failed to convert PrecursorCharge to u8: {:?}", e);
                SkylinePrecursorParsingError::IonOverCapacity
            })?;
    let is_decoy = first_row.is_decoy()?;

    let mut fragment_mzs = Vec::with_capacity(fragment_rows.len());
    buffers.fragment_labels.clear();
    let mut relative_intensities = Vec::with_capacity(fragment_rows.len());
    let mut num_unknown_losses = 0u8;

    for (i, row) in fragment_rows.iter().enumerate() {
        let fragment_mz = row.product_mz;
        let frag_charge: u8 =
            row.product_charge
                .try_into()
                .map_err(|e: std::num::TryFromIntError| {
                    error!("Failed to convert Product Charge to u8: {:?}", e);
                    SkylinePrecursorParsingError::IonOverCapacity
                })?;
        // Skyline transition lists often lack library intensities (#N/A).
        // Default to 1.0; downstream normalization handles relative scaling.
        let rel_intensity = row.library_intensity.unwrap_or(1.0);

        let frag_type = row.fragment_ion_type.trim();
        // Known backbone series: a/b/c/x/y/z. Anything else -> unknown ion.
        let frag_char_opt = frag_type
            .chars()
            .next()
            .filter(|c| matches!(c, 'a' | 'b' | 'c' | 'x' | 'y' | 'z'));

        let ion_annot = match frag_char_opt {
            Some(frag_char) => {
                let frag_num: u8 = row.fragment_ion_ordinal.try_into().map_err(|_| {
                    error!(
                        "Invalid fragment ion ordinal (expected < 255): {}",
                        row.fragment_ion_ordinal
                    );
                    SkylinePrecursorParsingError::IonOverCapacity
                })?;
                IonAnnot::try_new(frag_char, Some(frag_num), frag_charge as i8, 0)?
            }
            None => {
                warn!(
                    "Unsupported Skyline fragment ion type '{}' at row {}; \
                     falling back to unknown ion",
                    frag_type, i
                );
                num_unknown_losses = num_unknown_losses.saturating_add(1);
                IonAnnot::try_new('?', Some(num_unknown_losses), frag_charge as i8, 0)?
            }
        };

        buffers.fragment_labels.push(ion_annot);
        fragment_mzs.push(fragment_mz);
        relative_intensities.push((ion_annot, rel_intensity));
    }

    let modified_peptide = first_row.peptide_modified_sequence.clone();
    let stripped_peptide = strip_modifications(&modified_peptide);

    let precursor_extras = SkylinePrecursorExtras {
        modified_peptide,
        stripped_peptide,
        protein_id: first_row.protein_name.clone(),
        is_decoy,
        relative_intensities,
    };

    let eg = TimsElutionGroup::builder()
        .id(id)
        .mobility_ook0(0.0)
        .rt_seconds(0.0)
        .fragment_labels(buffers.fragment_labels.as_slice().into())
        .fragment_mzs(fragment_mzs)
        .precursor_labels(tiny_vec![0])
        .precursor(precursor_mz, precursor_charge)
        .try_build()
        .unwrap();

    debug!(
        "Parsed Skyline precursor group {}: {} fragments, decoy={}",
        id,
        buffers.fragment_labels.len(),
        is_decoy
    );

    Ok(Some((eg, precursor_extras)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn fixture_path() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("skyline_io_files")
            .join("sample_transition_list.csv")
    }

    #[test]
    fn test_strip_modifications() {
        assert_eq!(strip_modifications("PEPTIDE"), "PEPTIDE");
        assert_eq!(strip_modifications("C[+57.021]AM"), "CAM");
        assert_eq!(strip_modifications("P[UniMod:35]IDE"), "PIDE");
        assert_eq!(strip_modifications("[+42]AB"), "AB");
    }

    #[test]
    fn test_sniff_skyline_library_file() {
        let file_path = fixture_path();
        let result = sniff_skyline_library_file(&file_path);
        assert!(
            result.is_ok(),
            "File should be detected as Skyline library: {:?}",
            result.err()
        );

        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let non_skyline = PathBuf::from(manifest_dir).join("Cargo.toml");
        assert!(
            sniff_skyline_library_file(non_skyline).is_err(),
            "Cargo.toml should not be detected as Skyline library"
        );
    }

    #[test]
    fn test_sniff_diann_not_skyline() {
        let diann_file = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("diann_io_files")
            .join("sample_lib.txt");
        let result = sniff_skyline_library_file(diann_file);
        assert!(
            result.is_err(),
            "DIA-NN library should not be detected as Skyline library"
        );
    }

    #[test]
    fn test_read_library_file() {
        let elution_groups =
            read_library_file(fixture_path()).expect("Failed to read Skyline library");

        // Fixture has 14 distinct PRTC peptides
        assert_eq!(elution_groups.len(), 14, "Expected 14 elution groups");

        for (eg, extras) in &elution_groups {
            assert_eq!(extras.protein_id, "PRTC");
            assert!(!extras.is_decoy);
            assert_eq!(extras.modified_peptide, extras.stripped_peptide);
            assert!(eg.fragment_count() > 0, "Each precursor should have fragments");
            assert_eq!(
                extras.relative_intensities.len(),
                eg.fragment_count(),
                "Relative intensities should match fragment count"
            );
        }
    }

    #[test]
    fn test_precursor_isotope_rows_are_skipped() {
        let elution_groups = read_library_file(fixture_path()).expect("Failed to read library");

        // SSAAPPPPPR has 3 precursor rows + 4 y fragments in the fixture
        let ssaa = elution_groups
            .iter()
            .find(|(_, extras)| extras.stripped_peptide == "SSAAPPPPPR")
            .expect("SSAAPPPPPR should be present");
        assert_eq!(
            ssaa.0.fragment_count(),
            4,
            "SSAAPPPPPR should have 4 product-ion fragments (precursor rows skipped)"
        );
    }

    #[test]
    fn test_library_intensity_na_handling() {
        // All intensities in the fixture are #N/A -> default to 1.0
        let elution_groups = read_library_file(fixture_path()).expect("Failed to read library");
        for (_, extras) in &elution_groups {
            for (_lab, intensity) in &extras.relative_intensities {
                assert!(
                    (*intensity - 1.0).abs() < f32::EPSILON,
                    "Expected default intensity 1.0, got {}",
                    intensity
                );
            }
        }
    }
}
