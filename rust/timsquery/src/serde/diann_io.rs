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
    info,
    warn,
};

#[derive(Debug)]
pub enum DiannReadingError {
    IoError(std::io::Error),
    CsvError(csv::Error),
    DiannPrecursorParsingError(DiannPrecursorParsingError),
    UnableToParseElutionGroups,
}

#[derive(Debug)]
pub enum DiannPrecursorParsingError {
    IonParsingError(IonParsingError),
    IonAnnotParseError(String),
    Other(String),
}

impl From<IonParsingError> for DiannPrecursorParsingError {
    fn from(err: IonParsingError) -> Self {
        DiannPrecursorParsingError::IonParsingError(err)
    }
}

impl From<DiannPrecursorParsingError> for DiannReadingError {
    fn from(err: DiannPrecursorParsingError) -> Self {
        DiannReadingError::DiannPrecursorParsingError(err)
    }
}

impl From<csv::Error> for DiannReadingError {
    fn from(err: csv::Error) -> Self {
        DiannReadingError::CsvError(err)
    }
}

impl From<std::io::Error> for DiannReadingError {
    fn from(err: std::io::Error) -> Self {
        DiannReadingError::IoError(err)
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
        return Err(DiannPrecursorParsingError::Other(
            "Empty precursor group".to_string(),
        ));
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
                DiannPrecursorParsingError::Other(format!(
                    "Failed to convert PrecursorCharge to u8: {:?}",
                    e
                ))
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
                    DiannPrecursorParsingError::Other(format!(
                        "Failed to convert FragmentCharge to u8: {:?}",
                        e
                    ))
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

        let frag_char =
            row.fragment_type.chars().next().ok_or_else(|| {
                DiannPrecursorParsingError::Other("Empty fragment type".to_string())
            })?;

        let frag_num = row.fragment_number.try_into().map_err(|_| {
            DiannPrecursorParsingError::Other(format!(
                "Invalid fragment number: {}",
                row.fragment_number
            ))
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
            return Err(DiannPrecursorParsingError::Other(format!(
                "Unexpected Decoy value: {}",
                other
            )));
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
        assert!(is_diann, "speclib2 TSV file should be detected as DIA-NN library");

        // This test mainly checks that the name aliases wotk correctly
        let _elution_groups = read_library_file(file_path).expect("Failed to read library");
    }
}
