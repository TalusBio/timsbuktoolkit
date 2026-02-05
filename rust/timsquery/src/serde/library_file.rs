pub use super::diann_io::DiannPrecursorExtras;
use super::diann_io::{
    read_library_file as read_diann_tsv,
    read_parquet_library_file as read_diann_parquet,
    sniff_diann_library_file,
    sniff_diann_parquet_library_file,
};
use super::elution_group_inputs::{
    ElutionGroupInput,
    ElutionGroupInputError,
};
pub use super::spectronaut_io::SpectronautPrecursorExtras;
use super::spectronaut_io::{
    read_library_file as read_spectronaut_tsv,
    sniff_spectronaut_library_file,
};
use crate::TimsElutionGroup;
use crate::ion::IonAnnot;
use std::path::Path;
use tracing::{
    debug,
    info,
    warn,
};

#[derive(Debug)]
pub enum LibraryReadingError {
    IoError(std::io::Error),
    SerdeJsonError(serde_json::Error),
    ElutionGroupInputError(ElutionGroupInputError),
    UnableToParseElutionGroups,
}

impl From<serde_json::Error> for LibraryReadingError {
    fn from(err: serde_json::Error) -> Self {
        LibraryReadingError::SerdeJsonError(err)
    }
}
impl From<ElutionGroupInputError> for LibraryReadingError {
    fn from(err: ElutionGroupInputError) -> Self {
        LibraryReadingError::ElutionGroupInputError(err)
    }
}

#[derive(Debug)]
pub enum FileReadingExtras {
    Diann(Vec<DiannPrecursorExtras>),
    Spectronaut(Vec<SpectronautPrecursorExtras>),
}

#[derive(Debug)]
pub enum ElutionGroupCollection {
    StringLabels(Vec<TimsElutionGroup<String>>, Option<FileReadingExtras>),
    MzpafLabels(Vec<TimsElutionGroup<IonAnnot>>, Option<FileReadingExtras>),
    TinyIntLabels(Vec<TimsElutionGroup<u8>>, Option<FileReadingExtras>),
    IntLabels(Vec<TimsElutionGroup<u32>>, Option<FileReadingExtras>),
}

impl ElutionGroupCollection {
    pub fn len(&self) -> usize {
        match self {
            ElutionGroupCollection::StringLabels(egs, _) => egs.len(),
            ElutionGroupCollection::MzpafLabels(egs, _) => egs.len(),
            ElutionGroupCollection::TinyIntLabels(egs, _) => egs.len(),
            ElutionGroupCollection::IntLabels(egs, _) => egs.len(),
        }
    }

    fn try_read_json(path: &Path) -> Result<Self, LibraryReadingError> {
        let file_content = std::fs::read_to_string(path).map_err(LibraryReadingError::IoError)?;
        info!("Read file content from {}", path.display());
        // First try direct deserialization
        if let Ok(egs) = Self::try_deser_direct(&file_content) {
            info!("Successfully deserialized elution groups directly");
            return Ok(egs);
        }
        // Next try deserialization via inputed format
        match Self::try_deser_inputed(&file_content) {
            Ok(egs) => {
                info!("Successfully deserialized elution groups via inputed format");
                Ok(egs)
            }
            Err(_) => Err(LibraryReadingError::UnableToParseElutionGroups),
        }
    }

    fn try_deser_inputed(content: &str) -> Result<Self, LibraryReadingError> {
        // We can try from smallest to largest overhead
        // Here we try to do the deser into ElutionGroupInput variants first
        debug!("Attempting deserialization of elution group inputs");
        debug!("Attempting to deserialize elution group inputs with tiny int labels");
        if let Ok(eg_inputs) = serde_json::from_str::<Vec<ElutionGroupInput<u8>>>(content) {
            // Here we can handle filling the inputs if they are needed...
            let eg_inputs = if eg_inputs.first().is_some_and(|x| x.needs_fragment_labels()) {
                debug!("Filling missing fragment labels with tiny int labels");
                eg_inputs
                    .into_iter()
                    .map(|x| x.try_fill_labels_u8())
                    .collect::<Result<_, _>>()
            } else {
                Ok(eg_inputs)
            };

            let out: Result<Vec<TimsElutionGroup<u8>>, ElutionGroupInputError> = eg_inputs?
                .into_iter()
                .map(<ElutionGroupInput<u8> as TryInto<TimsElutionGroup<u8>>>::try_into)
                .collect();
            return Ok(ElutionGroupCollection::TinyIntLabels(out?, None));
        }
        debug!("Attempting to deserialize elution group inputs with int labels");
        if let Ok(eg_inputs) = serde_json::from_str::<Vec<ElutionGroupInput<u32>>>(content) {
            let out: Result<Vec<TimsElutionGroup<u32>>, ElutionGroupInputError> =
                eg_inputs.into_iter().map(|x| x.try_into()).collect();
            return Ok(ElutionGroupCollection::IntLabels(out?, None));
        }
        debug!("Attempting to deserialize elution group inputs with mzpaf labels");
        if let Ok(eg_inputs) = serde_json::from_str::<Vec<ElutionGroupInput<IonAnnot>>>(content) {
            let out: Result<Vec<TimsElutionGroup<IonAnnot>>, ElutionGroupInputError> =
                eg_inputs.into_iter().map(|x| x.try_into()).collect();
            return Ok(ElutionGroupCollection::MzpafLabels(out?, None));
        }
        debug!("Attempting to deserialize elution group inputs with string labels");
        if let Ok(eg_inputs) = serde_json::from_str::<Vec<ElutionGroupInput<String>>>(content) {
            let out: Result<Vec<TimsElutionGroup<String>>, ElutionGroupInputError> =
                eg_inputs.into_iter().map(|x| x.try_into()).collect();
            return Ok(ElutionGroupCollection::StringLabels(out?, None));
        }
        Err(LibraryReadingError::UnableToParseElutionGroups)
    }

    fn try_deser_direct(content: &str) -> Result<Self, LibraryReadingError> {
        // We can try from smallest to largest overhead
        // Here we try to do the direct deser into ElutionGroupCollection variants
        // u8 -> u32 -> IonAnnot -> String
        debug!("Attempting direct deserialization of elution groups");
        debug!("Attempting to deserialize elution groups with tiny int labels");
        if let Ok(egs) = serde_json::from_str::<Vec<TimsElutionGroup<u8>>>(content) {
            return Ok(ElutionGroupCollection::TinyIntLabels(egs, None));
        }
        debug!("Attempting to deserialize elution groups with int labels");
        if let Ok(egs) = serde_json::from_str::<Vec<TimsElutionGroup<u32>>>(content) {
            return Ok(ElutionGroupCollection::IntLabels(egs, None));
        }
        debug!("Attempting to deserialize elution groups with mzpaf labels");
        if let Ok(egs) = serde_json::from_str::<Vec<TimsElutionGroup<IonAnnot>>>(content) {
            return Ok(ElutionGroupCollection::MzpafLabels(egs, None));
        }
        debug!("Attempting to deserialize elution groups with string labels");
        if let Ok(egs) = serde_json::from_str::<Vec<TimsElutionGroup<String>>>(content) {
            return Ok(ElutionGroupCollection::StringLabels(egs, None));
        }
        Err(LibraryReadingError::UnableToParseElutionGroups)
    }

    fn try_read_diann(path: &Path) -> Result<Self, LibraryReadingError> {
        // First check if it's a parquet file
        let is_diann_parquet = sniff_diann_parquet_library_file(path);
        if is_diann_parquet {
            info!("Detected DIA-NN parquet library file (DiaNN 2.2+)");
            let egs = match read_diann_parquet(path) {
                Ok(egs) => egs,
                Err(e) => {
                    warn!("Failed to read DIA-NN parquet library file: {:?}", e);
                    return Err(LibraryReadingError::UnableToParseElutionGroups);
                }
            };
            let (egs, extras): (Vec<_>, Vec<_>) = egs.into_iter().unzip();
            info!("Successfully read DIA-NN parquet library file");
            return Ok(ElutionGroupCollection::MzpafLabels(
                egs,
                Some(FileReadingExtras::Diann(extras)),
            ));
        }

        // Then check if it's a TSV file
        let is_diann_tsv = sniff_diann_library_file(path);
        if is_diann_tsv {
            info!("Detected DIA-NN TSV library file");
            let egs = match read_diann_tsv(path) {
                Ok(egs) => egs,
                Err(e) => {
                    warn!("Failed to read DIA-NN TSV library file: {:?}", e);
                    return Err(LibraryReadingError::UnableToParseElutionGroups);
                }
            };
            let (egs, extras): (Vec<_>, Vec<_>) = egs.into_iter().unzip();
            info!("Successfully read DIA-NN TSV library file");
            return Ok(ElutionGroupCollection::MzpafLabels(
                egs,
                Some(FileReadingExtras::Diann(extras)),
            ));
        }

        Err(LibraryReadingError::UnableToParseElutionGroups)
    }

    fn try_read_spectronaut(path: &Path) -> Result<Self, LibraryReadingError> {
        match sniff_spectronaut_library_file(path) {
            Ok(()) => {
                info!("Detected Spectronaut TSV library file");
                let egs = match read_spectronaut_tsv(path) {
                    Ok(egs) => egs,
                    Err(e) => {
                        warn!("Failed to read Spectronaut TSV library file: {:?}", e);
                        return Err(LibraryReadingError::UnableToParseElutionGroups);
                    }
                };
                let (egs, extras): (Vec<_>, Vec<_>) = egs.into_iter().unzip();
                info!("Successfully read Spectronaut TSV library file");
                Ok(ElutionGroupCollection::MzpafLabels(
                    egs,
                    Some(FileReadingExtras::Spectronaut(extras)),
                ))
            }
            Err(e) => {
                debug!("File is not Spectronaut format: {:?}", e);
                Err(LibraryReadingError::UnableToParseElutionGroups)
            }
        }
    }
}

pub fn read_library_file<T: AsRef<Path>>(
    path: T,
) -> Result<ElutionGroupCollection, LibraryReadingError> {
    // Try DIA-NN first
    let diann_attempt = ElutionGroupCollection::try_read_diann(path.as_ref());
    if let Ok(egs) = diann_attempt {
        return Ok(egs);
    }

    // Try Spectronaut next
    let spectronaut_attempt = ElutionGroupCollection::try_read_spectronaut(path.as_ref());
    if let Ok(egs) = spectronaut_attempt {
        return Ok(egs);
    }

    // Fall back to JSON
    ElutionGroupCollection::try_read_json(path.as_ref())
}
