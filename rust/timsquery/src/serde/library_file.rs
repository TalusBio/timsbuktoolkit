use super::diann_io::{
    read_library_file as read_diann,
    sniff_diann_library_file,
};
use super::elution_group_inputs::{
    ElutionGroupInput,
    ElutionGroupInputError,
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
pub enum ElutionGroupCollection {
    StringLabels(Vec<TimsElutionGroup<String>>),
    MzpafLabels(Vec<TimsElutionGroup<IonAnnot>>),
    TinyIntLabels(Vec<TimsElutionGroup<u8>>),
    IntLabels(Vec<TimsElutionGroup<u32>>),
}

impl ElutionGroupCollection {
    pub fn len(&self) -> usize {
        match self {
            ElutionGroupCollection::StringLabels(egs) => egs.len(),
            ElutionGroupCollection::MzpafLabels(egs) => egs.len(),
            ElutionGroupCollection::TinyIntLabels(egs) => egs.len(),
            ElutionGroupCollection::IntLabels(egs) => egs.len(),
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
            return Ok(ElutionGroupCollection::TinyIntLabels(out?));
        }
        debug!("Attempting to deserialize elution group inputs with int labels");
        if let Ok(eg_inputs) = serde_json::from_str::<Vec<ElutionGroupInput<u32>>>(content) {
            let out: Result<Vec<TimsElutionGroup<u32>>, ElutionGroupInputError> =
                eg_inputs.into_iter().map(|x| x.try_into()).collect();
            return Ok(ElutionGroupCollection::IntLabels(out?));
        }
        debug!("Attempting to deserialize elution group inputs with mzpaf labels");
        if let Ok(eg_inputs) = serde_json::from_str::<Vec<ElutionGroupInput<IonAnnot>>>(content) {
            let out: Result<Vec<TimsElutionGroup<IonAnnot>>, ElutionGroupInputError> =
                eg_inputs.into_iter().map(|x| x.try_into()).collect();
            return Ok(ElutionGroupCollection::MzpafLabels(out?));
        }
        debug!("Attempting to deserialize elution group inputs with string labels");
        if let Ok(eg_inputs) = serde_json::from_str::<Vec<ElutionGroupInput<String>>>(content) {
            let out: Result<Vec<TimsElutionGroup<String>>, ElutionGroupInputError> =
                eg_inputs.into_iter().map(|x| x.try_into()).collect();
            return Ok(ElutionGroupCollection::StringLabels(out?));
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
            return Ok(ElutionGroupCollection::TinyIntLabels(egs));
        }
        debug!("Attempting to deserialize elution groups with int labels");
        if let Ok(egs) = serde_json::from_str::<Vec<TimsElutionGroup<u32>>>(content) {
            return Ok(ElutionGroupCollection::IntLabels(egs));
        }
        debug!("Attempting to deserialize elution groups with mzpaf labels");
        if let Ok(egs) = serde_json::from_str::<Vec<TimsElutionGroup<IonAnnot>>>(content) {
            return Ok(ElutionGroupCollection::MzpafLabels(egs));
        }
        debug!("Attempting to deserialize elution groups with string labels");
        if let Ok(egs) = serde_json::from_str::<Vec<TimsElutionGroup<String>>>(content) {
            return Ok(ElutionGroupCollection::StringLabels(egs));
        }
        Err(LibraryReadingError::UnableToParseElutionGroups)
    }

    fn try_read_diann(path: &Path) -> Result<Self, LibraryReadingError> {
        let is_diann = sniff_diann_library_file(path);
        info!("Detected DIA-NN library file with label type",);
        if is_diann {
            let egs = match read_diann(path) {
                Ok(egs) => egs,
                Err(e) => {
                    warn!("Failed to read DIA-NN library file: {:?}", e);
                    return Err(LibraryReadingError::UnableToParseElutionGroups);
                }
            };
            let out = egs.into_iter().map(|x| x.0).collect();
            info!("Successfully read DIA-NN library file");
            return Ok(ElutionGroupCollection::MzpafLabels(out));
        }
        Err(LibraryReadingError::UnableToParseElutionGroups)
    }
}

pub fn read_library_file<T: AsRef<Path>>(
    path: T,
) -> Result<ElutionGroupCollection, LibraryReadingError> {
    let diann_attempt = ElutionGroupCollection::try_read_diann(path.as_ref());
    if let Ok(egs) = diann_attempt {
        return Ok(egs);
    }
    ElutionGroupCollection::try_read_json(path.as_ref())
}
