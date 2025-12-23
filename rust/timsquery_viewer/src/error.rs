use std::path::PathBuf;
use thiserror::Error;
use timsrust::readers::FrameReaderError;

#[derive(Error, Debug)]
pub enum ViewerError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("JSON parsing error: {0}")]
    Json(#[from] serde_json::Error),

    #[error("TIMS frame reader error: {0}")]
    FrameReader(#[from] FrameReaderError),

    #[error("Failed to load data from {path}: {source}")]
    DataLoading {
        path: PathBuf,
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },

    #[error("General error: {0}")]
    General(String),
}

impl From<timsquery::serde::LibraryReadingError> for ViewerError {
    fn from(err: timsquery::serde::LibraryReadingError) -> Self {
        match err {
            timsquery::serde::LibraryReadingError::IoError(e) => ViewerError::Io(e),
            timsquery::serde::LibraryReadingError::SerdeJsonError(e) => ViewerError::Json(e),
            timsquery::serde::LibraryReadingError::ElutionGroupInputError(e) => {
                ViewerError::General(format!("Elution group input error: {:?}", e))
            }
            timsquery::serde::LibraryReadingError::UnableToParseElutionGroups => {
                ViewerError::General("Unable to parse elution groups".to_string())
            }
        }
    }
}
