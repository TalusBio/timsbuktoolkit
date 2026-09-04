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

    #[error("General error: {0}")]
    General(String),
}

impl From<timsquery::serde::TargetReadingError> for ViewerError {
    fn from(err: timsquery::serde::TargetReadingError) -> Self {
        match err {
            timsquery::serde::TargetReadingError::IoError(e) => ViewerError::Io(e),
            timsquery::serde::TargetReadingError::SerdeJsonError(e) => ViewerError::Json(e),
            timsquery::serde::TargetReadingError::ElutionGroupInputError(e) => {
                ViewerError::General(format!("Elution group input error: {:?}", e))
            }
            timsquery::serde::TargetReadingError::UnableToParseElutionGroups => {
                ViewerError::General("Unable to parse elution groups".to_string())
            }
            timsquery::serde::TargetReadingError::UnsupportedSpeclibVersion(v) => {
                ViewerError::General(format!("Unsupported .speclib version: {v}"))
            }
            timsquery::serde::TargetReadingError::SpeclibParse(msg) => {
                ViewerError::General(format!(".speclib parse error: {msg}"))
            }
            timsquery::serde::TargetReadingError::SourceId(e) => {
                ViewerError::General(format!("source id error: {e}"))
            }
        }
    }
}
