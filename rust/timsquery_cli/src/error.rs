use std::path::PathBuf;
use thiserror::Error;
use timsrust::TimsTofPathError;
use timsrust::readers::FrameReaderError;

#[derive(Error, Debug)]
pub enum CliError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("JSON parsing error: {0}")]
    Json(#[from] serde_json::Error),

    #[error("Failed to load TIMS file at '{path}': {source}")]
    TimsFileLoad {
        path: PathBuf,
        source: TimsTofPathError,
    },

    #[error("TIMS frame reader error: {0}")]
    FrameReader(#[from] FrameReaderError),

    #[error("No non-zero intensities found for chromatogram id {0}")]
    EmptyChromatogram(u64),

    #[error("Data processing error: {0}")]
    DataProcessing(String),

    #[error("Data reading error: {0}")]
    DataReading(String),

    #[error("Non-recoverable internal error: {0}")]
    NonRecoverableError(String),
}
