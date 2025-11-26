use std::path::PathBuf;
use thiserror::Error;
use timsrust::TimsTofPathError;
use timsrust::readers::FrameReaderError;

#[derive(Error, Debug)]
pub enum ViewerError {
    #[error("IO error: {0}")]
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

    #[error("Failed to load data from {path}: {source}")]
    DataLoading {
        path: PathBuf,
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },

    #[error("General error: {0}")]
    General(String),
}
