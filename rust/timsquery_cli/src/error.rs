use thiserror::Error;
use timsrust::readers::FrameReaderError;

#[derive(Error, Debug)]
pub enum CliError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("JSON parsing error: {0}")]
    Json(#[from] serde_json::Error),

    #[error("TIMS frame reader error: {0}")]
    FrameReader(#[from] FrameReaderError),

    #[error("Data processing error: {0}")]
    DataProcessing(String),

    #[error("Data reading error: {0}")]
    DataReading(String),
}
