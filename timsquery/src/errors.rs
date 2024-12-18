use crate::models::frames::expanded_frame::FrameProcessingConfig;
use std::fmt::Display;
use timsrust::TimsRustError;

#[derive(Debug)]
pub enum TimsqueryError {
    DataReadingError(DataReadingError),
    DataProcessingError(DataProcessingError),
    Other(String),
}

pub type Result<T> = std::result::Result<T, TimsqueryError>;

impl Display for TimsqueryError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl TimsqueryError {
    pub fn custom(msg: impl Display) -> Self {
        Self::Other(msg.to_string())
    }
}

#[derive(Debug)]
pub enum DataReadingError {
    UnsupportedDataError(UnsupportedDataError),
    TimsRustError(TimsRustError), // Why doesnt timsrust error derive clone?
}

impl From<UnsupportedDataError> for DataReadingError {
    fn from(e: UnsupportedDataError) -> Self {
        DataReadingError::UnsupportedDataError(e)
    }
}

#[derive(Debug)]
pub enum UnsupportedDataError {
    NoMS2DataError,
}

#[derive(Debug)]
pub enum DataProcessingError {
    CentroidingError(FrameProcessingConfig),
    ExpectedVectorLength { real: usize, expected: usize },
    ExpectedNonEmptyData,
    InsufficientData { real: usize, expected: usize },
    ExpectedVectorSameLength(usize, usize),
    UnexpectedInfiniteError(usize),
    UnexpectedInfiniteErrors(Vec<(usize, f64)>),
}

impl From<DataProcessingError> for TimsqueryError {
    fn from(e: DataProcessingError) -> Self {
        TimsqueryError::DataProcessingError(e)
    }
}

impl<T: Into<DataReadingError>> From<T> for TimsqueryError {
    fn from(e: T) -> Self {
        TimsqueryError::DataReadingError(e.into())
    }
}

impl<T: Into<TimsRustError>> From<T> for DataReadingError {
    fn from(e: T) -> Self {
        DataReadingError::TimsRustError(e.into())
    }
}
