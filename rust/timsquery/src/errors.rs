use std::fmt::Display;
use timscentroid::serialization::SerializationError;
use timsrust::{
    TimsRustError,
    TimsTofPathError,
};

#[derive(Debug)]
pub enum TimsqueryError {
    DataReadingError(DataReadingError),
    DataProcessingError(DataProcessingError),
    Other(String),
}

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
    TimsTofPathError(TimsTofPathError),
    TimsRustError(TimsRustError), // Why doesnt timsrust error derive clone?
    SerializationError(SerializationError),
}

impl From<UnsupportedDataError> for DataReadingError {
    fn from(e: UnsupportedDataError) -> Self {
        DataReadingError::UnsupportedDataError(e)
    }
}

// Note: Can't implement From<SerializationError> due to blanket impl conflict
// Use DataReadingError::SerializationError(e) directly instead

#[derive(Debug)]
pub enum UnsupportedDataError {
    NoMS2DataError,
    CloudRawDataNotSupported { url: String, suggestion: String },
}

impl Display for UnsupportedDataError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NoMS2DataError => write!(f, "No MS2 data found"),
            Self::CloudRawDataNotSupported { url, suggestion } => {
                write!(
                    f,
                    "Cannot read raw .d files from cloud storage: {}\n\n{}",
                    url, suggestion
                )
            }
        }
    }
}

#[derive(Debug)]
pub enum DataProcessingError {
    ExpectedVectorLength { real: usize, expected: usize },
    ExpectedNonEmptyData,
    InsufficientData { real: usize, expected: usize },
    ExpectedVectorSameLength,
    IndexOutOfBoundsError(usize),
    UnexpectedInfiniteError(usize),
    UnexpectedInfiniteErrors(Vec<(usize, f64)>),
    KeyNotFound,
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
