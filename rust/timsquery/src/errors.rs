use array2d::Array2DError;
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
    /// Failure dispatching to / building from a raw-format reader (the registry
    /// + staging). Stringified to stay decoupled from the staging error type.
    RawReadError(String),
}

impl From<UnsupportedDataError> for DataReadingError {
    fn from(e: UnsupportedDataError) -> Self {
        DataReadingError::UnsupportedDataError(e)
    }
}
// No `From<ReadError>` impl: the blanket `From<T: Into<TimsRustError>>` below
// would conflict. Map explicitly with `DataReadingError::RawReadError` instead.

// Note: Can't implement From<SerializationError> due to blanket impl conflict
// Use DataReadingError::SerializationError(e) directly instead

#[derive(Debug)]
pub enum UnsupportedDataError {
    NoMS2DataError,
    CloudRawDataNotSupported { url: String, suggestion: String },
    LazyMaterializationUnsupported,
    CacheDisabled,
    InvalidCacheUrl { url: String },
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
            Self::LazyMaterializationUnsupported => write!(
                f,
                "Lazy -> eager materialization not yet implemented; \
                 load the index eagerly from source instead"
            ),
            Self::CacheDisabled => write!(
                f,
                "Cache is disabled; cannot convert eager index to lazy \
                 without a cache location"
            ),
            Self::InvalidCacheUrl { url } => {
                write!(f, "Cache URL is not a valid StorageLocation: {}", url)
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

impl From<Array2DError> for DataProcessingError {
    fn from(e: Array2DError) -> Self {
        match e {
            Array2DError::EmptyData => DataProcessingError::ExpectedNonEmptyData,
            Array2DError::DimensionMismatch => DataProcessingError::ExpectedVectorSameLength,
            Array2DError::IndexOutOfBounds(n) => DataProcessingError::IndexOutOfBoundsError(n),
        }
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
