use serde_json;
use timsquery::TimsqueryError;
use timsrust::TimsRustError;

#[derive(Debug)]
pub enum DataProcessingError {
    ExpectedSlicesSameLength {
        expected: usize,
        other: usize,
        context: String,
    },
    ExpectedNonEmptyData {
        context: String,
    },
    ExpectedFiniteNonNanData {
        context: String,
    },
}

impl DataProcessingError {
    pub fn append_to_context(mut self, context: &str) -> Self {
        match &mut self {
            DataProcessingError::ExpectedSlicesSameLength {
                context: owned_context,
                ..
            } => {
                owned_context.push_str(context);
            }
            DataProcessingError::ExpectedNonEmptyData {
                context: owned_context,
            } => {
                owned_context.push_str(context);
            }
            DataProcessingError::ExpectedFiniteNonNanData {
                context: owned_context,
            } => {
                owned_context.push_str(context);
            }
        }
        self
    }
}

#[derive(Debug)]
pub enum TimsSeekError {
    TimsRust(TimsRustError),
    Timsquery(TimsqueryError),
    Io {
        source: std::io::Error,
        path: Option<std::path::PathBuf>,
    },
    ParseError {
        msg: String,
    },
    DataProcessingError(DataProcessingError),
}

impl std::fmt::Display for TimsSeekError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

pub type Result<T> = std::result::Result<T, TimsSeekError>;

impl From<TimsRustError> for TimsSeekError {
    fn from(x: TimsRustError) -> Self {
        Self::TimsRust(x)
    }
}

impl From<TimsqueryError> for TimsSeekError {
    fn from(x: TimsqueryError) -> Self {
        Self::Timsquery(x)
    }
}

impl From<std::num::ParseIntError> for TimsSeekError {
    fn from(x: std::num::ParseIntError) -> Self {
        Self::ParseError { msg: x.to_string() }
    }
}

impl Into<TimsSeekError> for serde_json::Error {
    fn into(self) -> TimsSeekError {
        TimsSeekError::ParseError {
            msg: self.to_string(),
        }
    }
}

impl From<DataProcessingError> for TimsSeekError {
    fn from(x: DataProcessingError) -> Self {
        Self::DataProcessingError(x)
    }
}
