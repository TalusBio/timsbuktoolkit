use serde_json;
use std::path::PathBuf;
use timsquery::{
    DataProcessingError as TQDataProcessingError,
    TimsqueryError,
};
use timsrust::{
    TimsRustError,
    TimsTofPathError,
};

// TODO: break up ... the type system RN gives no info bc
// everything is a datada processing error...
#[derive(Debug)]
pub enum DataProcessingError {
    ExpectedSlicesSameLength {
        expected: usize,
        other: usize,
        context: String,
    },
    ExpectedNonEmptyData {
        context: Option<String>,
    },
    ExpectedFiniteNonNanData {
        context: String,
    },
    ExpectedSetField {
        field: &'static str,
        context: String,
    },
    TimsQueryDataProcessingError {
        error: TQDataProcessingError,
        context: String,
    },
}

#[derive(Debug)]
pub enum LibraryReadingError {
    SpeclibParsingError {
        source: serde_json::Error,
        context: &'static str,
    },
    FileReadingError {
        source: std::io::Error,
        context: &'static str,
        path: PathBuf,
    },
}

impl From<TQDataProcessingError> for DataProcessingError {
    fn from(x: TQDataProcessingError) -> Self {
        Self::TimsQueryDataProcessingError {
            error: x,
            context: "".to_string(),
        }
    }
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
            } => match owned_context {
                Some(x) => x.push_str(context),
                None => *owned_context = Some(context.to_string()),
            },
            DataProcessingError::ExpectedFiniteNonNanData {
                context: owned_context,
            } => {
                owned_context.push_str(context);
            }
            DataProcessingError::ExpectedSetField {
                field: _owned_field,
                context: owned_context,
            } => {
                owned_context.push_str(context);
            }
            DataProcessingError::TimsQueryDataProcessingError {
                context: owned_context,
                ..
            } => {
                owned_context.push_str(context);
            }
        }
        self
    }
}

#[derive(Debug)]
pub enum TimsSeekError {
    TimsTofPath(TimsTofPathError),
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
    LibraryReadingError(LibraryReadingError),
}

impl std::fmt::Display for TimsSeekError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

pub type Result<T> = std::result::Result<T, TimsSeekError>;

impl From<TimsTofPathError> for TimsSeekError {
    fn from(x: TimsTofPathError) -> Self {
        Self::TimsTofPath(x)
    }
}

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

impl From<serde_json::Error> for TimsSeekError {
    fn from(val: serde_json::Error) -> Self {
        TimsSeekError::ParseError {
            msg: val.to_string(),
        }
    }
}

impl From<DataProcessingError> for TimsSeekError {
    fn from(x: DataProcessingError) -> Self {
        Self::DataProcessingError(x)
    }
}

impl From<LibraryReadingError> for TimsSeekError {
    fn from(x: LibraryReadingError) -> Self {
        Self::LibraryReadingError(x)
    }
}

impl From<TQDataProcessingError> for TimsSeekError {
    fn from(x: TQDataProcessingError) -> Self {
        Self::DataProcessingError(DataProcessingError::TimsQueryDataProcessingError {
            error: x,
            context: "".to_string(),
        })
    }
}
