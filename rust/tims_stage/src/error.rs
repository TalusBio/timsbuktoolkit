//! Errors produced by the tims_stage crate.

use thiserror::Error;

#[derive(Debug, Error)]
pub enum StageError {
    #[error("unknown URI suffix: {0} (expected .idx, .d, .d/, or .tar)")]
    UnknownSuffix(String),

    #[error("URI malformed: {0}")]
    InvalidUri(String),

    #[error("prefix '{prefix}' yields more than {cap} entries; refusing to list")]
    PrefixCapExceeded { prefix: String, cap: usize },

    #[error(".d source missing required basenames {missing:?}; found {found:?}")]
    MissingRequiredFiles {
        missing: Vec<String>,
        found: Vec<String>,
    },

    #[error(
        "GNU long-name tar entries are not supported; rebuild with `tar --format=ustar` or upload as an S3 prefix"
    )]
    UnsupportedTarFeature,

    #[error("tar payload short read: expected {expected} bytes, got {actual}")]
    ShortRead { expected: u64, actual: u64 },

    #[error("shape mismatch: {0}")]
    ShapeMismatch(String),

    #[error("storage transport error at {uri}: {source}")]
    Transport {
        uri: String,
        #[source]
        source: timscentroid::serialization::SerializationError,
    },

    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
}

/// Strip query strings from a URI before embedding in an error message —
/// presigned-style URLs carry credentials in the query.
pub(crate) fn redact_uri(uri: &str) -> String {
    match uri.split_once('?') {
        Some((base, _)) => base.to_string(),
        None => uri.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn redact_uri_strips_query_string() {
        assert_eq!(redact_uri("s3://bkt/key?X-Sig=abc"), "s3://bkt/key");
        assert_eq!(redact_uri("/tmp/a.bin"), "/tmp/a.bin");
    }

    #[test]
    fn stage_error_display_has_context() {
        let e = StageError::UnknownSuffix("/foo.txt".into());
        let s = format!("{e}");
        assert!(s.contains("/foo.txt"));
    }
}
