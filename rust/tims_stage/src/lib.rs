//! S3-aware input resolution and staging for TimsTOF pipelines.
//!
//! This crate implements Layers 1 (resolution) and 2 (staging policy) of the
//! design at `docs/superpowers/specs/2026-04-21-s3-staging-design.md`. Layer 0
//! primitives live in `timscentroid::storage`. Layer 3 orchestration lives in
//! `timsquery`.

pub mod backend;
pub(crate) mod common;
pub mod download;
pub mod error;
pub mod open;
pub mod prefix;
pub mod resolve;
pub mod tar;
pub mod upload;
pub mod uri;

pub use backend::{
    PerRunTempdir,
    StagedDotD,
    StagingBackend,
    StagingConfig,
};
pub use download::download_to_file;
pub use error::StageError;
pub use open::{
    DEFAULT_IN_MEMORY_CAP,
    open_reader,
    open_reader_with_cap,
    uri_exists,
};
pub use resolve::{
    Resolved,
    SourceSpec,
    resolve,
};
pub use upload::{
    DEFAULT_UPLOAD_CAP,
    upload_file,
    upload_file_with_cap,
};
pub use uri::{
    canonical_uri,
    is_remote_uri,
    sidecar_of,
    split_uri,
};
