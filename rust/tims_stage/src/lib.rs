//! S3-aware input resolution and staging for TimsTOF pipelines.
//!
//! This crate implements Layers 1 (resolution) and 2 (staging policy) of the
//! design at `docs/superpowers/specs/2026-04-21-s3-staging-design.md`. Layer 0
//! primitives live in `timscentroid::storage`. Layer 3 orchestration lives in
//! `timsquery`.

pub mod backend;
pub mod error;
pub mod open;
pub mod prefix;
pub mod resolve;
pub mod tar;
pub mod upload;
pub mod uri;

pub use error::StageError;
