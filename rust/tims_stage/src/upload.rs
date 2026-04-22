//! Upload a local file to a destination URI (local path or S3).
//!
//! The whole file is read into memory before the PUT. **Only suitable for
//! small outputs** (result JSONs, parquet shards, run reports). For the
//! default cap see [`DEFAULT_UPLOAD_CAP`]. To stream a large `.d` bundle
//! instead, bundle it into a tar and use the staging pipeline in reverse
//! (not yet implemented — avoid pointing `upload_file` at multi-GB payloads).

use crate::common::transport_err;
use crate::error::StageError;
use crate::uri::{
    is_remote_uri,
    split_uri,
};
use std::path::Path;
use timscentroid::StorageProvider;

/// Default ceiling on the size of a payload uploaded via [`upload_file`].
/// Matches [`crate::open::DEFAULT_IN_MEMORY_CAP`] (4 GiB).
pub const DEFAULT_UPLOAD_CAP: u64 = 4 * 1024 * 1024 * 1024;

pub fn upload_file(local: &Path, dest_uri: &str) -> Result<(), StageError> {
    upload_file_with_cap(local, dest_uri, DEFAULT_UPLOAD_CAP)
}

pub fn upload_file_with_cap(
    local: &Path,
    dest_uri: &str,
    max_bytes: u64,
) -> Result<(), StageError> {
    let size = std::fs::metadata(local).map_err(StageError::Io)?.len();
    if size > max_bytes {
        return Err(StageError::PayloadTooLarge {
            uri: crate::error::redact_uri(dest_uri).to_string(),
            size,
            cap: max_bytes,
        });
    }
    let bytes = std::fs::read(local).map_err(StageError::Io)?;
    if is_remote_uri(dest_uri) {
        let (loc, key) = split_uri(dest_uri)?;
        // Writes go through `new` (creates parent dir for local).
        let provider = StorageProvider::new(loc).map_err(transport_err(dest_uri))?;
        provider
            .write_bytes(&key, bytes)
            .map_err(transport_err(dest_uri))?;
    } else {
        if let Some(parent) = Path::new(dest_uri).parent() {
            std::fs::create_dir_all(parent).map_err(StageError::Io)?;
        }
        std::fs::write(dest_uri, &bytes).map_err(StageError::Io)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn copies_local_to_local() {
        let dir = TempDir::new().unwrap();
        let src = dir.path().join("src.bin");
        let dst = dir.path().join("subdir/dst.bin");
        std::fs::write(&src, b"payload").unwrap();
        upload_file(&src, dst.to_str().unwrap()).unwrap();
        assert_eq!(std::fs::read(&dst).unwrap(), b"payload");
    }

    #[test]
    fn rejects_oversize_payload() {
        let dir = TempDir::new().unwrap();
        let src = dir.path().join("src.bin");
        let dst = dir.path().join("dst.bin");
        std::fs::write(&src, b"1234567890").unwrap();
        let result = upload_file_with_cap(&src, dst.to_str().unwrap(), 5);
        assert!(matches!(result, Err(StageError::PayloadTooLarge { .. })));
    }
}
