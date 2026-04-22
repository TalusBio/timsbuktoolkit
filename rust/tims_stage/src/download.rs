//! Stream a URI to a local file.
//!
//! For URIs that may exceed RAM (speclib parquet, large FASTAs), this is the
//! right primitive — it streams via `StorageProvider::get_to_file` for remote
//! and is a plain file copy for local. Unlike [`crate::open_reader`] it does
//! not buffer the whole payload in memory, so no size cap is imposed.

use crate::common::transport_err;
use crate::error::StageError;
use crate::uri::{
    is_remote_uri,
    split_uri,
};
use std::path::Path;
use timscentroid::StorageProvider;

pub fn download_to_file(uri: &str, dst: &Path) -> Result<(), StageError> {
    if let Some(parent) = dst.parent() {
        std::fs::create_dir_all(parent).map_err(StageError::Io)?;
    }
    if is_remote_uri(uri) {
        let (loc, key) = split_uri(uri)?;
        let provider = StorageProvider::open(loc).map_err(transport_err(uri))?;
        let bar = indicatif::ProgressBar::hidden();
        provider
            .get_to_file(&key, dst, &bar)
            .map_err(transport_err(uri))?;
    } else {
        std::fs::copy(uri, dst).map_err(StageError::Io)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn copies_local() {
        let dir = TempDir::new().unwrap();
        let src = dir.path().join("src.bin");
        let dst = dir.path().join("nested/dst.bin");
        std::fs::write(&src, b"payload").unwrap();
        download_to_file(src.to_str().unwrap(), &dst).unwrap();
        assert_eq!(std::fs::read(&dst).unwrap(), b"payload");
    }
}
