//! Generic URI → sync `Read` for payloads that fit in RAM.
//!
//! Remote URIs are fetched whole via `StorageProvider::get_bytes` and wrapped
//! in `std::io::Cursor<Bytes>`. Callers see a boxed `Read`. Local paths get a
//! `BufReader<File>`. Rationale: the `object_store` GetResult body stream is
//! async; bridging it to a sync `Read` via `block_on_or_in_place` on every
//! `read()` call is reentrancy-hostile under rayon, so we front-load into a
//! single buffer.
//!
//! **This path is only appropriate for small inputs** (speclibs, FASTAs,
//! peptide lists — typically MBs to low-GB). The default cap is
//! [`DEFAULT_IN_MEMORY_CAP`] (4 GiB); remote payloads that exceed this abort
//! with [`StageError::PayloadTooLarge`] *before* the full GET. For streaming
//! reads of larger objects (e.g. raw `.d` bundles), use the tar / prefix
//! staging path instead.
//!
//! Callers that know their payload is bounded by some stricter limit can call
//! [`open_reader_with_cap`] directly.
//!
//! Local files are not size-checked: they do not count against the remote
//! cap because `File::open` is cheap and callers already control local disk.

use crate::common::transport_err;
use crate::error::StageError;
use crate::uri::{
    is_remote_uri,
    split_uri,
};
use bytes::Bytes;
use std::io::{
    Cursor,
    Read,
};
use timscentroid::StorageProvider;

/// Default ceiling on the size of a remote payload loaded into RAM via
/// [`open_reader`]. Intended to catch "someone pointed this at a raw .d" bugs
/// early, not to be a production resource limit.
pub const DEFAULT_IN_MEMORY_CAP: u64 = 4 * 1024 * 1024 * 1024;

pub fn open_reader(uri: &str) -> Result<Box<dyn Read + Send>, StageError> {
    open_reader_with_cap(uri, DEFAULT_IN_MEMORY_CAP)
}

/// Cheap existence probe for a URI. Remote URIs dispatch a HEAD through the
/// provider; local paths call `Path::exists`. Network failures surface as
/// [`StageError::Transport`] — only `NotFound` translates to `Ok(false)`.
pub fn uri_exists(uri: &str) -> Result<bool, StageError> {
    if is_remote_uri(uri) {
        let (loc, key) = split_uri(uri)?;
        let provider = StorageProvider::open(loc).map_err(transport_err(uri))?;
        provider.exists(&key).map_err(transport_err(uri))
    } else {
        Ok(std::path::Path::new(uri).exists())
    }
}

pub fn open_reader_with_cap(uri: &str, max_bytes: u64) -> Result<Box<dyn Read + Send>, StageError> {
    if is_remote_uri(uri) {
        let bytes = fetch_remote_bytes(uri, max_bytes)?;
        Ok(Box::new(Cursor::new(bytes)))
    } else {
        let f = std::fs::File::open(uri).map_err(StageError::Io)?;
        Ok(Box::new(std::io::BufReader::new(f)))
    }
}

fn fetch_remote_bytes(uri: &str, max_bytes: u64) -> Result<Bytes, StageError> {
    let (loc, key) = split_uri(uri)?;
    let provider = StorageProvider::open(loc).map_err(transport_err(uri))?;
    let meta = provider.head(&key).map_err(transport_err(uri))?;
    let size = meta.size;
    if size > max_bytes {
        return Err(StageError::PayloadTooLarge {
            uri: crate::error::redact_uri(uri).to_string(),
            size,
            cap: max_bytes,
        });
    }
    provider.get_bytes(&key).map_err(transport_err(uri))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::TempDir;

    #[test]
    fn reads_local_file() {
        let dir = TempDir::new().unwrap();
        let p = dir.path().join("a.txt");
        std::fs::write(&p, b"hello").unwrap();
        let mut r = open_reader(p.to_str().unwrap()).unwrap();
        let mut s = String::new();
        r.read_to_string(&mut s).unwrap();
        assert_eq!(s, "hello");
    }

    #[test]
    fn local_file_bypasses_cap() {
        // Local files are not size-checked; pick a trivially small cap and
        // show it doesn't matter for local paths.
        let dir = TempDir::new().unwrap();
        let p = dir.path().join("a.txt");
        std::fs::write(&p, b"hello").unwrap();
        let mut r = open_reader_with_cap(p.to_str().unwrap(), 1).unwrap();
        let mut s = String::new();
        r.read_to_string(&mut s).unwrap();
        assert_eq!(s, "hello");
    }

    #[test]
    fn errors_for_missing_local_file() {
        let result = open_reader("/definitely/not/a/real/file.txt");
        assert!(result.is_err());
        if let Err(StageError::Io(_)) = result {
            // success
        } else {
            panic!("expected StageError::Io");
        }
    }
}
