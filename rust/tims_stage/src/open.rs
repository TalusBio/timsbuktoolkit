//! Generic URI → sync `Read`.
//!
//! Remote URIs are fetched whole via `StorageProvider::get_bytes` and wrapped
//! in `std::io::Cursor<Bytes>`. Callers see a boxed `Read`. Local paths get a
//! `BufReader<File>`. Rationale: the `object_store` GetResult body stream is
//! async; bridging it to a sync `Read` via `block_on_or_in_place` on every
//! `read()` call is reentrancy-hostile under rayon, so we front-load into a
//! single buffer. Speclibs and FASTAs are MBs to a few GB max — fits RAM.

use crate::error::{
    StageError,
    redact_uri,
};
use crate::uri::split_uri;
use bytes::Bytes;
use std::io::{
    Cursor,
    Read,
};
use timscentroid::StorageProvider;

pub fn open_reader(uri: &str) -> Result<Box<dyn Read + Send>, StageError> {
    if is_remote(uri) {
        let bytes = fetch_remote_bytes(uri)?;
        Ok(Box::new(Cursor::new(bytes)))
    } else {
        let f = std::fs::File::open(uri).map_err(StageError::Io)?;
        Ok(Box::new(std::io::BufReader::new(f)))
    }
}

fn is_remote(uri: &str) -> bool {
    uri.starts_with("s3://") || uri.starts_with("gs://") || uri.starts_with("az://")
}

fn fetch_remote_bytes(uri: &str) -> Result<Bytes, StageError> {
    let (loc, key) = split_uri(uri)?;
    let provider = StorageProvider::open(loc).map_err(|e| StageError::Transport {
        uri: redact_uri(uri),
        source: e,
    })?;
    provider.get_bytes(&key).map_err(|e| StageError::Transport {
        uri: redact_uri(uri),
        source: e,
    })
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
