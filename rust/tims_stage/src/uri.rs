//! URI shape parsing + sidecar convention + single canonical splitter.

use crate::error::StageError;
use std::path::Path;
use timscentroid::StorageLocation;

/// Where an URI points — local filesystem or remote object store.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum LocKind {
    Local,
    Remote,
}

/// What kind of artifact the URI names, by suffix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum NameKind {
    Idx,
    DotD,
    Tar,
}

#[derive(Debug, Clone)]
pub(crate) struct UriShape {
    pub loc: LocKind,
    pub name: NameKind,
    #[allow(dead_code)]
    pub raw: String,
}

/// Classify a URI string by scheme + suffix. Pure logic, no I/O.
pub(crate) fn parse_uri_shape(uri: &str) -> Result<UriShape, StageError> {
    let trimmed = uri.trim_end_matches('/');
    let loc = if is_remote_uri(uri) {
        LocKind::Remote
    } else {
        LocKind::Local
    };
    let name = if trimmed.ends_with(".idx") {
        NameKind::Idx
    } else if trimmed.ends_with(".tar") {
        NameKind::Tar
    } else if trimmed.ends_with(".d") {
        NameKind::DotD
    } else {
        return Err(StageError::UnknownSuffix(uri.to_string()));
    };
    Ok(UriShape {
        loc,
        name,
        raw: uri.to_string(),
    })
}

/// Cheap test for a supported remote URI scheme.
/// Exported so consumers don't re-implement `uri.starts_with("s3://") || ...`.
pub fn is_remote_uri(uri: &str) -> bool {
    remote_scheme_prefix(uri).is_some()
}

/// Expand `~/` (and bare `~`) in local-path URIs to `$HOME`. Remote URIs
/// pass through unchanged. `Path::exists()` and `File::open()` don't
/// recognise `~`; callers pulling paths from config files / non-shell
/// sources should run this before any filesystem probe.
pub fn expand_local_uri(uri: &str) -> String {
    if is_remote_uri(uri) {
        return uri.to_string();
    }
    let home = std::env::var_os("HOME");
    match (uri, &home) {
        ("~", Some(h)) => h.to_string_lossy().into_owned(),
        (u, Some(h)) if u.starts_with("~/") => {
            let mut p = std::path::PathBuf::from(h);
            p.push(&u[2..]);
            p.to_string_lossy().into_owned()
        }
        _ => uri.to_string(),
    }
}

fn remote_scheme_prefix(uri: &str) -> Option<&'static str> {
    if uri.starts_with("s3://") {
        Some("s3://")
    } else if uri.starts_with("gs://") {
        Some("gs://")
    } else if uri.starts_with("az://") {
        Some("az://")
    } else {
        None
    }
}

/// Canonical `.idx` sidecar URI. Strips one trailing `/` so `sample.d/` →
/// `sample.d.idx`, not `sample.d/.idx`. Matches the existing index_serde
/// convention for local paths.
pub fn sidecar_of(uri: &str) -> String {
    let trimmed = uri.strip_suffix('/').unwrap_or(uri);
    format!("{trimmed}.idx")
}

/// Split `s3://bucket/path/to/obj` into
/// (`StorageLocation::from_url("s3://bucket")`, `"path/to/obj"`).
/// For local paths, returns
/// (`StorageLocation::from_path(parent)`, `file_name`). One splitter used by
/// all downstream call sites.
pub fn split_uri(uri: &str) -> Result<(StorageLocation, String), StageError> {
    let trimmed = uri.trim_end_matches('/');
    if let Some(scheme_prefix) = remote_scheme_prefix(trimmed) {
        let rest = &trimmed[scheme_prefix.len()..];
        let (bucket, key) = rest.split_once('/').unwrap_or((rest, ""));
        let loc = StorageLocation::from_url(format!("{scheme_prefix}{bucket}"))
            .map_err(|e| StageError::InvalidUri(format!("{uri}: {e}")))?;
        return Ok((loc, key.to_string()));
    }
    let p = Path::new(trimmed);
    let parent = p.parent().unwrap_or_else(|| Path::new("."));
    let tail = p
        .file_name()
        .ok_or_else(|| StageError::InvalidUri(format!("missing filename in {uri}")))?
        .to_string_lossy()
        .to_string();
    Ok((StorageLocation::from_path(parent), tail))
}

/// Canonical URI used as a cache key / sidecar-origin identifier.
///
/// - Strips exactly one trailing `/` if present (so `foo.d/` and `foo.d`
///   yield the same key).
/// - Preserves scheme + host + key verbatim for remote URIs (S3 keys are
///   case-sensitive).
/// - For local paths, runs `std::fs::canonicalize` when the path exists;
///   otherwise returns the trimmed input unchanged.
pub fn canonical_uri(uri: &str) -> String {
    let trimmed = uri.trim_end_matches('/').to_string();
    if is_remote_uri(&trimmed) {
        return trimmed;
    }
    match std::fs::canonicalize(&trimmed) {
        Ok(p) => p.to_string_lossy().to_string(),
        Err(_) => trimmed,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recognizes_local_dotd() {
        let s = parse_uri_shape("/tmp/sample.d").unwrap();
        assert_eq!(s.loc, LocKind::Local);
        assert_eq!(s.name, NameKind::DotD);
    }

    #[test]
    fn recognizes_local_dotd_with_trailing_slash() {
        let s = parse_uri_shape("/tmp/sample.d/").unwrap();
        assert_eq!(s.name, NameKind::DotD);
    }

    #[test]
    fn recognizes_s3_tar() {
        let s = parse_uri_shape("s3://bkt/sample.d.tar").unwrap();
        assert_eq!(s.loc, LocKind::Remote);
        assert_eq!(s.name, NameKind::Tar);
    }

    #[test]
    fn recognizes_s3_idx() {
        let s = parse_uri_shape("s3://bkt/sample.d.idx").unwrap();
        assert_eq!(s.name, NameKind::Idx);
    }

    #[test]
    fn errors_on_unknown_suffix() {
        let err = parse_uri_shape("/tmp/some.txt").unwrap_err();
        assert!(matches!(err, StageError::UnknownSuffix(_)));
    }

    #[test]
    fn sidecar_strips_trailing_slash() {
        assert_eq!(sidecar_of("/tmp/sample.d"), "/tmp/sample.d.idx");
        assert_eq!(sidecar_of("/tmp/sample.d/"), "/tmp/sample.d.idx");
    }

    #[test]
    fn sidecar_of_s3_prefix() {
        assert_eq!(sidecar_of("s3://bkt/sample.d/"), "s3://bkt/sample.d.idx");
    }

    #[test]
    fn sidecar_of_s3_tar() {
        assert_eq!(
            sidecar_of("s3://bkt/sample.d.tar"),
            "s3://bkt/sample.d.tar.idx"
        );
    }

    #[test]
    fn splits_s3_uri() {
        let (_loc, key) = split_uri("s3://bkt/prefix/a.bin").unwrap();
        assert_eq!(key, "prefix/a.bin");
    }

    #[test]
    fn splits_s3_bucket_root() {
        let (_loc, key) = split_uri("s3://bkt/a.bin").unwrap();
        assert_eq!(key, "a.bin");
    }

    #[test]
    fn splits_local_path() {
        let (_loc, key) = split_uri("/tmp/a.bin").unwrap();
        assert_eq!(key, "a.bin");
    }

    #[test]
    fn canonical_uri_strips_trailing_slash_on_remote() {
        assert_eq!(canonical_uri("s3://bkt/foo/"), "s3://bkt/foo");
    }

    #[test]
    fn canonical_uri_local_existing_path_is_absolute() {
        let d = tempfile::TempDir::new().unwrap();
        let p = d.path().join("x");
        std::fs::write(&p, b"").unwrap();
        let got = canonical_uri(p.to_str().unwrap());
        assert!(std::path::Path::new(&got).is_absolute());
    }

    #[test]
    fn canonical_uri_local_nonexistent_returns_trimmed() {
        let got = canonical_uri("/tmp/definitely-does-not-exist-xyz123/");
        assert_eq!(got, "/tmp/definitely-does-not-exist-xyz123");
    }

    #[test]
    fn expand_local_uri_tilde() {
        // SAFETY: test-only env mutation; tests are serial within this module.
        unsafe {
            std::env::set_var("HOME", "/home/test");
        }
        assert_eq!(expand_local_uri("~/data/a.txt"), "/home/test/data/a.txt");
        assert_eq!(expand_local_uri("~"), "/home/test");
        assert_eq!(expand_local_uri("/abs/path"), "/abs/path");
        assert_eq!(expand_local_uri("relative/path"), "relative/path");
        // Remote URIs pass through even if they start with ~ by accident.
        assert_eq!(expand_local_uri("s3://bkt/key"), "s3://bkt/key");
    }
}
