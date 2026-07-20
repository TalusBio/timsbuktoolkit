//! Decision-tree resolution from URI string to Resolved branch.
//!
//! See `docs/superpowers/specs/2026-04-21-s3-staging-design.md` §"Resolution
//! decision tree".

use crate::error::{
    StageError,
    redact_uri,
};
use crate::uri::{
    LocKind,
    NameKind,
    parse_uri_shape,
    sidecar_of,
};
use std::path::PathBuf;
use timscentroid::{
    StorageLocation,
    StorageProvider,
};

#[derive(Debug)]
pub enum Resolved {
    /// A prebuilt `.idx` (local or remote) — load directly, no reader dispatch.
    Idx { loc: StorageLocation },
    /// A raw vendor artifact (local path or `s3://…`). Handed to
    /// [`crate::load::load_raw`], which sniffs the reader and manifest-stages
    /// remote inputs. Carries the canonical URI verbatim.
    Raw { uri: String },
    /// A tarred `.d` container — extracted to a tempdir first, then read as a
    /// local `.d`. (Bruker-specific transport.)
    Tar { spec: SourceSpec },
}

#[derive(Debug, Clone)]
pub enum SourceSpec {
    S3Tar { loc: StorageLocation, key: String },
    LocalTar { path: PathBuf },
}

/// Resolve a URI to a concrete Resolved variant. May perform at most one
/// remote round-trip: a HEAD for `<sidecar_uri>/metadata.json`. No staging
/// or downloading happens here.
pub fn resolve(uri: &str) -> Result<Resolved, StageError> {
    let shape = parse_uri_shape(uri)?;

    // If user already pointed at an .idx, no sidecar hunt.
    if shape.name == NameKind::Idx {
        return Ok(Resolved::Idx {
            loc: idx_location_for(uri)?,
        });
    }

    // Check for an .idx sidecar. The `.idx` is a DIRECTORY containing
    // `metadata.json` + parquet shards; probe by existence of that metadata
    // file inside the directory.
    let sidecar = sidecar_of(uri);
    if sidecar_ready(&sidecar)? {
        return Ok(Resolved::Idx {
            loc: idx_location_for(&sidecar)?,
        });
    }

    // No sidecar. `.tar` is a container (extract first); everything else is a
    // raw artifact handed to the reader registry via `load_raw` — no vendor
    // shape decided here.
    match (shape.loc, shape.name) {
        (_, NameKind::Raw) => Ok(Resolved::Raw {
            uri: uri.trim_end_matches('/').to_string(),
        }),
        (LocKind::Local, NameKind::Tar) => Ok(Resolved::Tar {
            spec: SourceSpec::LocalTar {
                path: PathBuf::from(uri),
            },
        }),
        (LocKind::Remote, NameKind::Tar) => {
            let (loc, key) = crate::uri::split_uri(uri)?;
            Ok(Resolved::Tar {
                spec: SourceSpec::S3Tar { loc, key },
            })
        }
        (_, NameKind::Idx) => unreachable!("handled above"),
    }
}

/// Build a `StorageLocation` that ROOTS at the `.idx` directory itself.
/// For local paths: `StorageLocation::Local(<path>)`. For remote URIs:
/// `StorageLocation::Url(<uri>)`.
fn idx_location_for(idx_uri: &str) -> Result<StorageLocation, StageError> {
    let trimmed = idx_uri.trim_end_matches('/');
    if crate::uri::is_remote_uri(trimmed) {
        StorageLocation::from_url(trimmed)
            .map_err(|e| StageError::InvalidUri(format!("{idx_uri}: {e}")))
    } else {
        Ok(StorageLocation::from_path(trimmed))
    }
}

/// Probe whether an `.idx` directory is ready by checking for its
/// `metadata.json`. Uses the non-creating `StorageProvider::open` so we
/// don't inadvertently materialize local directories during resolution.
fn sidecar_ready(sidecar_uri: &str) -> Result<bool, StageError> {
    let loc = match idx_location_for(sidecar_uri) {
        Ok(l) => l,
        Err(_) => return Ok(false),
    };
    let provider = match StorageProvider::open(loc) {
        Ok(p) => p,
        Err(_) => return Ok(false),
    };
    provider
        .exists("metadata.json")
        .map_err(|e| StageError::Transport {
            uri: redact_uri(sidecar_uri),
            source: e,
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    /// Build a fake `.idx` directory that passes the sidecar probe: it must
    /// contain a `metadata.json` file. Content doesn't matter for resolve().
    fn fake_idx_dir(root: &std::path::Path, name: &str) -> std::path::PathBuf {
        let p = root.join(name);
        std::fs::create_dir(&p).unwrap();
        std::fs::write(p.join("metadata.json"), b"{}").unwrap();
        p
    }

    #[test]
    fn local_idx_resolves_directly() {
        let dir = TempDir::new().unwrap();
        let idx = fake_idx_dir(dir.path(), "sample.d.idx");
        let r = resolve(idx.to_str().unwrap()).unwrap();
        assert!(matches!(r, Resolved::Idx { .. }));
    }

    #[test]
    fn local_raw_with_sidecar_becomes_idx() {
        let dir = TempDir::new().unwrap();
        let dotd = dir.path().join("sample.d");
        std::fs::create_dir(&dotd).unwrap();
        fake_idx_dir(dir.path(), "sample.d.idx");
        let r = resolve(dotd.to_str().unwrap()).unwrap();
        assert!(matches!(r, Resolved::Idx { .. }));
    }

    #[test]
    fn local_raw_without_sidecar_becomes_raw() {
        let dir = TempDir::new().unwrap();
        let dotd = dir.path().join("sample.d");
        std::fs::create_dir(&dotd).unwrap();
        let r = resolve(dotd.to_str().unwrap()).unwrap();
        assert!(matches!(r, Resolved::Raw { .. }));
    }

    #[test]
    fn local_tar_without_sidecar_becomes_tar() {
        let dir = TempDir::new().unwrap();
        let tar = dir.path().join("sample.d.tar");
        std::fs::write(&tar, b"fake").unwrap();
        let r = resolve(tar.to_str().unwrap()).unwrap();
        assert!(matches!(
            r,
            Resolved::Tar {
                spec: SourceSpec::LocalTar { .. }
            }
        ));
    }

    #[test]
    fn local_dotd_with_empty_idx_dir_does_not_shortcut() {
        // An .idx dir with no metadata.json is NOT a valid sidecar.
        let dir = TempDir::new().unwrap();
        let dotd = dir.path().join("sample.d");
        std::fs::create_dir(&dotd).unwrap();
        std::fs::create_dir(dir.path().join("sample.d.idx")).unwrap();
        let r = resolve(dotd.to_str().unwrap()).unwrap();
        assert!(
            matches!(r, Resolved::Raw { .. }),
            "empty .idx dir should not shortcut; got {:?}",
            r
        );
    }
}
