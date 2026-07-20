//! The single raw-index composition: URI → `IndexedTimstofPeaks`.
//!
//! `load_raw` is the one place raw dispatch happens. It sniffs the URI via the
//! reader registry, asks the chosen reader for its [`Manifest`] (the files it
//! needs), materializes those — in place for local inputs, or by fetching
//! exactly the declared files for remote ones — and calls `read`. Transport
//! never guesses vendor shape: it fetches what the reader declared, by name.
//!
//! Lives in `tims_stage` (one crate above `timscentroid`, where the registry
//! lives) so it can compose `pick` + `manifest` + `read` with the fetch engine
//! without a dependency cycle.

use std::ffi::OsString;
use std::path::Path;

use http::Uri;
use timscentroid::reader::{
    Manifest,
    ReadError,
    ReaderRegistry,
    ResolvedSource,
    local_uri,
};
use timscentroid::{
    CentroidingConfig,
    IndexedTimstofPeaks,
    StorageProvider,
};

use crate::backend::StagingBackend;
use crate::common::{
    make_bar,
    transport_err,
};
use crate::error::StageError;
use crate::uri::{
    is_remote_uri,
    split_uri,
};

/// A built index plus provenance: which reader built it, and whether that reader
/// caches to `.idx` (so callers can decide whether to persist a sidecar).
pub struct RawRead {
    pub index: IndexedTimstofPeaks,
    pub reader_name: &'static str,
    pub caches_to_idx: bool,
}

#[derive(Debug, thiserror::Error)]
pub enum LoadRawError {
    #[error(transparent)]
    Read(#[from] ReadError),
    #[error(transparent)]
    Stage(#[from] StageError),
}

/// Build an index from a raw URI (local path or `s3://…`), dispatching through
/// the reader registry. Remote inputs stage exactly the reader's declared
/// manifest into a tempdir (kept alive across the read); local inputs are read
/// in place.
pub fn load_raw(
    uri: &str,
    backend: &dyn StagingBackend,
    cfg: &CentroidingConfig,
) -> Result<RawRead, LoadRawError> {
    let registry = ReaderRegistry::with_builtins();

    if is_remote_uri(uri) {
        let parsed = uri.parse::<Uri>().map_err(|source| ReadError::UriParse {
            uri: uri.to_string(),
            source,
        })?;
        // No reader returns `Maybe` yet, so the magic-byte peek is never needed;
        // remote byte-sniffing (a range-GET here) is future work.
        let reader = registry.pick(&parsed, || None)?;
        let manifest = reader.manifest(&parsed);
        // `staged` owns the tempdir and must outlive the read below.
        let staged = stage_manifest(backend, &manifest)?;
        let index = reader.read(staged.source(), cfg)?;
        return Ok(RawRead {
            index,
            reader_name: reader.name(),
            caches_to_idx: reader.caches_to_idx(),
        });
    }

    // Local: canonicalize first so RELATIVE paths resolve — `local_uri`
    // (sniff) and `local_in_place` (read) both require an absolute path. This
    // is the single place all entry points funnel through, so relative inputs
    // work uniformly (`read_index`, `load_index_auto`, the pyo3 binding, …).
    let abs = std::fs::canonicalize(uri)
        .map_err(|e| ReadError::Build(format!("cannot resolve local raw path {uri:?}: {e}")))?;
    let parsed = local_uri(&abs)?;
    let reader = registry.pick(&parsed, || None)?;
    let src = ResolvedSource::local_in_place(&abs)?;
    let index = reader.read(&src, cfg)?;
    Ok(RawRead {
        index,
        reader_name: reader.name(),
        caches_to_idx: reader.caches_to_idx(),
    })
}

/// A staged remote bundle: a tempdir holding the materialized manifest.
#[derive(Debug)]
pub struct StagedBundle {
    _tempdir: tempfile::TempDir,
    src: ResolvedSource,
}

impl StagedBundle {
    pub fn source(&self) -> &ResolvedSource {
        &self.src
    }
}

/// Fetch a reader's declared manifest into one tempdir, preserving each file's
/// layout relative to the entry's parent, so the reader opens `entry` and finds
/// its siblings/children alongside. `required` errors if absent; `optional` is
/// best-effort.
pub fn stage_manifest(
    backend: &dyn StagingBackend,
    manifest: &Manifest,
) -> Result<StagedBundle, StageError> {
    let step = timscentroid::TimedStep::begin("Staging manifest");

    // Layout is derived from URI PATHS (transport-agnostic), independent of how
    // `split_uri` factors local vs remote into (location, key). The entry's
    // parent is the layout root: a `.d` reconstructs as
    // `<tempdir>/sample.d/analysis.tdf`, a single file as `<tempdir>/foo.mzML`.
    let entry_path = manifest.entry.path();
    let (parent_path, entry_name) = match entry_path.rsplit_once('/') {
        Some((p, name)) => (p, name),
        None => ("", entry_path),
    };
    if entry_name.is_empty() {
        return Err(StageError::ShapeMismatch(format!(
            "manifest entry has no file name: {entry_path:?}"
        )));
    }

    let tempdir = backend.new_run_tempdir()?;
    std::fs::File::create(tempdir.path().join(".lock")).map_err(StageError::Io)?;

    for member in &manifest.required {
        fetch_member(member, parent_path, tempdir.path(), true)?;
    }
    for member in &manifest.optional {
        fetch_member(member, parent_path, tempdir.path(), false)?;
    }

    step.finish();
    let src = ResolvedSource::new(tempdir.path().to_path_buf(), OsString::from(entry_name));
    Ok(StagedBundle {
        _tempdir: tempdir,
        src,
    })
}

fn fetch_member(
    member: &Uri,
    parent_path: &str,
    tempdir: &Path,
    required: bool,
) -> Result<(), StageError> {
    // Destination = member path relative to the entry's parent. The member must
    // live under that parent and must not escape the tempdir (cheap
    // defense-in-depth; manifests are reader-built/trusted).
    let rel = member
        .path()
        .strip_prefix(parent_path)
        .map(|s| s.trim_start_matches('/'))
        .filter(|s| !s.is_empty())
        .ok_or_else(|| {
            StageError::ShapeMismatch(format!(
                "manifest member {:?} is not under the entry parent {parent_path:?}",
                member.path()
            ))
        })?;
    if rel.split('/').any(|c| c == "..") {
        return Err(StageError::ShapeMismatch(format!(
            "unsafe manifest member path: {rel:?}"
        )));
    }
    let dest = tempdir.join(rel);

    // Fetch from the member's own store + key (`split_uri` factors transport).
    let member_uri = member.to_string();
    let (loc, key) = split_uri(&member_uri)?;
    let provider = StorageProvider::open(loc).map_err(transport_err(&member_uri))?;

    if !required && !matches!(provider.exists(&key), Ok(true)) {
        return Ok(()); // best-effort optional: skip if absent
    }
    if let Some(parent) = dest.parent() {
        std::fs::create_dir_all(parent).map_err(StageError::Io)?;
    }
    let bar = make_bar(0, rel);
    provider
        .get_to_file(&key, &dest, &bar)
        .map_err(transport_err(&member_uri))?;
    bar.finish_and_clear();
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::backend::{
        PerRunTempdir,
        StagingConfig,
    };
    use timscentroid::reader::Manifest;

    fn u(p: &std::path::Path) -> Uri {
        // Raw path URI (tempdir paths are space-free, so no percent-encoding).
        p.to_str().unwrap().parse::<Uri>().unwrap()
    }

    fn dotd_manifest(dotd: &std::path::Path) -> Manifest {
        Manifest {
            entry: u(dotd),
            required: vec![
                u(&dotd.join("analysis.tdf")),
                u(&dotd.join("analysis.tdf_bin")),
            ],
            optional: vec![],
        }
    }

    #[test]
    fn stage_manifest_materializes_required_by_layout() {
        let src = tempfile::tempdir().unwrap();
        let dotd = src.path().join("sample.d");
        std::fs::create_dir(&dotd).unwrap();
        std::fs::write(dotd.join("analysis.tdf"), b"tdf").unwrap();
        std::fs::write(dotd.join("analysis.tdf_bin"), b"bin").unwrap();

        let backend = PerRunTempdir::new(StagingConfig::default()).unwrap();
        let staged = stage_manifest(&backend, &dotd_manifest(&dotd)).unwrap();
        let entry = staged.source().entry_path();

        // Reconstructed as `<tempdir>/sample.d/{analysis.tdf,analysis.tdf_bin}`.
        assert_eq!(entry.file_name().unwrap(), "sample.d");
        assert_eq!(std::fs::read(entry.join("analysis.tdf")).unwrap(), b"tdf");
        assert_eq!(
            std::fs::read(entry.join("analysis.tdf_bin")).unwrap(),
            b"bin"
        );
    }

    #[test]
    fn stage_manifest_errors_on_missing_required() {
        let src = tempfile::tempdir().unwrap();
        let dotd = src.path().join("sample.d");
        std::fs::create_dir(&dotd).unwrap();
        std::fs::write(dotd.join("analysis.tdf"), b"tdf").unwrap();
        // analysis.tdf_bin intentionally absent — must error, not silently skip.

        let backend = PerRunTempdir::new(StagingConfig::default()).unwrap();
        let err = stage_manifest(&backend, &dotd_manifest(&dotd)).unwrap_err();
        assert!(
            matches!(err, StageError::Transport { .. } | StageError::Io(_)),
            "expected a fetch error for the missing required member, got {err:?}"
        );
    }
}
