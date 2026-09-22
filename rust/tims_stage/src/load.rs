//! The single raw-index composition: URI → `IndexedTimstofPeaks`.
//!
//! `PreparedSource` retains reader selection from preflight; `load_raw` selects
//! on demand for callers without preflight. Both ask that reader for its
//! [`Manifest`], materialize it -- in place for local inputs, or by fetching
//! exactly the declared files for remote ones -- and call `read`. Transport
//! never guesses vendor shape: it fetches what the reader declared, by name.
//!
//! Lives in `tims_stage` (one crate above `timscentroid`, where the registry
//! lives) so it can compose `pick` + `manifest` + `read` with the fetch engine
//! without a dependency cycle.

use std::ffi::OsString;
use std::path::Path;
use std::sync::Arc;

use http::Uri;
use timscentroid::reader::{
    BrukerTdfReader,
    Manifest,
    RawReader,
    ReadError,
    ReaderRegistry,
    ResolvedSource,
    local_uri,
};
use timscentroid::{
    IndexedTimstofPeaks,
    IndexingCentroidingConfig,
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

/// Reader selection retained from preflight through naming and index loading.
/// Standalone indexes may have no raw reader; tar transport contains Bruker .d.
#[derive(Clone)]
pub struct PreparedSource {
    uri: String,
    sample_name: String,
    reader: Option<Arc<dyn RawReader>>,
}

impl std::fmt::Debug for PreparedSource {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PreparedSource")
            .field("uri", &self.uri)
            .field("sample_name", &self.sample_name)
            .field("reader", &self.reader.as_ref().map(|r| r.name()))
            .finish()
    }
}

impl PreparedSource {
    pub fn new(uri: &str) -> Result<Self, ReadError> {
        let uri = crate::canonical_uri(uri);
        let location = uri.trim_end_matches('/');
        let mut entry = location;
        let mut wrapped = false;
        let tar = location.to_ascii_lowercase().ends_with(".tar");
        let idx = location.to_ascii_lowercase().ends_with(".idx");
        loop {
            let lower = entry.to_ascii_lowercase();
            if lower.ends_with(".idx") || lower.ends_with(".tar") {
                entry = &entry[..entry.len() - 4];
                wrapped = true;
            } else {
                break;
            }
        }
        let parsed = raw_uri(entry)?;
        let reader = if tar {
            // The existing tar loader is specifically a .d transport, not a
            // second vendor-dispatch mechanism inferred from the archive name.
            Some(Arc::new(BrukerTdfReader) as Arc<dyn RawReader>)
        } else {
            match ReaderRegistry::with_builtins().pick(&parsed, || None) {
                Ok(reader) => Some(reader),
                Err(ReadError::UnknownFormat(_)) if idx => None,
                Err(error) => return Err(error),
            }
        };
        Self::with_reader(&uri, entry, wrapped, reader)
    }

    fn with_reader(
        uri: &str,
        entry: &str,
        wrapped: bool,
        reader: Option<Arc<dyn RawReader>>,
    ) -> Result<Self, ReadError> {
        let name = entry.rsplit('/').next().unwrap_or(entry);
        let sample_name = reader
            .as_ref()
            .and_then(|r| r.sample_name(name))
            .or_else(|| wrapped.then_some(name))
            .ok_or_else(|| ReadError::Build(format!("selected reader cannot name {uri:?}")))?
            .to_owned();
        Ok(Self {
            uri: uri.to_owned(),
            sample_name,
            reader,
        })
    }

    pub fn uri(&self) -> &str {
        &self.uri
    }

    pub fn sample_name(&self) -> &str {
        &self.sample_name
    }

    /// Read the raw artifact (or the .d extracted by tar staging) using the
    /// same reader which supplied its name. No registry dispatch here.
    pub fn read_raw(
        &self,
        uri: &str,
        backend: &dyn StagingBackend,
        cfg: &IndexingCentroidingConfig,
    ) -> Result<RawRead, LoadRawError> {
        let reader = self
            .reader
            .as_ref()
            .ok_or_else(|| ReadError::Build("cached index has no raw reader".into()))?;
        load_with_reader(uri, reader.as_ref(), backend, cfg)
    }
}

fn raw_uri(uri: &str) -> Result<Uri, ReadError> {
    if is_remote_uri(uri) {
        uri.parse().map_err(|source| ReadError::UriParse {
            uri: uri.to_owned(),
            source,
        })
    } else {
        local_uri(&std::path::absolute(uri)?)
    }
}

/// Build an index from a raw URI (local path or `s3://…`), dispatching through
/// the reader registry. Remote inputs stage exactly the reader's declared
/// manifest into a tempdir (kept alive across the read); local inputs are read
/// in place.
pub fn load_raw(
    uri: &str,
    backend: &dyn StagingBackend,
    cfg: &IndexingCentroidingConfig,
) -> Result<RawRead, LoadRawError> {
    let parsed = raw_uri(&crate::canonical_uri(uri))?;
    let reader = ReaderRegistry::with_builtins().pick(&parsed, || None)?;
    load_with_reader(uri, reader.as_ref(), backend, cfg)
}

fn load_with_reader(
    uri: &str,
    reader: &dyn RawReader,
    backend: &dyn StagingBackend,
    cfg: &IndexingCentroidingConfig,
) -> Result<RawRead, LoadRawError> {
    if is_remote_uri(uri) {
        let parsed = uri.parse::<Uri>().map_err(|source| ReadError::UriParse {
            uri: uri.to_string(),
            source,
        })?;
        // No reader returns `Maybe` yet, so the magic-byte peek is never needed;
        // remote byte-sniffing (a range-GET here) is future work.
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

    // Resolve the actual local artifact without reselecting its reader.
    let abs = std::fs::canonicalize(uri)
        .map_err(|e| ReadError::Build(format!("cannot resolve local raw path {uri:?}: {e}")))?;
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
        // Require a component boundary so `/data/run` doesn't match `/data/run123`.
        .filter(|r| parent_path.is_empty() || r.starts_with('/'))
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

    #[test]
    fn prepared_names_unwrap_transport() {
        for name in [
            "My-Run.D",
            "My-Run.d.IDX",
            "My-Run.D.TAR",
            "My-Run.idx",
            "My-Run.tar",
        ] {
            assert_eq!(PreparedSource::new(name).unwrap().sample_name(), "My-Run");
        }
        assert_eq!(
            PreparedSource::new("run.raw.idx").unwrap().sample_name(),
            "run.raw"
        );
        assert!(PreparedSource::new("run.raw").is_err());
        assert!(PreparedSource::new("run.mzML.gz").is_err());
    }

    #[test]
    fn prepared_source_names_and_reads_with_retained_reader() {
        struct SelectedReader;
        impl RawReader for SelectedReader {
            fn name(&self) -> &'static str {
                "selected"
            }

            fn sniff(&self, _: &Uri) -> timscentroid::reader::Sniff {
                panic!("must not redispatch after selection")
            }

            fn sample_name<'a>(&self, name: &'a str) -> Option<&'a str> {
                name.strip_suffix(".custom")
            }

            fn manifest(&self, _: &Uri) -> Manifest {
                unreachable!("local input")
            }

            fn read(
                &self,
                src: &ResolvedSource,
                _: &IndexingCentroidingConfig,
            ) -> Result<IndexedTimstofPeaks, ReadError> {
                assert_eq!(src.entry_path().file_name().unwrap(), "run.custom");
                Err(ReadError::Build("retained reader called".into()))
            }
        }
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("run.custom");
        std::fs::write(&path, []).unwrap();
        let uri = path.to_str().unwrap();
        let source =
            PreparedSource::with_reader(uri, uri, false, Some(Arc::new(SelectedReader))).unwrap();
        assert_eq!(source.sample_name(), "run");
        let backend = PerRunTempdir::new(StagingConfig::default()).unwrap();
        assert!(
            matches!(source.read_raw(uri, &backend, &IndexingCentroidingConfig::default()),
            Err(LoadRawError::Read(ReadError::Build(message))) if message == "retained reader called")
        );
    }

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
        // analysis.tdf_bin intentionally absent -- must error, not silently skip.

        let backend = PerRunTempdir::new(StagingConfig::default()).unwrap();
        let err = stage_manifest(&backend, &dotd_manifest(&dotd)).unwrap_err();
        assert!(
            matches!(err, StageError::Transport { .. } | StageError::Io(_)),
            "expected a fetch error for the missing required member, got {err:?}"
        );
    }
}
