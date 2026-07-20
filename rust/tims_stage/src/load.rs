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
    let parsed = parse_uri(uri)?;
    let registry = ReaderRegistry::with_builtins();
    // No reader returns `Maybe` yet, so the magic-byte peek is never needed;
    // remote byte-sniffing (a range-GET here) is future work.
    let reader = registry.pick(&parsed, || None)?;

    // `_staged` owns the tempdir (if any) and must outlive the read.
    let (src, _staged) = if is_remote_uri(uri) {
        let manifest = reader.manifest(&parsed);
        let staged = stage_manifest(backend, &manifest)?;
        (staged.source().clone(), Some(staged))
    } else {
        (ResolvedSource::local_in_place(Path::new(uri))?, None)
    };

    let index = reader.read(&src, cfg)?;
    Ok(RawRead {
        index,
        reader_name: reader.name(),
        caches_to_idx: reader.caches_to_idx(),
    })
}

fn parse_uri(uri: &str) -> Result<Uri, ReadError> {
    if is_remote_uri(uri) {
        uri.parse::<Uri>().map_err(|source| ReadError::UriParse {
            uri: uri.to_string(),
            source,
        })
    } else {
        local_uri(Path::new(uri))
    }
}

/// A staged remote bundle: a tempdir holding the materialized manifest.
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
    let entry_uri = manifest.entry.to_string();
    let (loc, entry_key) = split_uri(&entry_uri)?;
    let provider = StorageProvider::open(loc).map_err(transport_err(&entry_uri))?;

    // Path of everything is expressed relative to the entry's parent, so a `.d`
    // reconstructs as `<tempdir>/sample.d/analysis.tdf` and a single file as
    // `<tempdir>/foo.mzML`.
    let parent_prefix = match entry_key.rsplit_once('/') {
        Some((p, _)) => format!("{p}/"),
        None => String::new(),
    };
    let entry_name = entry_key[parent_prefix.len()..].to_string();

    let tempdir = backend.new_run_tempdir()?;
    std::fs::File::create(tempdir.path().join(".lock")).map_err(StageError::Io)?;

    for member in &manifest.required {
        fetch_member(&provider, member, &parent_prefix, tempdir.path(), true)?;
    }
    for member in &manifest.optional {
        fetch_member(&provider, member, &parent_prefix, tempdir.path(), false)?;
    }

    step.finish();
    let src = ResolvedSource::new(tempdir.path().to_path_buf(), OsString::from(entry_name));
    Ok(StagedBundle {
        _tempdir: tempdir,
        src,
    })
}

fn fetch_member(
    provider: &StorageProvider,
    member: &Uri,
    parent_prefix: &str,
    tempdir: &Path,
    required: bool,
) -> Result<(), StageError> {
    let member_uri = member.to_string();
    let (_loc, key) = split_uri(&member_uri)?;
    let rel = key.strip_prefix(parent_prefix).unwrap_or(&key);
    let dest = tempdir.join(rel);

    if !required {
        // Best-effort: skip a missing optional member without erroring.
        match provider.exists(&key) {
            Ok(true) => {}
            _ => return Ok(()),
        }
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
