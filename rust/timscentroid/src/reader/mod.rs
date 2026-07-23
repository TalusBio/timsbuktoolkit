//! Vendor-neutral raw-format dispatch.
//!
//! A [`RawReader`] is a self-describing backend: it declares which URIs it
//! claims ([`RawReader::sniff`]), which artifacts belong together
//! ([`RawReader::manifest`]), and how to build the in-memory index
//! ([`RawReader::read`]). The [`ReaderRegistry`] picks a backend for a URI by
//! sniffing — vendor suffix/scheme knowledge lives ONLY in each reader, never
//! in the staging/resolve layer. Adding a format is one `impl RawReader` + one
//! registry push.

#[cfg(feature = "mzdata")]
pub mod mzdata;
pub mod uri_util;

use std::ffi::OsString;
use std::path::{
    Path,
    PathBuf,
};

use http::Uri;
use url::Url;

use crate::centroiding::IndexingCentroidingConfig;
use crate::indexing::IndexedTimstofPeaks;

/// Outcome of a cheap URI-only claim check.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sniff {
    /// This reader claims the URI outright (suffix/scheme match).
    Yes,
    /// This reader does not claim the URI.
    No,
    /// The suffix is ambiguous; defer to a magic-byte peek (`sniff_bytes`).
    Maybe,
}

/// The exact artifacts a reader reads, and which one it opens.
///
/// Staging materializes `required` (error if any absent) + `optional`
/// (best-effort) into one local dir; `read` opens `entry` within it. Anything
/// else next to the source is ignored.
#[derive(Debug, Clone)]
pub struct Manifest {
    pub entry: Uri,
    pub required: Vec<Uri>,
    pub optional: Vec<Uri>,
}

/// One local dir holding the artifact + its co-located members.
///
/// `entry` is always a NAME within `dir` (never the dir itself), so
/// `entry_path()` is unambiguous for single-file AND directory-entry shapes.
/// For local inputs `dir` is the artifact's real parent; for staged remote
/// inputs it is a tempdir holding the materialized manifest.
#[derive(Debug, Clone)]
pub struct ResolvedSource {
    dir: PathBuf,
    entry: OsString,
}

impl ResolvedSource {
    /// Construct from a local directory + entry name.
    pub fn new(dir: PathBuf, entry: OsString) -> Self {
        Self { dir, entry }
    }

    /// Borrow a local artifact in place — no staging. `dir` is the artifact's
    /// parent, `entry` its own name (a `.d` dir or an `.mzML` file).
    pub fn local_in_place(path: &Path) -> Result<Self, ReadError> {
        let dir = path
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .to_path_buf();
        let entry = path
            .file_name()
            .ok_or_else(|| ReadError::Build(format!("path has no file name: {path:?}")))?
            .to_os_string();
        Ok(Self::new(dir, entry))
    }

    /// The artifact the reader opens (e.g. `foo.mzML`, or the `sample.d` dir).
    pub fn entry_path(&self) -> PathBuf {
        self.dir.join(&self.entry)
    }

    /// Resolve a co-located member (a sidecar like `foo.wiff.scan`) to its local
    /// path, if it exists next to the entry. Works for both in-place local
    /// inputs and staged tempdirs.
    pub fn path_of(&self, name: &str) -> Option<PathBuf> {
        let p = self.dir.join(name);
        p.exists().then_some(p)
    }
}

/// Errors from format dispatch and index construction.
#[derive(Debug, thiserror::Error)]
pub enum ReadError {
    #[error("no registered reader claims URI: {0}")]
    UnknownFormat(String),
    #[error("ambiguous format: both {first} and {second} claim the URI")]
    AmbiguousFormat {
        first: &'static str,
        second: &'static str,
    },
    #[error("could not parse URI {uri:?}: {source}")]
    UriParse {
        uri: String,
        source: http::uri::InvalidUri,
    },
    #[error("required manifest member absent: {0}")]
    MissingRequired(String),
    #[error("index build failed: {0}")]
    Build(String),
    #[error(transparent)]
    Io(#[from] std::io::Error),
}

/// A self-describing raw-format backend. See module docs.
pub trait RawReader: Send + Sync {
    /// Name for logs/errors.
    fn name(&self) -> &'static str;

    /// Cheap claim check from the URI alone (scheme/suffix). `Maybe` defers to
    /// a magic-byte peek.
    fn sniff(&self, uri: &Uri) -> Sniff;

    /// Confirm a `Maybe` by inspecting leading bytes. Only called when
    /// `sniff` == `Maybe`.
    fn sniff_bytes(&self, _head: &[u8]) -> bool {
        false
    }

    /// The exact artifact set (entry + required + optional). Vendor "which
    /// files belong together" knowledge lives here, nowhere else — the staging
    /// layer fetches exactly these, by name, so transport never guesses shape.
    fn manifest(&self, uri: &Uri) -> Manifest;

    /// Whether a built index should be persisted as an `.idx` sidecar for reuse.
    /// TDF caches; mzdata builds in memory each run (no sidecar).
    fn caches_to_idx(&self) -> bool {
        true
    }

    /// Extensions this reader claims for single-FILE inputs, for building
    /// file-dialog / help filters (no leading dot, e.g. `["mzML", "mzml"]`).
    /// Directory formats (`.d`) return `&[]` — they are picked as folders.
    fn file_extensions(&self) -> &'static [&'static str] {
        &[]
    }

    /// Build the index. `src` is one local dir holding the staged manifest;
    /// open `src.entry_path()`.
    fn read(
        &self,
        src: &ResolvedSource,
        cfg: &IndexingCentroidingConfig,
    ) -> Result<IndexedTimstofPeaks, ReadError>;
}

/// Ordered set of backends. First `Yes` wins.
pub struct ReaderRegistry(Vec<Box<dyn RawReader>>);

impl ReaderRegistry {
    /// The built-in readers. `BrukerTdfReader` is always present; format
    /// backends behind features are pushed under their `#[cfg]`.
    pub fn with_builtins() -> Self {
        #[allow(unused_mut)]
        let mut v: Vec<Box<dyn RawReader>> = vec![Box::new(BrukerTdfReader)];
        #[cfg(feature = "mzdata")]
        v.push(Box::new(crate::reader::mzdata::MzdataReader));
        Self(v)
    }

    /// The union of single-file extensions claimed by all registered readers —
    /// the single source of truth for "what raw files can we open" (dialog
    /// filters, help text). Directory formats (`.d`) are excluded.
    pub fn file_extensions(&self) -> Vec<&'static str> {
        self.0
            .iter()
            .flat_map(|r| r.file_extensions().iter().copied())
            .collect()
    }

    /// Select the backend for `uri`. First `Yes` wins; `Maybe`s are
    /// disambiguated by a single magic-byte peek (`head`); no match →
    /// `UnknownFormat`; two outright `Yes` → `AmbiguousFormat`.
    pub fn pick(
        &self,
        uri: &Uri,
        head: impl FnOnce() -> Option<Vec<u8>>,
    ) -> Result<&dyn RawReader, ReadError> {
        let mut winner: Option<&dyn RawReader> = None;
        let mut maybes: Vec<&dyn RawReader> = Vec::new();
        for r in &self.0 {
            match r.sniff(uri) {
                Sniff::Yes => match winner {
                    None => winner = Some(r.as_ref()),
                    Some(w) => {
                        return Err(ReadError::AmbiguousFormat {
                            first: w.name(),
                            second: r.name(),
                        });
                    }
                },
                Sniff::Maybe => maybes.push(r.as_ref()),
                Sniff::No => {}
            }
        }
        if let Some(w) = winner {
            return Ok(w);
        }
        if !maybes.is_empty()
            && let Some(bytes) = head()
        {
            for r in maybes {
                if r.sniff_bytes(&bytes) {
                    return Ok(r);
                }
            }
        }
        Err(ReadError::UnknownFormat(uri.to_string()))
    }
}

/// Lowercased path suffix test for a URI (scheme-agnostic; keys on the path).
fn path_ends_with(uri: &Uri, suffix_lower: &str) -> bool {
    uri.path()
        .trim_end_matches('/')
        .to_ascii_lowercase()
        .ends_with(suffix_lower)
}

/// Build a `file://` URI from an absolute local path, percent-encoding as
/// needed (spaces etc.). Used to hand a local artifact to the registry for
/// sniff/manifest; the actual file access goes through [`ResolvedSource`], which
/// keeps the original `PathBuf` (no decode round-trip).
pub fn local_uri(path: &Path) -> Result<Uri, ReadError> {
    let url = Url::from_file_path(path)
        .map_err(|()| ReadError::Build(format!("not an absolute local path: {path:?}")))?;
    // Parse the percent-encoded PATH only (schemeless): `http::Uri` rejects the
    // empty authority in `file:///...`, and sniff/manifest key on the path.
    let encoded = url.path();
    encoded
        .parse::<Uri>()
        .map_err(|source| ReadError::UriParse {
            uri: encoded.to_string(),
            source,
        })
}

/// Bruker TDF `.d` reader — wraps [`IndexedTimstofPeaks::from_timstof_file`].
///
/// MVP regression-safety: the manifest declares the minimal read set
/// (`analysis.tdf` + `analysis.tdf_bin`) for documentation and future
/// selective-fetch, but `read` opens the whole `.d` dir exactly as before.
pub struct BrukerTdfReader;

impl RawReader for BrukerTdfReader {
    fn name(&self) -> &'static str {
        "bruker-tdf"
    }

    fn sniff(&self, uri: &Uri) -> Sniff {
        if path_ends_with(uri, ".d") {
            Sniff::Yes
        } else {
            Sniff::No
        }
    }

    fn manifest(&self, uri: &Uri) -> Manifest {
        Manifest {
            entry: uri.clone(),
            required: vec![
                uri_util::child(uri, "analysis.tdf"),
                uri_util::child(uri, "analysis.tdf_bin"),
            ],
            optional: vec![],
        }
    }

    fn read(
        &self,
        src: &ResolvedSource,
        cfg: &IndexingCentroidingConfig,
    ) -> Result<IndexedTimstofPeaks, ReadError> {
        let path = src.entry_path();
        let path_str = path
            .to_str()
            .ok_or_else(|| ReadError::Build(format!("path is not valid UTF-8: {path:?}")))?;
        let tt = crate::indexing::TimsTofPath::new(path_str)
            .map_err(|e| ReadError::Build(format!("failed to open {path:?}: {e:?}")))?;
        let (idx, _stats) = IndexedTimstofPeaks::from_timstof_file(&tt, *cfg);
        Ok(idx)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bruker_sniffs_dotd_yes_else_no() {
        let b = BrukerTdfReader;
        assert_eq!(b.sniff(&"/data/sample.d".parse().unwrap()), Sniff::Yes);
        assert_eq!(b.sniff(&"/data/sample.d/".parse().unwrap()), Sniff::Yes);
        assert_eq!(b.sniff(&"s3://bkt/sample.d".parse().unwrap()), Sniff::Yes);
        assert_eq!(b.sniff(&"/data/foo.mzML".parse().unwrap()), Sniff::No);
    }

    #[test]
    fn bruker_manifest_declares_minimal_read_set() {
        let b = BrukerTdfReader;
        let m = b.manifest(&"/data/sample.d".parse().unwrap());
        assert_eq!(m.entry.path(), "/data/sample.d");
        let req: Vec<_> = m.required.iter().map(|u| u.path()).collect();
        assert_eq!(
            req,
            vec![
                "/data/sample.d/analysis.tdf",
                "/data/sample.d/analysis.tdf_bin"
            ]
        );
    }

    #[cfg(feature = "mzdata")]
    #[test]
    fn registry_file_extensions_include_mzml_but_not_dotd() {
        // Real invariant: single-file formats surface in the file-dialog filter;
        // directory formats (.d) never do.
        let exts = ReaderRegistry::with_builtins().file_extensions();
        assert!(exts.contains(&"mzML"), "expected mzML in {exts:?}");
        assert!(!exts.iter().any(|e| e.eq_ignore_ascii_case("d")));
    }

    #[test]
    fn registry_picks_bruker_for_dotd() {
        let reg = ReaderRegistry::with_builtins();
        let name = reg
            .pick(&"/data/sample.d".parse().unwrap(), || None)
            .map(|r| r.name())
            .expect("bruker claims .d");
        assert_eq!(name, "bruker-tdf");
    }

    #[test]
    fn registry_unknown_format_for_unclaimed_suffix() {
        let reg = ReaderRegistry::with_builtins();
        // .txt is claimed by nobody.
        match reg.pick(&"/data/foo.txt".parse().unwrap(), || None) {
            Err(ReadError::UnknownFormat(_)) => {}
            other => panic!("expected UnknownFormat, got {:?}", other.map(|r| r.name())),
        }
    }

    #[test]
    fn local_uri_percent_encodes_spaces_and_sniffs() {
        let u = local_uri(Path::new("/data/my runs/foo.mzML")).unwrap();
        assert!(u.path().ends_with("/foo.mzML"));
        assert!(u.path().contains("%20"));
        // The registry can still sniff it once mzdata is present; Bruker says No.
        assert_eq!(BrukerTdfReader.sniff(&u), Sniff::No);
        let d = local_uri(Path::new("/data/sample.d")).unwrap();
        assert_eq!(BrukerTdfReader.sniff(&d), Sniff::Yes);
    }

    #[test]
    fn resolved_source_local_in_place() {
        let src = ResolvedSource::local_in_place(Path::new("/data/sample.d")).unwrap();
        assert_eq!(src.entry_path(), PathBuf::from("/data/sample.d"));
        let src = ResolvedSource::local_in_place(Path::new("/data/foo.mzML")).unwrap();
        assert_eq!(src.entry_path(), PathBuf::from("/data/foo.mzML"));
    }

    #[test]
    fn resolved_source_entry_and_path_of() {
        // path_of checks the real filesystem (works for staged + in-place).
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("run.wiff"), b"x").unwrap();
        std::fs::write(dir.path().join("run.wiff.scan"), b"y").unwrap();
        let src = ResolvedSource::new(dir.path().to_path_buf(), OsString::from("run.wiff"));
        assert_eq!(src.entry_path(), dir.path().join("run.wiff"));
        assert_eq!(
            src.path_of("run.wiff.scan"),
            Some(dir.path().join("run.wiff.scan"))
        );
        assert_eq!(src.path_of("missing.bin"), None);
    }
}
