//! StagingBackend trait, PerRunTempdir implementation, StagedDotD handle,
//! and the `run_on_staged` ergonomic wrapper.

use crate::error::StageError;
use crate::resolve::SourceSpec;
use std::path::{
    Path,
    PathBuf,
};

pub trait StagingBackend {
    fn stage(&self, spec: &SourceSpec) -> Result<StagedDotD, StageError>;
}

#[derive(Debug, Clone)]
pub struct StagingConfig {
    /// None = system default temp dir.
    pub tempdir_root: Option<PathBuf>,
    /// Hard cap on prefix listing size.
    pub max_prefix_keys: usize,
    /// Age threshold for the startup sweep. `0` disables the sweep.
    pub stale_sweep_age_hours: u64,
}

impl Default for StagingConfig {
    fn default() -> Self {
        Self {
            tempdir_root: None,
            max_prefix_keys: 256,
            stale_sweep_age_hours: 24,
        }
    }
}

/// RAII handle over a per-run tempdir + its inner `sample.d/` path.
///
/// Callers MUST keep `StagedDotD` alive for as long as they use the path —
/// the tempdir is deleted on drop, and `AsRef<Path>` does not extend the
/// borrow beyond the handle's lifetime. Do NOT store `&Path` derived from
/// `as_ref()` past the end of the scope that holds the handle.
pub struct StagedDotD {
    /// The tempdir owning `dotd`. `None` when the source is a local `.d`
    /// that we borrow in place (no tempdir to drop).
    _tempdir: Option<tempfile::TempDir>,
    dotd: PathBuf,
}

impl StagedDotD {
    pub(crate) fn owned(tempdir: tempfile::TempDir, dotd: PathBuf) -> Self {
        Self {
            _tempdir: Some(tempdir),
            dotd,
        }
    }

    #[allow(dead_code)]
    pub(crate) fn borrowed(dotd: PathBuf) -> Self {
        Self {
            _tempdir: None,
            dotd,
        }
    }
}

impl AsRef<Path> for StagedDotD {
    fn as_ref(&self) -> &Path {
        &self.dotd
    }
}

/// MVP-only staging backend: per-run tempdir, deleted when the
/// `StagedDotD` drops. The Phase 9 startup sweep (added later) will reclaim
/// any tempdirs leaked by abnormal termination.
pub struct PerRunTempdir {
    cfg: StagingConfig,
}

impl PerRunTempdir {
    pub fn new(cfg: StagingConfig) -> Result<Self, StageError> {
        if cfg.stale_sweep_age_hours > 0 {
            let root = cfg.tempdir_root.clone().unwrap_or_else(std::env::temp_dir);
            sweep_stale(&root, cfg.stale_sweep_age_hours)?;
        }
        Ok(Self { cfg })
    }

    pub fn config(&self) -> &StagingConfig {
        &self.cfg
    }

    /// Create a namespaced tempdir under the configured root. The
    /// `timsseek-staging-` prefix is what the startup sweep uses to
    /// identify our tempdirs.
    pub(crate) fn new_tempdir(&self) -> Result<tempfile::TempDir, StageError> {
        let builder = tempfile::Builder::new()
            .prefix("timsseek-staging-")
            .to_owned();
        let td = match &self.cfg.tempdir_root {
            Some(root) => builder.tempdir_in(root).map_err(StageError::Io)?,
            None => builder.tempdir().map_err(StageError::Io)?,
        };
        Ok(td)
    }
}

impl StagingBackend for PerRunTempdir {
    fn stage(&self, spec: &SourceSpec) -> Result<StagedDotD, StageError> {
        match spec {
            SourceSpec::S3Tar { .. } => crate::tar::stage_s3_tar(self, spec),
            SourceSpec::S3Prefix { .. } => crate::prefix::stage_s3_prefix(self, spec),
            SourceSpec::LocalTar { .. } => crate::tar::stage_local_tar(self, spec),
        }
    }
}

/// Ergonomic wrapper — stage, run the supplied closure against the local
/// `.d` path, drop the tempdir on return. Both success and error paths
/// drop the handle. Free function (not a trait method) so `StagingBackend`
/// stays dyn-safe.
///
/// Abnormal mid-operation termination (SIGKILL, crash) leaves the tempdir;
/// the startup sweep reclaims it on the next run.
pub fn run_on_staged<B, W, F>(backend: &B, spec: &SourceSpec, op: F) -> Result<W, StageError>
where
    B: StagingBackend + ?Sized,
    F: FnOnce(&Path) -> Result<W, StageError>,
{
    let staged = backend.stage(spec)?;
    op(staged.as_ref())
    // staged drops here regardless of Ok or Err path
}

fn sweep_stale(root: &Path, age_hours: u64) -> Result<(), StageError> {
    let threshold = std::time::SystemTime::now()
        .checked_sub(std::time::Duration::from_secs(age_hours * 3600))
        .unwrap_or(std::time::UNIX_EPOCH);
    let rd = match std::fs::read_dir(root) {
        Ok(x) => x,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(()),
        Err(e) => return Err(StageError::Io(e)),
    };
    for entry in rd.flatten() {
        let name = entry.file_name().to_string_lossy().to_string();
        if !name.starts_with("timsseek-staging-") {
            continue;
        }
        let meta = match entry.metadata() {
            Ok(m) => m,
            Err(_) => continue,
        };
        if !meta.is_dir() {
            continue;
        }
        let path = entry.path();
        let has_lock = path.join(".lock").exists();
        let mtime = meta.modified().ok();
        if has_lock {
            // Live session. Only force-remove if it's very stale (2x threshold)
            // AND still claims to be locked — best-effort for crashed runs.
            if let Some(mt) = mtime {
                let hard_threshold = std::time::SystemTime::now()
                    .checked_sub(std::time::Duration::from_secs(age_hours * 2 * 3600))
                    .unwrap_or(std::time::UNIX_EPOCH);
                if mt < hard_threshold {
                    tracing::info!(?path, "removing very-stale locked tempdir");
                    let _ = std::fs::remove_dir_all(&path);
                }
            }
            continue;
        }
        if let Some(mt) = mtime {
            if mt < threshold {
                tracing::info!(?path, "sweeping stale tempdir");
                let _ = std::fs::remove_dir_all(&path);
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_config_has_256_cap_and_24h_sweep() {
        let c = StagingConfig::default();
        assert_eq!(c.max_prefix_keys, 256);
        assert_eq!(c.stale_sweep_age_hours, 24);
    }

    #[test]
    fn tempdir_name_has_namespace_prefix() {
        let backend = PerRunTempdir::new(StagingConfig::default()).unwrap();
        let td = backend.new_tempdir().unwrap();
        let name = td.path().file_name().unwrap().to_string_lossy().to_string();
        assert!(name.starts_with("timsseek-staging-"), "got {}", name);
    }
}
