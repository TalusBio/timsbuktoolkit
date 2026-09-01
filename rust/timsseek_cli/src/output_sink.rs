//! Where a run's artifacts land.
//!
//! A destination is either a local directory or a remote bucket prefix, and the
//! rest of the run should not have to care which: it writes into a working
//! directory and asks the sink to finalize. Remote destinations write into a
//! tempdir and upload; local ones write in place.

use tims_stage::is_remote_uri;

use crate::errors;

/// Local-vs-remote output routing: a remote `dest_uri` (s3://, gs://, az://)
/// writes into a tempdir and uploads per-sample; a local one writes directly.
pub(crate) struct OutputSink {
    dest_uri: String,
    working_dir: std::path::PathBuf,
    remote: bool,
    _tempdir: Option<tempfile::TempDir>,
}

impl OutputSink {
    pub(crate) fn new(dest_uri: &str) -> Result<Self, errors::CliError> {
        if is_remote_uri(dest_uri) {
            let td = tempfile::Builder::new()
                .prefix("timsseek-output-")
                .tempdir()
                .map_err(|e| errors::CliError::Io {
                    source: format!("output tempdir: {e}"),
                    path: None,
                })?;
            let working_dir = td.path().to_path_buf();
            Ok(Self {
                dest_uri: dest_uri.to_string(),
                working_dir,
                remote: true,
                _tempdir: Some(td),
            })
        } else {
            std::fs::create_dir_all(dest_uri).map_err(|e| errors::CliError::Io {
                source: format!("create output dir: {e}"),
                path: Some(dest_uri.to_string()),
            })?;
            Ok(Self {
                dest_uri: dest_uri.to_string(),
                working_dir: std::path::PathBuf::from(dest_uri),
                remote: false,
                _tempdir: None,
            })
        }
    }

    fn sample_dir(&self, sample: &str) -> std::path::PathBuf {
        self.working_dir.join(sample)
    }

    pub(crate) fn root(&self) -> &std::path::Path {
        &self.working_dir
    }

    /// Where a sample's files *end up*, not the working tempdir. Use for
    /// user-facing output so users don't see a tempdir that will be wiped.
    pub(crate) fn dest_uri_for_sample(&self, sample: &str) -> String {
        if self.remote {
            format!("{}/{}", self.dest_uri.trim_end_matches('/'), sample)
        } else {
            self.sample_dir(sample).to_string_lossy().into_owned()
        }
    }

    /// Upload and remove a per-sample subdir after the sample has finished
    /// writing; no-op for local destinations.
    pub(crate) fn finalize_sample(&self, sample: &str) -> Result<(), errors::CliError> {
        if !self.remote {
            return Ok(());
        }
        let local = self.sample_dir(sample);
        let sample_dest = self.dest_uri_for_sample(sample);
        for entry in std::fs::read_dir(&local).map_err(|e| errors::CliError::Io {
            source: format!("read sample dir: {e}"),
            path: Some(local.to_string_lossy().to_string()),
        })? {
            let entry = entry.map_err(|e| errors::CliError::Io {
                source: format!("read dir entry: {e}"),
                path: None,
            })?;
            let bn = entry.file_name().to_string_lossy().to_string();
            let dest = format!("{sample_dest}/{bn}");
            tims_stage::upload_file(&entry.path(), &dest).map_err(|e| errors::CliError::Io {
                source: format!("upload {dest}: {e}"),
                path: None,
            })?;
        }
        std::fs::remove_dir_all(&local).map_err(|e| errors::CliError::Io {
            source: format!("cleanup sample dir: {e}"),
            path: Some(local.to_string_lossy().to_string()),
        })?;
        Ok(())
    }

    /// Upload named top-level files (run_report.json, config_used.json)
    /// that exist in the working dir; no-op for local destinations.
    pub(crate) fn finalize_run(&self, files: &[&str]) -> Result<(), errors::CliError> {
        if !self.remote {
            return Ok(());
        }
        for bn in files {
            let local = self.working_dir.join(bn);
            if !local.exists() {
                continue;
            }
            let dest = format!("{}/{}", self.dest_uri.trim_end_matches('/'), bn);
            tims_stage::upload_file(&local, &dest).map_err(|e| errors::CliError::Io {
                source: format!("upload {dest}: {e}"),
                path: None,
            })?;
        }
        Ok(())
    }
}

/// Join an artifact path onto a base output URI. Remote URIs get a plain
/// `base/rel` concat; local paths go through `PathBuf::join` so OS-specific
/// separators are respected.
pub(crate) fn join_output_uri(base: &str, rel: &str) -> String {
    if is_remote_uri(base) {
        format!("{}/{}", base.trim_end_matches('/'), rel)
    } else {
        std::path::Path::new(base)
            .join(rel)
            .to_string_lossy()
            .to_string()
    }
}

/// Cheap existence probe: local `Path::exists` or remote HEAD.
pub(crate) fn probe_uri_exists(uri: &str) -> Result<bool, errors::CliError> {
    tims_stage::uri_exists(uri).map_err(|e| errors::CliError::Io {
        source: format!("existence probe {uri}: {e}"),
        path: Some(uri.to_string()),
    })
}

/// `sample.d.tar`, `sample.d/`, `sample.d.idx/` all collapse to `sample`.
pub(crate) fn sample_name_from_uri(uri: &str) -> Option<String> {
    let trimmed = uri.trim_end_matches('/');
    let mut stem = trimmed.rsplit('/').next()?;
    // Loop so chained suffixes collapse fully. Order matters: `.idx`/`.tar`
    // come off before `.d` so they can't leave a bare `.d` behind.
    loop {
        let before = stem;
        for ext in [".idx", ".tar", ".d"] {
            if let Some(s) = stem.strip_suffix(ext) {
                stem = s;
            }
        }
        if stem == before {
            break;
        }
    }
    if stem.is_empty() {
        None
    } else {
        Some(stem.to_string())
    }
}

#[cfg(test)]
mod sample_name_tests {
    use super::sample_name_from_uri;
    #[test]
    fn local_dotd_plain() {
        assert_eq!(sample_name_from_uri("/data/run.d").as_deref(), Some("run"));
    }
    #[test]
    fn local_dotd_trailing_slash() {
        assert_eq!(sample_name_from_uri("/data/run.d/").as_deref(), Some("run"));
    }
    #[test]
    fn s3_tar_collapses_both_suffixes() {
        assert_eq!(
            sample_name_from_uri("s3://bkt/run.d.tar").as_deref(),
            Some("run")
        );
    }
    #[test]
    fn s3_idx_directory() {
        assert_eq!(
            sample_name_from_uri("s3://bkt/run.d.idx/").as_deref(),
            Some("run")
        );
    }
}
