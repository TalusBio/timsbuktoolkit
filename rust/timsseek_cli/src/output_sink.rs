//! Where a run's artifacts land.
//!
//! A destination is either a local directory or a remote bucket prefix, and the
//! rest of the run should not have to care which: it writes into a working
//! directory and asks the sink to finalize. Remote destinations write into a
//! tempdir and upload; local ones write in place.

use tims_stage::is_remote_uri;

use crate::errors;
/// Local-vs-remote output routing: a remote destination (s3://, gs://, az://)
/// writes into a tempdir and uploads per-sample; a local one writes directly.
///
/// The tempdir is owned by the `Remote` variant so it outlives the uploads and
/// so a local destination has no tempdir to accidentally consult.
enum Destination {
    Local(std::path::PathBuf),
    Remote {
        uri: String,
        tempdir: tempfile::TempDir,
    },
}

pub(crate) struct OutputSink {
    dest: Destination,
}

impl OutputSink {
    pub(crate) fn new(dest_uri: &str) -> Result<Self, errors::CliError> {
        let dest = if is_remote_uri(dest_uri) {
            let tempdir = tempfile::Builder::new()
                .prefix("timsseek-output-")
                .tempdir()
                .map_err(|e| errors::CliError::Io {
                    source: format!("output tempdir: {e}"),
                    path: None,
                })?;
            Destination::Remote {
                uri: dest_uri.to_string(),
                tempdir,
            }
        } else {
            std::fs::create_dir_all(dest_uri).map_err(|e| errors::CliError::Io {
                source: format!("create output dir: {e}"),
                path: Some(dest_uri.to_string()),
            })?;
            Destination::Local(std::path::PathBuf::from(dest_uri))
        };
        Ok(Self { dest })
    }

    pub(crate) fn sample_dir(&self, sample: &str) -> std::path::PathBuf {
        self.root().join(sample)
    }

    /// Where the run writes: the destination itself when it is local, the
    /// staging tempdir when it is remote.
    pub(crate) fn root(&self) -> &std::path::Path {
        match &self.dest {
            Destination::Local(dir) => dir,
            Destination::Remote { tempdir, .. } => tempdir.path(),
        }
    }

    /// Where the run's artifacts *end up*, not the working tempdir. Use for
    /// user-facing output so users don't see a tempdir that will be wiped.
    pub(crate) fn dest_root(&self) -> String {
        match &self.dest {
            Destination::Local(dir) => dir.to_string_lossy().into_owned(),
            Destination::Remote { uri, .. } => uri.trim_end_matches('/').to_string(),
        }
    }

    /// The only place a relative artifact path is joined onto the destination,
    /// so local and remote spellings cannot drift apart.
    pub(crate) fn dest_uri_for(&self, rel: &str) -> String {
        join_output_uri(&self.dest_root(), rel)
    }

    /// Remove whatever an earlier run left in a sample's directory, so fresh
    /// results are never mixed with older ones.
    pub(crate) fn clear_existing(&self, sample: &str) -> Result<(), errors::CliError> {
        let dir = self.sample_dir(sample);
        for artifact in crate::artifacts::per_sample_artifacts(true) {
            let path = dir.join(artifact);
            if !path.exists() {
                continue;
            }
            std::fs::remove_file(&path).map_err(|e| errors::CliError::Io {
                source: format!("Failed to remove existing {artifact}: {e}"),
                path: Some(path.to_string_lossy().to_string()),
            })?;
        }
        Ok(())
    }

    /// Upload and remove a per-sample subdir after the sample has finished
    /// writing; no-op for local destinations, which wrote in place.
    pub(crate) fn finalize_sample(&self, sample: &str) -> Result<(), errors::CliError> {
        match &self.dest {
            Destination::Local(_) => Ok(()),
            Destination::Remote { .. } => {
                let local = self.sample_dir(sample);
                for entry in std::fs::read_dir(&local).map_err(|e| errors::CliError::Io {
                    source: format!("read sample dir: {e}"),
                    path: Some(local.to_string_lossy().to_string()),
                })? {
                    let entry = entry.map_err(|e| errors::CliError::Io {
                        source: format!("read dir entry: {e}"),
                        path: None,
                    })?;
                    let bn = entry.file_name().to_string_lossy().to_string();
                    let dest = self.dest_uri_for(&format!("{sample}/{bn}"));
                    tims_stage::upload_file(&entry.path(), &dest).map_err(|e| {
                        errors::CliError::Io {
                            source: format!("upload {dest}: {e}"),
                            path: None,
                        }
                    })?;
                }
                std::fs::remove_dir_all(&local).map_err(|e| errors::CliError::Io {
                    source: format!("cleanup sample dir: {e}"),
                    path: Some(local.to_string_lossy().to_string()),
                })?;
                Ok(())
            }
        }
    }

    /// Upload named top-level files (run_report.json, config_used.json)
    /// that exist in the working dir; no-op for local destinations.
    pub(crate) fn finalize_run(&self, files: &[&str]) -> Result<(), errors::CliError> {
        match &self.dest {
            Destination::Local(_) => Ok(()),
            Destination::Remote { .. } => {
                for bn in files {
                    let local = self.root().join(bn);
                    if !local.exists() {
                        continue;
                    }
                    let dest = self.dest_uri_for(bn);
                    tims_stage::upload_file(&local, &dest).map_err(|e| errors::CliError::Io {
                        source: format!("upload {dest}: {e}"),
                        path: None,
                    })?;
                }
                Ok(())
            }
        }
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

#[cfg(test)]
mod destination_tests {
    use super::OutputSink;

    #[test]
    fn a_local_sink_writes_where_it_says_it_writes() {
        let dir = tempfile::tempdir().unwrap();
        let sink = OutputSink::new(dir.path().to_str().unwrap()).unwrap();

        assert_eq!(sink.dest_root(), sink.root().to_string_lossy());
        assert_eq!(
            sink.dest_uri_for("run"),
            sink.root().join("run").to_string_lossy()
        );
    }

    #[test]
    fn a_remote_sink_reports_the_remote_uri_rather_than_its_tempdir() {
        let sink = OutputSink::new("s3://bkt/out/").unwrap();

        assert_eq!(sink.dest_root(), "s3://bkt/out");
        assert_eq!(sink.dest_uri_for("run"), "s3://bkt/out/run");
        assert!(
            sink.root().is_dir() && !sink.root().starts_with("s3:"),
            "a remote run stages locally first: {:?}",
            sink.root()
        );
    }
}

#[cfg(test)]
mod cleanup_tests {
    use super::OutputSink;
    use crate::artifacts::FEATURE_STATS_TSV;

    #[test]
    fn clearing_a_sample_removes_a_sidecar_this_run_would_not_have_written() {
        let dir = tempfile::tempdir().unwrap();
        let sink = OutputSink::new(dir.path().to_str().unwrap()).unwrap();
        let sample_dir = sink.sample_dir("run");
        std::fs::create_dir_all(&sample_dir).unwrap();
        let sidecar = sample_dir.join(FEATURE_STATS_TSV);
        std::fs::write(&sidecar, "name\tmean\n").unwrap();

        sink.clear_existing("run").unwrap();

        assert!(
            !sidecar.exists(),
            "a stale sidecar next to fresh results reads as though it described them"
        );
    }
}
