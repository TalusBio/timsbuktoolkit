//! The files a run writes, named once.
//!
//! Every writer, the pre-run collision probe and the overwrite cleanup read
//! these constants, so a filename cannot be spelled one way where it is written
//! and another way where it is looked for.

use crate::errors::CliError;
use crate::output_sink::{
    join_output_uri,
    probe_uri_exists,
    sample_name_from_uri,
};

pub(crate) const RESULTS_PARQUET: &str = "results.parquet";
pub(crate) const PERFORMANCE_REPORT: &str = "performance_report.json";
pub(crate) const FEATURE_STATS_TSV: &str = "results.feature_stats.tsv";
pub(crate) const FEATURE_IMPORTANCE_TSV: &str = "results.feature_importance.tsv";
pub(crate) const RUN_REPORT: &str = "run_report.json";
pub(crate) const CONFIG_USED: &str = "config_used.json";

/// Written once per run, at the root of the output destination.
pub(crate) const RUN_ARTIFACTS: &[&str] = &[RUN_REPORT, CONFIG_USED];

/// Written into a per-sample subdirectory of the output destination. The two
/// feature sidecars appear only when the run is rescoring with feature stats
/// enabled.
///
/// The overwrite cleanup asks for them unconditionally, whatever the current
/// run does: a sidecar from an earlier run describes that run's model, and one
/// left sitting next to fresh results reads as though it described those.
pub(crate) fn per_sample_artifacts(feature_stats: bool) -> Vec<&'static str> {
    let mut artifacts = vec![RESULTS_PARQUET, PERFORMANCE_REPORT];
    if feature_stats {
        artifacts.push(FEATURE_STATS_TSV);
        artifacts.push(FEATURE_IMPORTANCE_TSV);
    }
    artifacts
}

/// Probe every artifact the run is about to write, returning the URIs of the
/// ones already there. Raw inputs rather than sample names, because the sample
/// name is derived from the URI and a URI with no derivable name is an error
/// worth reporting here.
pub(crate) fn probe_collisions(
    output_uri: &str,
    raw_inputs: &[String],
    feature_stats: bool,
) -> Result<Vec<String>, CliError> {
    let mut collisions: Vec<String> = Vec::new();
    for raw_uri in raw_inputs {
        let sample = sample_name_from_uri(raw_uri).ok_or_else(|| CliError::Io {
            source: "Unable to extract file stem".to_string(),
            path: Some(raw_uri.clone()),
        })?;
        for artifact in per_sample_artifacts(feature_stats) {
            let uri = join_output_uri(output_uri, &format!("{sample}/{artifact}"));
            if probe_uri_exists(&uri)? {
                collisions.push(uri);
            }
        }
    }
    for artifact in RUN_ARTIFACTS {
        let uri = join_output_uri(output_uri, artifact);
        if probe_uri_exists(&uri)? {
            collisions.push(uri);
        }
    }
    Ok(collisions)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The sidecars are as much a result as the parquet is, so a rerun that
    /// would write them has to see the ones already on disk.
    #[test]
    fn the_probe_reports_a_feature_stats_sidecar_from_an_earlier_run() {
        let dir = tempfile::tempdir().unwrap();
        let sample_dir = dir.path().join("run");
        std::fs::create_dir_all(&sample_dir).unwrap();
        std::fs::write(sample_dir.join(FEATURE_STATS_TSV), "name\tmean\n").unwrap();

        let raw_inputs = vec!["/data/run.d".to_string()];
        let output_uri = dir.path().to_string_lossy().to_string();

        let collisions = probe_collisions(&output_uri, &raw_inputs, true).unwrap();
        assert_eq!(
            collisions,
            vec![
                sample_dir
                    .join(FEATURE_STATS_TSV)
                    .to_string_lossy()
                    .to_string()
            ]
        );

        assert!(
            probe_collisions(&output_uri, &raw_inputs, false)
                .unwrap()
                .is_empty(),
            "a run that writes no sidecar cannot collide with one"
        );
    }
}
