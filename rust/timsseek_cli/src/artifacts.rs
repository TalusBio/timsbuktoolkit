//! The files a run writes, named once.
//!
//! Every writer, the pre-run collision probe and the overwrite cleanup read
//! these constants, so a filename cannot be spelled one way where it is written
//! and another way where it is looked for.

use std::path::{
    Path,
    PathBuf,
};

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
/// feature sidecars appear only when feature stats are enabled.
///
/// The collision probe and the overwrite cleanup both read the whole list,
/// whatever the current run is configured to write: a sidecar from an earlier
/// run describes that run's model, and one left sitting next to fresh results
/// reads as though it described those.
pub(crate) const PER_SAMPLE_ARTIFACTS: &[&str] = &[
    RESULTS_PARQUET,
    PERFORMANCE_REPORT,
    FEATURE_STATS_TSV,
    FEATURE_IMPORTANCE_TSV,
];

/// The suffix a library predicted by a search is written under: mzSpecLib text,
/// gzipped. msspeculator picks the format from the suffix, and this is the
/// spelling `build-library` recommends.
pub(crate) const BUILT_LIBRARY_SUFFIX: &str = ".mzspeclib.txt.gz";

/// Where `--build-if-missing` writes when it was given no URI to write to:
/// `<output>/<fasta stem>.mzspeclib.txt.gz`.
///
/// The stem is the file name with one extension off, and a `.gz` before it, so
/// `proteome.fasta` and `proteome.fasta.gz` name the same library. Not
/// [`sample_name_from_uri`], which strips the `.d`/`.tar`/`.idx` a raw input
/// carries and would leave a FASTA's own suffix in front of `.mzspeclib`.
pub(crate) fn built_library_path(output_uri: &str, fasta: &Path) -> PathBuf {
    let name = fasta
        .file_name()
        .map(|name| name.to_string_lossy())
        .unwrap_or_default();
    let stem = name.strip_suffix(".gz").unwrap_or(&name);
    let stem = stem.rsplit_once('.').map_or(stem, |(stem, _)| stem);
    Path::new(output_uri).join(format!("{stem}{BUILT_LIBRARY_SUFFIX}"))
}

/// Probe every artifact a run can write, returning the URIs of the ones already
/// there. Raw inputs rather than sample names, because the sample name is
/// derived from the URI and a URI with no derivable name is an error worth
/// reporting here.
///
/// `built_library` is the library the run is about to predict, which is an
/// output of the run like any other: one already sitting there cost minutes to
/// predict and describes whatever FASTA and settings produced it.
pub(crate) fn probe_collisions(
    output_uri: &str,
    raw_inputs: &[String],
    built_library: Option<&Path>,
) -> Result<Vec<String>, CliError> {
    let mut collisions: Vec<String> = Vec::new();
    if let Some(library) = built_library {
        let uri = library.to_string_lossy().into_owned();
        if probe_uri_exists(&uri)? {
            collisions.push(uri);
        }
    }
    for raw_uri in raw_inputs {
        let sample = sample_name_from_uri(raw_uri).ok_or_else(|| CliError::Io {
            source: "Unable to extract file stem".to_string(),
            path: Some(raw_uri.clone()),
        })?;
        for artifact in PER_SAMPLE_ARTIFACTS {
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

    /// The sidecars are as much a result as the parquet is, and a run that
    /// writes no sidecar of its own leaves the earlier one describing results
    /// it never produced, so the probe has to report it either way.
    #[test]
    fn the_probe_reports_a_feature_stats_sidecar_from_an_earlier_run() {
        let dir = tempfile::tempdir().unwrap();
        let sample_dir = dir.path().join("run");
        std::fs::create_dir_all(&sample_dir).unwrap();
        std::fs::write(sample_dir.join(FEATURE_STATS_TSV), "name\tmean\n").unwrap();

        let raw_inputs = vec!["/data/run.d".to_string()];
        let output_uri = dir.path().to_string_lossy().to_string();

        let collisions = probe_collisions(&output_uri, &raw_inputs, None).unwrap();
        assert_eq!(
            collisions,
            vec![
                sample_dir
                    .join(FEATURE_STATS_TSV)
                    .to_string_lossy()
                    .to_string()
            ]
        );
    }

    /// A predicted library costs minutes, so a run that would write over one
    /// says so instead of replacing it in silence.
    #[test]
    fn the_probe_reports_a_library_an_earlier_run_predicted() {
        let dir = tempfile::tempdir().unwrap();
        let library = dir.path().join("proteome.mzspeclib.txt.gz");
        std::fs::write(&library, b"library").unwrap();

        let collisions = probe_collisions(
            &dir.path().to_string_lossy(),
            &["/data/run.d".to_string()],
            Some(&library),
        )
        .unwrap();
        assert_eq!(collisions, vec![library.to_string_lossy().to_string()]);
    }

    /// The name a run derives when it was told to build a library but not where
    /// to put it, for each way a FASTA is spelled.
    #[test]
    fn a_derived_library_is_named_after_the_fasta_it_is_predicted_from() {
        for fasta in ["proteome.fasta", "proteome.fasta.gz", "proteome.fa"] {
            assert_eq!(
                built_library_path("out", Path::new(fasta)),
                Path::new("out").join("proteome.mzspeclib.txt.gz"),
                "{fasta}"
            );
        }
        assert_eq!(
            built_library_path("out", Path::new("/seq/human.2024.fasta")),
            Path::new("out").join("human.2024.mzspeclib.txt.gz"),
            "only the last extension comes off"
        );
    }
}
