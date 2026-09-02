//! `timsseek build-library`, which predicts a library locally and writes it.
//!
//! This crate owns the flag surface and none of the digestion, modification or
//! prediction logic: every setting is passed straight through to msspeculator,
//! which is also where the output format is decided, from the `--out` suffix.
//!
//! Resolution is the same rule as everywhere else. A flag beats the
//! configuration file's `[library]` section, which beats msspeculator's own
//! default, and "beats" means the flag was given at all.

use std::path::{
    Path,
    PathBuf,
};

use msspeculator_inference::{
    BuiltinModel,
    Difference,
    LibraryCheck,
    LibraryOptions,
    LibraryProvenance,
    ModelSource,
    ProgressFn,
    StreamOptions,
    check_against,
    check_library,
    stream_library,
    write_library,
};
use timsseek::DecoyPolicy;
use tracing::info;

use crate::build_progress::BuildProgress;
use crate::cli::BuildLibraryArgs;
use crate::config::{
    BuildConfig,
    LibraryConfig,
};
use crate::errors::CliError;
use crate::predicted_library::{
    self,
    PredictedLibrary,
};

/// msspeculator's defaults, restated here only so a partially-specified
/// `[library]` section does not have to name every field to change one.
const DEFAULT_MODEL: &str = "builtin:small-v0";
/// The one default that is this project's rather than msspeculator's, whose
/// `--decoys` is off.
///
/// Every library built here exists to be searched, and a search needs a decoy
/// for every target to put an FDR on. Predicting them costs twice the
/// prediction; the alternative costs more, because a library that ships none
/// gets ±CH2 mass-shift decoys derived at load, and a mass shift is a far
/// weaker null than a pseudo-reversed sequence run through the same model.
/// Those exist for libraries nothing can regenerate -- a vendor export -- not
/// for one being predicted right here.
const DEFAULT_DECOYS: bool = true;
const DEFAULT_MISSED_CLEAVAGES: usize = 2;
const DEFAULT_MIN_LENGTH: usize = 7;
const DEFAULT_MAX_LENGTH: usize = 30;
const DEFAULT_MIN_CHARGE: i64 = 2;
const DEFAULT_MAX_CHARGE: i64 = 4;
const DEFAULT_MAX_VARIABLE_MODS: usize = 1;
const DEFAULT_MIN_INTENSITY: f64 = 0.01;
const DEFAULT_FIXED_MOD: &str = "C[UNIMOD:4]";
const DEFAULT_VARIABLE_MOD: &str = "M[UNIMOD:35]";

/// Every setting that decides what a library *contains*, with the command line
/// already resolved against the configuration file. Owned rather than borrowed
/// so the caller holds one value per setting rather than three sources to
/// consult.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedPrediction {
    pub fasta: PathBuf,
    pub model: String,
    pub missed_cleavages: usize,
    pub min_length: usize,
    pub max_length: usize,
    pub min_charge: i64,
    pub max_charge: i64,
    pub fixed_mods: Vec<String>,
    pub variable_mods: Vec<String>,
    pub max_variable_mods: usize,
    pub min_intensity: f64,
    pub max_fragments: Option<usize>,
    pub decoys: bool,
}

/// What to predict, plus where `build-library` puts it.
///
/// Two types rather than one because the output half means nothing to a caller
/// that predicts into memory, and a placeholder `out` is a field whose value no
/// reader could trust.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedBuild {
    pub prediction: ResolvedPrediction,
    pub out: PathBuf,
    pub config_out: Option<PathBuf>,
    pub overwrite: bool,
}

/// Fold the command line over the configuration file's `[library]` section.
pub fn resolve_build(args: &BuildLibraryArgs, config: &BuildConfig) -> ResolvedBuild {
    let library = library_overrides(args, config.library.clone().unwrap_or_default());
    ResolvedBuild {
        prediction: resolve_prediction(args.fasta.clone(), &library),
        out: args.out.clone(),
        config_out: sidecar_path(args),
        // No configuration counterpart: replacing a file is a decision about
        // this invocation, not about how to predict.
        overwrite: args.overwrite,
    }
}

/// The `[library]` section the command line would have written, so "the flag was
/// given at all" is decided in one place and the defaults are applied in
/// another.
fn library_overrides(args: &BuildLibraryArgs, config: LibraryConfig) -> LibraryConfig {
    LibraryConfig {
        model: args.model.clone().or(config.model),
        missed_cleavages: args.missed_cleavages.or(config.missed_cleavages),
        min_length: args.min_length.or(config.min_length),
        max_length: args.max_length.or(config.max_length),
        min_charge: args.min_charge.or(config.min_charge),
        max_charge: args.max_charge.or(config.max_charge),
        fixed_mods: args.fixed_mods().or(config.fixed_mods),
        variable_mods: args.variable_mods.clone().or(config.variable_mods),
        max_variable_mods: args.max_variable_mods.or(config.max_variable_mods),
        min_intensity: args.min_intensity.or(config.min_intensity),
        max_fragments: args.max_fragments.or(config.max_fragments),
        decoys: args.decoys().or(config.decoys),
    }
}

/// Fill in msspeculator's defaults wherever the section said nothing.
///
/// The only place a default is applied, so a `build-library` run and a search
/// that predicts its own library cannot read one `[library]` section two ways.
pub fn resolve_prediction(fasta: PathBuf, library: &LibraryConfig) -> ResolvedPrediction {
    ResolvedPrediction {
        fasta,
        model: library
            .model
            .clone()
            .unwrap_or_else(|| DEFAULT_MODEL.to_string()),
        missed_cleavages: library.missed_cleavages.unwrap_or(DEFAULT_MISSED_CLEAVAGES),
        min_length: library.min_length.unwrap_or(DEFAULT_MIN_LENGTH),
        max_length: library.max_length.unwrap_or(DEFAULT_MAX_LENGTH),
        min_charge: library.min_charge.unwrap_or(DEFAULT_MIN_CHARGE),
        max_charge: library.max_charge.unwrap_or(DEFAULT_MAX_CHARGE),
        fixed_mods: library
            .fixed_mods
            .clone()
            .unwrap_or_else(|| vec![DEFAULT_FIXED_MOD.to_string()]),
        variable_mods: library
            .variable_mods
            .clone()
            .unwrap_or_else(|| vec![DEFAULT_VARIABLE_MOD.to_string()]),
        max_variable_mods: library
            .max_variable_mods
            .unwrap_or(DEFAULT_MAX_VARIABLE_MODS),
        min_intensity: library.min_intensity.unwrap_or(DEFAULT_MIN_INTENSITY),
        max_fragments: library.max_fragments,
        decoys: library.decoys.unwrap_or(DEFAULT_DECOYS),
    }
}

/// The prediction settings for a search that named no library, resolved from the
/// same `[library]` section a `build-library` run reads.
///
/// The decoy policy does not enter here. `[library] decoys` decides what the
/// library carries and the policy decides what happens to it at seal, exactly as
/// for a library on disk. Predicting them is the default, so the usual run seals
/// with real decoys and derives nothing; a `decoys = false` gets mass-shift
/// decoys under `if-missing` and none at all under `never`, which
/// `ReferenceLibrary::try_from` warns about. Reading the policy here
/// instead would override a `decoys = false` the user wrote down.
pub fn resolve_search_prediction(
    fasta: PathBuf,
    library: Option<&LibraryConfig>,
) -> ResolvedPrediction {
    resolve_prediction(fasta, &library.cloned().unwrap_or_default())
}

/// The build a `--build-if-missing` search performs before it opens the library
/// it just wrote.
///
/// The sidecar is written here as it is for `build-library`: a library on disk
/// and no record of the settings that produced it is a library nobody can
/// reproduce, and the two routes writing the same pair of files means a library
/// is the same artifact whichever command made it.
pub(crate) fn resolve_search_build(
    out: PathBuf,
    fasta: PathBuf,
    library: Option<&LibraryConfig>,
    overwrite: bool,
) -> ResolvedBuild {
    ResolvedBuild {
        prediction: resolve_search_prediction(fasta, library),
        config_out: Some(default_sidecar(&out)),
        out,
        overwrite,
    }
}

/// A build reads and writes the filesystem directly, so a remote URI is
/// rejected by name rather than reaching `File::open` and surfacing as "No such
/// file or directory".
///
/// The staging machinery a search uses is right here in this crate, so wiring it
/// through is a small change rather than an impossible one; until then the
/// limitation is stated instead of discovered.
fn reject_remote_paths(resolved: &ResolvedBuild) -> Result<(), CliError> {
    for (flag, path) in [
        ("--fasta", &resolved.prediction.fasta),
        ("--out", &resolved.out),
    ] {
        let uri = path.to_string_lossy();
        if tims_stage::is_remote_uri(&uri) {
            return Err(CliError::Config {
                source: format!(
                    "{flag} {uri} is a remote URI, and build-library reads and writes local \
                     paths only. Build locally and copy, or run a search with --speclib-uri, \
                     which does stage remote inputs."
                ),
            });
        }
    }
    Ok(())
}

/// Refuse to replace an existing library unless asked.
///
/// Checked before the model loads, so a mistyped `--out` costs a second rather
/// than the minutes a proteome takes to predict and then discard.
fn reject_existing_output(resolved: &ResolvedBuild) -> Result<(), CliError> {
    if resolved.overwrite {
        return Ok(());
    }
    let existing: Vec<_> = std::iter::once(&resolved.out)
        .chain(resolved.config_out.as_ref())
        .filter(|path| path.exists())
        .collect();
    if existing.is_empty() {
        Ok(())
    } else {
        Err(CliError::Config {
            source: format!(
                "artifact(s) already exist: {}; pass --overwrite/-O to replace them",
                existing
                    .iter()
                    .map(|path| path.display().to_string())
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
        })
    }
}

/// `--config-out <path>`, `<out>.config.json` by default, `None` when
/// `--no-config-out` was given. A dedicated negative flag rather than an empty
/// string, so "skip the sidecar" cannot be confused with a path.
fn sidecar_path(args: &BuildLibraryArgs) -> Option<PathBuf> {
    if args.no_config_out {
        return None;
    }
    args.config_out
        .clone()
        .or_else(|| Some(default_sidecar(&args.out)))
}

/// The sidecar that belongs to a library, named after the file rather than a
/// stem, because the suffix is what picks the format.
pub(crate) fn default_sidecar(out: &Path) -> PathBuf {
    let mut path = out.to_path_buf().into_os_string();
    path.push(".config.json");
    PathBuf::from(path)
}

/// The only place a [`StreamOptions`] is built.
///
/// One construction site, because a search that predicts its own library has to
/// ask msspeculator for exactly what `build-library` would have written, and two
/// literals agreeing today is not the same as their agreeing next year.
fn stream_options<'a>(
    prediction: &'a ResolvedPrediction,
    model: ModelSource,
    report: &'a ProgressFn<'a>,
) -> StreamOptions<'a> {
    StreamOptions {
        model,
        progress: Some(report),
        fasta: &prediction.fasta,
        // No acquisition, activation or chromatography context: the model
        // artifact's own defaults. Selecting a different one is a decision
        // about the model, and this command chooses a model rather than
        // reconfiguring one.
        activation: None,
        ms_context: None,
        chrom_context: None,
        min_intensity: prediction.min_intensity,
        missed_cleavages: prediction.missed_cleavages,
        min_length: prediction.min_length,
        max_length: prediction.max_length,
        min_charge: prediction.min_charge,
        max_charge: prediction.max_charge,
        fixed_mods: &prediction.fixed_mods,
        variable_mods: &prediction.variable_mods,
        max_variable_mods: prediction.max_variable_mods,
        max_fragments: prediction.max_fragments,
        generate_decoys: prediction.decoys,
    }
}

fn differences_text(differences: &[Difference]) -> String {
    differences
        .iter()
        .map(|difference| format!("\n  {difference}"))
        .collect::<Vec<_>>()
        .join("")
}

/// Compare a file-backed search library with the prediction settings supplied
/// alongside it. Advisory: a difference is visible, but never prevents search.
fn search_library_check(
    library: &Path,
    prediction: &ResolvedPrediction,
) -> Result<LibraryCheck, CliError> {
    let model = parse_model_source(&prediction.model)?;
    let no_progress = |_progress| {};
    let options = stream_options(prediction, model, &no_progress);
    check_library(library, None, &options).map_err(|error| CliError::Config {
        source: format!("checking how {} was built: {error:#}", library.display()),
    })
}

pub(crate) fn check_search_library(
    library: &Path,
    display_name: &str,
    prediction: &ResolvedPrediction,
) -> Result<(), CliError> {
    match search_library_check(library, prediction)? {
        LibraryCheck::Same => {
            info!("Library {display_name} was built with the supplied prediction settings")
        }
        LibraryCheck::Different(differences) => tracing::warn!(
            "Library {display_name} was built with different prediction settings:{}",
            differences_text(&differences),
        ),
        LibraryCheck::Unknown => {
            info!("Library {display_name} carries no comparable msspeculator provenance")
        }
    }
    Ok(())
}

fn report_rebuild_check(library: &Path, sidecar: Option<&Path>, provenance: &LibraryProvenance) {
    match check_against(library, sidecar, &provenance.settings) {
        LibraryCheck::Same => info!(
            "Rebuilding {}; it was built with these settings",
            library.display()
        ),
        LibraryCheck::Different(differences) => tracing::warn!(
            "Overwriting {}, which was built with different settings:{}",
            library.display(),
            differences_text(&differences),
        ),
        LibraryCheck::Unknown => tracing::warn!(
            "Overwriting {}, which carries no comparable msspeculator provenance",
            library.display()
        ),
    }
}

/// Predict a library straight into the arena a search scores, writing nothing.
///
/// The same settings, the same [`StreamOptions`] and the same progress rendering
/// [`run`] uses; all that differs is where the rows land.
pub(crate) fn predict_in_memory(
    prediction: &ResolvedPrediction,
    decoys: DecoyPolicy,
) -> Result<PredictedLibrary, CliError> {
    let model = parse_model_source(&prediction.model)?;

    info!(
        "Predicting a library from {} with {}",
        prediction.fasta.display(),
        prediction.model
    );
    // Held in a binding: the callback borrows it, so a temporary would be
    // dropped before `stream_library` is called.
    let progress = BuildProgress::new();
    let report = progress.callback();
    let (handle, sink) = predicted_library::sink();
    let stats = stream_library(&stream_options(prediction, model, &report), sink).map_err(|e| {
        CliError::LibraryBuild {
            source: format!("predicting from {}: {e:#}", prediction.fasta.display()),
        }
    })?;
    // Before anything else writes: an open bar is a line the next line out would
    // land on the end of.
    progress.finish();

    // The digestion counts exist nowhere else on this path -- there is no
    // sidecar and no library file to inspect afterwards.
    info!(
        "{} proteins -> {} peptides -> {} precursors ({} decoys) -> {} fragments",
        stats.proteins, stats.peptides, stats.precursors, stats.decoys, stats.fragments,
    );
    handle.into_library(&stats, decoys)
}

/// Predict a library and write it, with no network and no server.
pub fn run(resolved: &ResolvedBuild) -> Result<(), CliError> {
    reject_remote_paths(resolved)?;
    reject_existing_output(resolved)?;
    let model = parse_model_source(&resolved.prediction.model)?;

    info!(
        "Predicting a library from {} with {}",
        resolved.prediction.fasta.display(),
        resolved.prediction.model
    );
    // Held in a binding: the callback borrows it, so a temporary would be
    // dropped before `write_library` is called.
    let progress = BuildProgress::new();
    let report = progress.callback();
    let pending_out = pending_path(&resolved.out);
    let pending_config = resolved.config_out.as_deref().map(pending_path);
    let mut cleanup = PendingArtifacts::new(
        std::iter::once(pending_out.clone())
            .chain(pending_config.iter().cloned())
            .collect(),
    );
    let report_rebuild = |provenance: &LibraryProvenance| {
        if resolved.out.exists()
            || resolved
                .config_out
                .as_ref()
                .is_some_and(|path| path.exists())
        {
            report_rebuild_check(&resolved.out, resolved.config_out.as_deref(), provenance);
        }
    };
    let stats = write_library(&LibraryOptions {
        out: &pending_out,
        config_out: pending_config.as_deref(),
        stream: stream_options(&resolved.prediction, model, &report),
        before_writing: Some(&report_rebuild),
    })
    .map_err(|e| CliError::LibraryBuild {
        source: format!("{}: {e:#}", resolved.out.display()),
    })?;
    // Before anything else writes: an open bar is a line the summary below
    // would land on the end of.
    progress.finish();

    if let Some(path) = pending_config.as_deref() {
        rewrite_sidecar_output(path, &resolved.out)?;
    }
    let mut artifacts = vec![(pending_out.as_path(), resolved.out.as_path())];
    if let (Some(pending), Some(final_path)) =
        (pending_config.as_deref(), resolved.config_out.as_deref())
    {
        artifacts.push((pending, final_path));
    }
    publish_artifacts(&artifacts, resolved.overwrite)?;
    cleanup.disarm();

    println!(
        "{} proteins -> {} peptides -> {} precursors ({} decoys) -> {} fragments -> {}",
        stats.proteins,
        stats.peptides,
        stats.precursors,
        stats.decoys,
        stats.fragments,
        resolved.out.display()
    );
    Ok(())
}

fn pending_path(final_path: &Path) -> PathBuf {
    use std::sync::atomic::{
        AtomicU64,
        Ordering,
    };
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let name = final_path.file_name().unwrap_or_default().to_string_lossy();
    final_path.with_file_name(format!(
        ".timsseek-{}-{}-{name}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed),
    ))
}

/// Remove unpublished artifacts on every error path.
struct PendingArtifacts(Vec<PathBuf>);

impl PendingArtifacts {
    fn new(paths: Vec<PathBuf>) -> Self {
        Self(paths)
    }

    fn disarm(&mut self) {
        self.0.clear();
    }
}

impl Drop for PendingArtifacts {
    fn drop(&mut self) {
        for path in self.0.drain(..) {
            let _ = std::fs::remove_file(path);
        }
    }
}

/// Publish all generated files together. Existing artifacts are first moved to
/// same-directory backups, so any failed rename restores the previous pair.
fn publish_artifacts(artifacts: &[(&Path, &Path)], overwrite: bool) -> Result<(), CliError> {
    if !overwrite && let Some((_, path)) = artifacts.iter().find(|(_, path)| path.exists()) {
        return Err(CliError::Config {
            source: format!(
                "{} appeared while the library was building; refusing to overwrite it",
                path.display()
            ),
        });
    }

    let mut backups: Vec<(PathBuf, PathBuf)> = Vec::new();
    if overwrite {
        for (_, final_path) in artifacts {
            if final_path.exists() {
                let backup = pending_path(final_path);
                if let Err(error) = std::fs::rename(final_path, &backup) {
                    restore_backups(&backups);
                    return Err(io_error(final_path, error));
                }
                backups.push((final_path.to_path_buf(), backup));
            }
        }
    }

    let mut published: Vec<PathBuf> = Vec::new();
    for (pending, final_path) in artifacts {
        if let Err(error) = std::fs::rename(pending, final_path) {
            for path in published.iter().rev() {
                let _ = std::fs::remove_file(path);
            }
            restore_backups(&backups);
            return Err(io_error(final_path, error));
        }
        published.push(final_path.to_path_buf());
    }
    for (_, backup) in backups {
        if let Err(error) = std::fs::remove_file(&backup) {
            tracing::warn!(
                "failed to remove replaced-artifact backup {}: {error}",
                backup.display()
            );
        }
    }
    Ok(())
}

fn restore_backups(backups: &[(PathBuf, PathBuf)]) {
    for (final_path, backup) in backups.iter().rev() {
        let _ = std::fs::rename(backup, final_path);
    }
}

fn io_error(path: &Path, error: std::io::Error) -> CliError {
    CliError::Io {
        source: error.to_string(),
        path: Some(path.display().to_string()),
    }
}

fn rewrite_sidecar_output(path: &Path, final_out: &Path) -> Result<(), CliError> {
    let bytes = std::fs::read(path).map_err(|e| CliError::Io {
        source: e.to_string(),
        path: Some(path.display().to_string()),
    })?;
    let mut sidecar: serde_json::Value =
        serde_json::from_slice(&bytes).map_err(|e| CliError::LibraryBuild {
            source: format!("invalid generated sidecar {}: {e}", path.display()),
        })?;
    if let Some(output) = sidecar
        .get_mut("output")
        .and_then(serde_json::Value::as_object_mut)
    {
        output.insert(
            "path".to_string(),
            serde_json::Value::String(final_out.display().to_string()),
        );
    }
    let bytes = serde_json::to_vec_pretty(&sidecar).map_err(|e| CliError::LibraryBuild {
        source: format!("serializing generated sidecar {}: {e}", path.display()),
    })?;
    std::fs::write(path, bytes).map_err(|e| CliError::Io {
        source: e.to_string(),
        path: Some(path.display().to_string()),
    })
}

/// `builtin:NAME` for a model compiled into msspeculator, anything else a path
/// to a `.safetensors` artifact.
///
/// The prefix and the set of names come from msspeculator, so a model added
/// there is accepted here without an edit.
fn parse_model_source(spec: &str) -> Result<ModelSource, CliError> {
    let Some(name) = spec.strip_prefix(msspeculator_core::builtin::BUILTIN_PREFIX) else {
        return Ok(ModelSource::File(PathBuf::from(spec)));
    };
    [BuiltinModel::SmallV0]
        .into_iter()
        .find(|model| model.name() == name)
        .map(ModelSource::Builtin)
        .ok_or_else(|| CliError::Config {
            source: format!(
                "--model: unknown builtin model {name:?}; this build carries: {}",
                msspeculator_core::builtin::names().join(", ")
            ),
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    fn args(argv: &[&str]) -> BuildLibraryArgs {
        #[derive(Parser)]
        struct Wrapper {
            #[command(flatten)]
            args: BuildLibraryArgs,
        }
        Wrapper::parse_from(
            [
                "build-library",
                "--fasta",
                "p.fasta",
                "--out",
                "lib.mzspeclib.txt.gz",
            ]
            .into_iter()
            .chain(argv.iter().copied()),
        )
        .args
    }

    /// A library-only configuration file, with no search sections at all. That
    /// it parses is half of what these tests cover.
    fn build_config(toml: &str) -> BuildConfig {
        toml::from_str(toml).expect("a library-only configuration must parse")
    }

    /// The literals, not the constants. `x.unwrap_or(K) == K` holds whatever `K`
    /// says, so a default drifting from msspeculator's would not show; spelled
    /// out, it shows in the diff.
    ///
    /// `decoys` is the exception and is asserted the other way round, because it
    /// is deliberately not msspeculator's answer.
    #[test]
    fn an_absent_setting_falls_through_to_msspeculators_default() {
        let resolved = resolve_build(&args(&[]), &BuildConfig::default()).prediction;
        assert_eq!(resolved.model, "builtin:small-v0");
        assert_eq!(resolved.missed_cleavages, 2);
        assert_eq!(resolved.min_length, 7);
        assert_eq!(resolved.max_length, 30);
        assert_eq!(resolved.min_charge, 2);
        assert_eq!(resolved.max_charge, 4);
        assert_eq!(resolved.max_variable_mods, 1);
        assert_eq!(resolved.min_intensity, 0.01);
        assert_eq!(resolved.fixed_mods, ["C[UNIMOD:4]"]);
        assert_eq!(resolved.variable_mods, ["M[UNIMOD:35]"]);
        assert_eq!(resolved.max_fragments, None);
        assert!(
            resolved.decoys,
            "predicted decoys beat the mass-shift ones a search would derive"
        );
    }

    #[test]
    fn prediction_streams_directly_into_a_search_library() {
        let fasta = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/test_data/tiny.fasta");
        let mut prediction = resolve_prediction(fasta, &LibraryConfig::default());
        prediction.decoys = false;
        prediction.max_fragments = Some(4);

        let predicted = predict_in_memory(&prediction, DecoyPolicy::IfMissing)
            .expect("builtin model predicts directly into the arena");
        assert!(!predicted.library.is_empty());
        assert!(predicted.provenance.is_object());

        let dir = tempfile::tempdir().unwrap();
        let library = dir.path().join("library.tsv");
        std::fs::write(&library, "PrecursorMz\tProductMz\n100.0\t200.0\n").unwrap();
        std::fs::write(
            default_sidecar(&library),
            serde_json::to_vec(&predicted.provenance).unwrap(),
        )
        .unwrap();
        assert_eq!(
            search_library_check(&library, &prediction).unwrap(),
            LibraryCheck::Same
        );

        let changed = ResolvedPrediction {
            max_fragments: Some(3),
            ..prediction
        };
        let LibraryCheck::Different(differences) =
            search_library_check(&library, &changed).unwrap()
        else {
            panic!("changed prediction settings must be reported");
        };
        assert_eq!(
            differences
                .iter()
                .map(|difference| difference.key.as_str())
                .collect::<Vec<_>>(),
            ["fragments.max_fragments"]
        );
    }

    /// A flag has to be able to turn a configured setting *off*, which a bare
    /// `bool` and a defaulted list could not express: `--decoys` could only ever
    /// add decoys, and an empty `--fixed-mod` list cannot be spelled by
    /// repeating a flag zero times.
    #[test]
    fn a_flag_can_turn_a_configured_setting_off() {
        let config = build_config("[library]\ndecoys = true\nfixed_mods = [\"C[UNIMOD:4]\"]\n");

        let on = resolve_build(&args(&[]), &config).prediction;
        assert!(on.decoys, "the configuration is in force with no flag");
        assert_eq!(on.fixed_mods, ["C[UNIMOD:4]"]);

        let off = resolve_build(&args(&["--no-decoys", "--no-fixed-mods"]), &config).prediction;
        assert!(!off.decoys);
        assert!(
            off.fixed_mods.is_empty(),
            "--no-fixed-mods means none, not the default"
        );
    }

    #[test]
    fn the_configuration_beats_the_default() {
        let config = build_config(
            r#"
[library]
missed_cleavages = 3
decoys = true
fixed_mods = ["C[UNIMOD:4]", "K[UNIMOD:737]"]
"#,
        );
        let resolved = resolve_build(&args(&[]), &config).prediction;
        assert_eq!(resolved.missed_cleavages, 3);
        assert!(resolved.decoys);
        assert_eq!(resolved.fixed_mods.len(), 2);
        // Untouched fields still fall through.
        assert_eq!(resolved.min_length, 7);
    }

    /// One file serves the build and the searches that read what it built, so a
    /// configuration carrying search sections has to be accepted here and its
    /// search half ignored.
    #[test]
    fn search_sections_are_ignored_rather_than_rejected() {
        let config = build_config(
            r#"
[input]
type = "speclib"
uri = "lib.mzspeclib.txt.gz"

[analysis]
chunk_size = 20000

[library]
missed_cleavages = 3
"#,
        );
        assert_eq!(
            resolve_build(&args(&[]), &config)
                .prediction
                .missed_cleavages,
            3
        );
    }

    /// Unknown *top-level* keys are tolerated, but a typo inside `[library]` is
    /// still an error -- that is the half worth keeping strict.
    #[test]
    fn a_typo_inside_the_library_section_is_rejected() {
        let err = toml::from_str::<BuildConfig>(
            r#"
[library]
missed_cleavges = 3
"#,
        )
        .expect_err("a misspelled library field must not be silently ignored")
        .to_string();
        assert!(
            err.contains("missed_cleavges") || err.contains("unknown field"),
            "got: {err}"
        );
    }

    #[test]
    fn a_flag_beats_the_configuration() {
        let config = build_config(
            r#"
[library]
missed_cleavages = 3
model = "builtin:from-config"
"#,
        );
        let resolved = resolve_build(
            &args(&["--missed-cleavages", "1", "--model", "builtin:from-flag"]),
            &config,
        )
        .prediction;
        assert_eq!(resolved.missed_cleavages, 1);
        assert_eq!(resolved.model, "builtin:from-flag");
    }

    /// The suffix decides the format, so the sidecar is named after the file
    /// rather than after a stem the suffix rules would have to agree with.
    #[test]
    fn the_sidecar_is_named_after_the_library() {
        let resolved = resolve_build(&args(&[]), &BuildConfig::default());
        assert_eq!(
            resolved.config_out.as_deref(),
            Some(std::path::Path::new("lib.mzspeclib.txt.gz.config.json"))
        );
    }

    #[test]
    fn no_config_out_skips_the_sidecar() {
        let resolved = resolve_build(&args(&["--no-config-out"]), &BuildConfig::default());
        assert_eq!(resolved.config_out, None);
    }

    /// `ResolvedBuild` for the guards, which read only the two paths.
    fn paths(fasta: &str, out: &str) -> ResolvedBuild {
        let mut resolved = resolve_build(&args(&[]), &BuildConfig::default());
        resolved.prediction.fasta = PathBuf::from(fasta);
        resolved.out = PathBuf::from(out);
        resolved
    }

    /// A remote URI is named, not left to surface as a missing file. A search
    /// stages these; a build does not.
    #[test]
    fn a_remote_path_is_rejected_by_name() {
        for (flag, resolved) in [
            ("--fasta", paths("s3://bkt/p.fasta", "lib.mzspeclib.txt")),
            ("--out", paths("p.fasta", "s3://bkt/lib.mzspeclib.txt")),
        ] {
            let err = reject_remote_paths(&resolved).expect_err("a remote URI must be rejected");
            let msg = err.to_string();
            assert!(msg.contains("remote URI"), "{flag}: {msg}");
            assert!(msg.contains(flag), "names the flag: {msg}");
        }
    }

    /// Predicting a proteome takes minutes, so replacing one by accident costs
    /// the old library and the wait to rebuild it.
    #[test]
    fn an_existing_library_is_not_replaced_without_asking() {
        let dir = tempfile::tempdir().expect("tempdir");
        let out = dir.path().join("lib.mzspeclib.txt.gz");
        std::fs::write(&out, b"existing").expect("write fixture");

        let mut resolved = paths("p.fasta", &out.to_string_lossy());
        let err = reject_existing_output(&resolved)
            .expect_err("an existing library must not be clobbered");
        assert!(err.to_string().contains("--overwrite"), "{err}");

        resolved.overwrite = true;
        assert!(reject_existing_output(&resolved).is_ok());

        // A path that does not exist needs no flag.
        let fresh = paths(
            "p.fasta",
            &dir.path().join("new.mzspeclib.txt.gz").to_string_lossy(),
        );
        assert!(reject_existing_output(&fresh).is_ok());
    }

    #[test]
    fn an_existing_sidecar_is_not_replaced_without_asking() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("lib.mzspeclib.txt.gz");
        let sidecar = default_sidecar(&out);
        std::fs::write(&sidecar, b"existing").unwrap();
        let mut resolved = paths("p.fasta", &out.to_string_lossy());
        resolved.config_out = Some(sidecar);

        assert!(reject_existing_output(&resolved).is_err());
        resolved.overwrite = true;
        assert!(reject_existing_output(&resolved).is_ok());
    }

    #[test]
    fn publishing_one_failed_artifact_restores_the_previous_pair() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("library");
        let sidecar = dir.path().join("sidecar");
        let pending_out = dir.path().join("pending-library");
        let missing_pending_sidecar = dir.path().join("missing-pending-sidecar");
        std::fs::write(&out, b"old library").unwrap();
        std::fs::write(&sidecar, b"old sidecar").unwrap();
        std::fs::write(&pending_out, b"new library").unwrap();

        assert!(
            publish_artifacts(
                &[(&pending_out, &out), (&missing_pending_sidecar, &sidecar),],
                true,
            )
            .is_err()
        );
        assert_eq!(std::fs::read(out).unwrap(), b"old library");
        assert_eq!(std::fs::read(sidecar).unwrap(), b"old sidecar");
    }

    /// The "same targets" criterion, as far as it can be asserted without a
    /// model: both routes resolve one `[library]` section to the same settings.
    ///
    /// The `StreamOptions` half is not asserted, because there is one
    /// construction site for it and comparing two calls to `stream_options`
    /// would restate that rather than test it.
    #[test]
    fn a_search_resolves_the_prediction_settings_a_build_would() {
        let config = build_config(
            r#"
[library]
model = "builtin:small-v0"
missed_cleavages = 3
min_charge = 1
max_fragments = 12
decoys = true
"#,
        );
        let built = resolve_build(&args(&[]), &config).prediction;
        let searched = resolve_search_prediction(PathBuf::from("p.fasta"), config.library.as_ref());
        assert_eq!(built, searched);
    }

    /// A library a search built and one `build-library` wrote are the same pair
    /// of files, produced from the same settings.
    #[test]
    fn a_search_that_builds_a_library_writes_what_build_library_would() {
        let built = resolve_build(&args(&[]), &BuildConfig::default());
        let searched = resolve_search_build(
            built.out.clone(),
            built.prediction.fasta.clone(),
            None,
            false,
        );
        assert_eq!(searched.prediction, built.prediction);
        assert_eq!(
            searched.config_out, built.config_out,
            "a library with no record of the settings that produced it is not reproducible"
        );
    }

    /// What a predicted library carries is `[library] decoys` and nothing else.
    /// Reading the policy here as well would mean a search predicted a different
    /// library than the `build-library` run it is supposed to match, and would
    /// override a `decoys = false` the user wrote down.
    #[test]
    fn the_decoy_policy_does_not_change_what_a_prediction_generates() {
        for configured in [false, true] {
            let library = LibraryConfig {
                decoys: Some(configured),
                ..Default::default()
            };
            let prediction = resolve_search_prediction(PathBuf::from("p.fasta"), Some(&library));
            assert_eq!(prediction.decoys, configured);
        }
    }

    /// A search flag under `build-library` is an error, not an argument that is
    /// accepted and then ignored.
    #[test]
    fn a_search_flag_is_not_accepted_by_build_library() {
        use crate::cli::Cli;
        let err = Cli::try_parse_from(["timsseek", "build-library", "--raw-inputs", "a.d"])
            .expect_err("--raw-inputs is a search flag");
        assert_eq!(err.kind(), clap::error::ErrorKind::UnknownArgument);
    }
}
