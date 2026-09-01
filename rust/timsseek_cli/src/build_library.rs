//! `timsseek build-library`, which predicts a library locally and writes it.
//!
//! This crate owns the flag surface and none of the digestion, modification or
//! prediction logic: every setting is passed straight through to msspeculator,
//! which is also where the output format is decided, from the `--out` suffix.
//!
//! Resolution is the same rule as everywhere else. A flag beats the
//! configuration file's `[library]` section, which beats msspeculator's own
//! default, and "beats" means the flag was given at all.

use std::path::PathBuf;

use msspeculator_inference::{
    BuiltinModel,
    LibraryOptions,
    ModelSource,
    ProgressFn,
    StreamOptions,
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
        decoys: library.decoys.unwrap_or(false),
    }
}

/// The prediction settings for a search that named no library, resolved from the
/// same `[library]` section a `build-library` run reads.
///
/// The decoy policy does not enter here. `[library] decoys` decides what the
/// library carries and the policy decides what happens to it at seal, exactly as
/// for a library on disk: a predicted library carrying none gets mass-shift
/// decoys derived under `if-missing`, has them derived over its own under
/// `force`, and scores without any under `never` -- which
/// `ReferenceLibrary::from_sealed_arena` warns about, because a run with nothing
/// to estimate FDR from is worth saying out loud rather than silently correcting.
/// Deriving the setting from the policy instead would override a `decoys = false`
/// the user wrote down.
pub fn resolve_search_prediction(
    fasta: PathBuf,
    library: Option<&LibraryConfig>,
) -> ResolvedPrediction {
    resolve_prediction(fasta, &library.cloned().unwrap_or_default())
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
    if resolved.overwrite || !resolved.out.exists() {
        return Ok(());
    }
    Err(CliError::Config {
        source: format!(
            "{} already exists; pass --overwrite/-O to replace it",
            resolved.out.display()
        ),
    })
}

/// `--config-out <path>`, `<out>.config.json` by default, `None` when
/// `--no-config-out` was given. A dedicated negative flag rather than an empty
/// string, so "skip the sidecar" cannot be confused with a path.
fn sidecar_path(args: &BuildLibraryArgs) -> Option<PathBuf> {
    if args.no_config_out {
        return None;
    }
    args.config_out.clone().or_else(|| {
        let mut path = args.out.clone().into_os_string();
        path.push(".config.json");
        Some(PathBuf::from(path))
    })
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
    let stats = write_library(&LibraryOptions {
        out: &resolved.out,
        config_out: resolved.config_out.as_deref(),
        stream: stream_options(&resolved.prediction, model, &report),
    })
    .map_err(|e| CliError::LibraryBuild {
        source: format!("{}: {e:#}", resolved.out.display()),
    })?;
    // Before anything else writes: an open bar is a line the summary below
    // would land on the end of.
    progress.finish();

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
        assert!(!resolved.decoys);
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
