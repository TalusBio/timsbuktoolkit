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
    MsContext,
    StreamOptions,
    write_library,
};
use tracing::info;

use crate::cli::BuildLibraryArgs;
use crate::config::{
    BuildConfig,
    LibraryConfig,
};
use crate::errors::CliError;

/// What `ModelSource::spec()` prefixes a compiled-in model with.
const BUILTIN_PREFIX: &str = "builtin:";

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

/// Every prediction setting, with the command line already resolved against the
/// configuration file. Owned rather than borrowed so the caller holds one value
/// per setting rather than three sources to consult.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedBuild {
    pub fasta: PathBuf,
    pub out: PathBuf,
    pub config_out: Option<PathBuf>,
    pub model: String,
    pub missed_cleavages: usize,
    pub min_length: usize,
    pub max_length: usize,
    pub min_charge: i64,
    pub max_charge: i64,
    pub fixed_mods: Vec<String>,
    pub variable_mods: Vec<String>,
    pub max_variable_mods: usize,
    pub ms_context: Option<String>,
    pub nce: Option<f32>,
    pub chrom_context: Option<String>,
    pub min_intensity: f64,
    pub max_fragments: Option<usize>,
    pub decoys: bool,
}

/// Fold the command line over the configuration file's `[library]` section.
pub fn resolve_build(args: &BuildLibraryArgs, config: &BuildConfig) -> ResolvedBuild {
    let library = config.library.clone().unwrap_or_default();
    let LibraryConfig {
        model,
        missed_cleavages,
        min_length,
        max_length,
        min_charge,
        max_charge,
        fixed_mods,
        variable_mods,
        max_variable_mods,
        ms_context,
        nce,
        chrom_context,
        min_intensity,
        max_fragments,
        decoys,
    } = library;

    ResolvedBuild {
        fasta: args.fasta.clone(),
        out: args.out.clone(),
        config_out: sidecar_path(args),
        model: args
            .model
            .clone()
            .or(model)
            .unwrap_or_else(|| DEFAULT_MODEL.to_string()),
        missed_cleavages: args
            .missed_cleavages
            .or(missed_cleavages)
            .unwrap_or(DEFAULT_MISSED_CLEAVAGES),
        min_length: args.min_length.or(min_length).unwrap_or(DEFAULT_MIN_LENGTH),
        max_length: args.max_length.or(max_length).unwrap_or(DEFAULT_MAX_LENGTH),
        min_charge: args.min_charge.or(min_charge).unwrap_or(DEFAULT_MIN_CHARGE),
        max_charge: args.max_charge.or(max_charge).unwrap_or(DEFAULT_MAX_CHARGE),
        fixed_mods: args
            .fixed_mods
            .clone()
            .or(fixed_mods)
            .unwrap_or_else(|| vec![DEFAULT_FIXED_MOD.to_string()]),
        variable_mods: args
            .variable_mods
            .clone()
            .or(variable_mods)
            .unwrap_or_else(|| vec![DEFAULT_VARIABLE_MOD.to_string()]),
        max_variable_mods: args
            .max_variable_mods
            .or(max_variable_mods)
            .unwrap_or(DEFAULT_MAX_VARIABLE_MODS),
        ms_context: args.ms_context.clone().or(ms_context),
        nce: args.nce.or(nce),
        chrom_context: args.chrom_context.clone().or(chrom_context),
        min_intensity: args
            .min_intensity
            .or(min_intensity)
            .unwrap_or(DEFAULT_MIN_INTENSITY),
        max_fragments: args.max_fragments.or(max_fragments),
        // A boolean flag cannot say "not given", so the configuration is what
        // turns decoys on when the flag is absent.
        decoys: args.decoys || decoys.unwrap_or(false),
    }
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

/// Predict a library and write it, with no network and no server.
pub fn run(resolved: &ResolvedBuild) -> Result<(), CliError> {
    let model = parse_model_source(&resolved.model)?;
    let ms_context = resolved
        .ms_context
        .as_deref()
        .map(parse_ms_context)
        .transpose()?
        .or_else(|| {
            resolved.nce.map(|energy| MsContext::Factors {
                instrument: String::new(),
                detector: String::new(),
                fragmentation: String::new(),
                energy: Some(energy),
            })
        });

    info!(
        "Predicting a library from {} with {}",
        resolved.fasta.display(),
        resolved.model
    );
    let stats = write_library(&LibraryOptions {
        out: &resolved.out,
        config_out: resolved.config_out.as_deref(),
        stream: StreamOptions {
            model,
            fasta: &resolved.fasta,
            activation: None,
            ms_context: ms_context.as_ref(),
            chrom_context: resolved.chrom_context.as_deref(),
            min_intensity: resolved.min_intensity,
            missed_cleavages: resolved.missed_cleavages,
            min_length: resolved.min_length,
            max_length: resolved.max_length,
            min_charge: resolved.min_charge,
            max_charge: resolved.max_charge,
            fixed_mods: &resolved.fixed_mods,
            variable_mods: &resolved.variable_mods,
            max_variable_mods: resolved.max_variable_mods,
            max_fragments: resolved.max_fragments,
            generate_decoys: resolved.decoys,
        },
    })
    .map_err(|e| CliError::DataReading {
        source: format!("building library {}: {e:#}", resolved.out.display()),
    })?;

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
fn parse_model_source(spec: &str) -> Result<ModelSource, CliError> {
    match spec.strip_prefix(BUILTIN_PREFIX) {
        Some("small-v0") => Ok(ModelSource::Builtin(BuiltinModel::SmallV0)),
        Some(name) => Err(CliError::Config {
            source: format!("--model: unknown builtin model {name:?}"),
        }),
        None => Ok(ModelSource::File(PathBuf::from(spec))),
    }
}

/// A named setup fitted into the artifact, or the four acquisition factors
/// spelled out. The `::` separator is what tells them apart: a setup name is a
/// label and nothing about it parses as a factor list.
fn parse_ms_context(spec: &str) -> Result<MsContext, CliError> {
    let invalid = |detail: &str| CliError::Config {
        source: format!("--ms-context {spec:?}: {detail}"),
    };
    if !spec.contains("::") {
        if spec.trim().is_empty() {
            return Err(invalid(
                "expected a setup name or INSTRUMENT::DETECTOR::FRAGMENTATION::ENERGY",
            ));
        }
        return Ok(MsContext::Named(spec.to_string()));
    }
    let parts: Vec<&str> = spec.split("::").collect();
    let [instrument, detector, fragmentation, energy] = parts.as_slice() else {
        return Err(invalid(
            "expected INSTRUMENT::DETECTOR::FRAGMENTATION::ENERGY",
        ));
    };
    Ok(MsContext::Factors {
        instrument: (*instrument).to_string(),
        detector: (*detector).to_string(),
        fragmentation: (*fragmentation).to_string(),
        energy: Some(
            energy
                .parse()
                .map_err(|_| invalid(&format!("energy {energy:?} is not a number")))?,
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
        let resolved = resolve_build(&args(&[]), &BuildConfig::default());
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
        let resolved = resolve_build(&args(&[]), &config);
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
        assert_eq!(resolve_build(&args(&[]), &config).missed_cleavages, 3);
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
        );
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
