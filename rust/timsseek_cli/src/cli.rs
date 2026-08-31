use clap::{
    Parser,
    ValueEnum,
};
use std::path::PathBuf;
use timsseek::DecoyPolicy;
use timsseek::ml::RescoreModel;

/// Clap mirror of [`timsseek::ml::RescoreModel`]: `ValueEnum` is a foreign
/// trait, so it cannot be implemented for the lib-owned type from this crate.
///
/// `rename_all` matches `RescoreModel`'s serde spelling so model names are
/// identical on the command line and in configuration files.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
#[clap(rename_all = "snake_case")]
pub enum CliRescoreModel {
    Gbm,
    Lda,
    Hybrid,
    Mlp,
}

impl From<CliRescoreModel> for RescoreModel {
    fn from(v: CliRescoreModel) -> Self {
        match v {
            CliRescoreModel::Gbm => RescoreModel::Gbm,
            CliRescoreModel::Lda => RescoreModel::Lda,
            CliRescoreModel::Hybrid => RescoreModel::Hybrid,
            CliRescoreModel::Mlp => RescoreModel::Mlp,
        }
    }
}

#[cfg(test)]
impl From<RescoreModel> for CliRescoreModel {
    fn from(v: RescoreModel) -> Self {
        match v {
            RescoreModel::Gbm => CliRescoreModel::Gbm,
            RescoreModel::Lda => CliRescoreModel::Lda,
            RescoreModel::Hybrid => CliRescoreModel::Hybrid,
            RescoreModel::Mlp => CliRescoreModel::Mlp,
        }
    }
}

/// `timsseek` with no subcommand is a search, so the search arguments are
/// flattened at the top level and the subcommand is optional. Every existing
/// invocation keeps parsing to the same values.
///
/// `args_conflicts_with_subcommands` is what makes the flattened arguments safe
/// alongside a subcommand: two copies of [`SearchArgs`] exist and only one is
/// read, so mixing the two has to be a parse error rather than a silent choice
/// between them.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
#[command(args_conflicts_with_subcommands = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Command>,

    #[command(flatten)]
    pub search: SearchArgs,
}

#[derive(clap::Subcommand, Debug)]
pub enum Command {
    /// Deprecated spelling of a bare invocation; use `timsseek <SEARCH ARGS>`.
    /// Removed after 2026-12-31.
    Search(Box<SearchArgs>),
    /// Predict a spectral library from a FASTA and write it.
    BuildLibrary(Box<BuildLibraryArgs>),
}

/// Subcommand spelling put on a clock, warned about by
/// [`warn_deprecated_spellings`].
pub const DEPRECATED_SUBCOMMANDS: &[(&str, &str)] =
    &[("search", "run timsseek with no subcommand")];

/// The date the deprecated spellings stop being accepted, in the first release
/// after it. Filed as its own issue so the removal is scheduled work rather
/// than a comment nobody reads.
pub const DEPRECATION_REMOVAL_DATE: &str = "2026-12-31";

impl Cli {
    /// The search arguments this invocation resolved to, whether they were
    /// given bare or under the deprecated `search`.
    ///
    /// Reading one and discarding the other is only safe because
    /// `args_conflicts_with_subcommands` makes both-populated unrepresentable.
    pub fn search_args(&self) -> &SearchArgs {
        match &self.command {
            Some(Command::Search(args)) => args,
            _ => &self.search,
        }
    }
}

/// Warn once per deprecated spelling actually typed.
///
/// clap resolves an alias to its canonical name and does not report which
/// spelling produced it, so the raw arguments are what has to be scanned. That
/// covers `--flag value` and `--flag=value` alike, and a subcommand name, which
/// clap likewise reports only as the variant it selected.
///
/// Called before the subscriber is installed, so this writes to stderr rather
/// than going through `tracing`.
pub fn warn_deprecated_spellings<I, S>(argv: I)
where
    I: IntoIterator<Item = S>,
    S: AsRef<str>,
{
    let mut seen: Vec<&str> = Vec::new();
    for arg in argv {
        // `--flag=value` and `--flag value` both have to match, so compare the
        // part left of any `=`.
        let arg = arg.as_ref();
        let spelling = arg.split('=').next().unwrap_or(arg);
        let Some((spelling, replacement)) = DEPRECATED_SUBCOMMANDS
            .iter()
            .find(|(deprecated, _)| *deprecated == spelling)
        else {
            continue;
        };
        if seen.contains(spelling) {
            continue;
        }
        seen.push(spelling);
        eprintln!(
            "warning: `{spelling}` is deprecated and will be removed after \
             {DEPRECATION_REMOVAL_DATE}; {replacement} instead"
        );
    }
}

#[derive(clap::Args, Debug)]
pub struct SearchArgs {
    /// Path to the log file.
    /// Defaults to {output_dir}/timsseek.log.
    /// Use "-" to send logs to stderr instead of a file.
    #[arg(long, value_name = "PATH")]
    pub log_path: Option<PathBuf>,

    /// Log level for the log file (default: info)
    #[arg(long, value_name = "LEVEL", default_value = "info")]
    pub log_level: String,

    /// Path to the configuration file, TOML or JSON (sniffed by extension). Optional; defaults are used if omitted. Use `--write-default-config` to scaffold one.
    #[arg(short, long)]
    pub config: Option<PathBuf>,

    /// Path to the .d file(s) (will over-write the config file)
    /// Can be specified multiple times for batch processing
    #[arg(long = "raw-inputs", alias = "dotd-files", value_name = "URI")]
    pub raw_inputs: Option<Vec<String>>,

    /// Path to the speclib file (will over-write the config file)
    #[arg(long = "speclib-uri", alias = "speclib-file")]
    pub speclib_uri: Option<String>,

    /// URI of a calibration library (optional).
    /// If provided, Phase 1 prescore uses this library for calibrant selection,
    /// while Phase 3 scoring uses the main speclib.
    /// If not provided, the main speclib is used for both phases.
    /// Accepts local paths or s3:// URIs.
    #[arg(long, value_name = "URI")]
    pub calib_lib: Option<String>,

    /// Path to the output directory
    #[arg(long = "output-uri", short = 'o', alias = "output-dir")]
    pub output_uri: Option<String>,

    /// Overwrite existing output directory if it exists
    #[arg(short = 'O', long)]
    pub overwrite: bool,

    /// Maximum q-value for output. Only results at or below this
    /// threshold are written to the Parquet file.
    /// Default: 0.5
    #[arg(long, default_value = "0.5")]
    pub max_qvalue: f32,

    /// Decoy generation strategy
    /// Options: if-missing (default), force, never
    /// - if-missing: Generate mass-shift decoys only if library has none
    /// - force: Drop library decoys and regenerate mass-shift decoys
    /// - never: Use library as-is without decoy generation
    #[arg(long, value_name = "STRATEGY")]
    pub decoy_strategy: Option<DecoyPolicy>,

    /// Rescore model; overrides the config file.
    #[arg(long, value_enum)]
    pub rescore_model: Option<CliRescoreModel>,

    /// Print the default TOML configuration to stdout and exit.
    #[arg(long, conflicts_with = "write_default_config")]
    pub print_default_config: bool,

    /// Write the default TOML configuration to the given path and exit.
    #[arg(long, value_name = "PATH")]
    pub write_default_config: Option<PathBuf>,

    /// Skip writing the post-rescore sidecars (results.feature_stats.tsv and
    /// results.feature_importance.tsv).
    #[arg(long)]
    pub no_feature_stats: bool,
}

/// Prediction settings, mirroring msspeculator's own library command rather
/// than inventing a second vocabulary for the same settings.
///
/// Every field is `Option`, because "the flag was given" is what makes a flag
/// beat the configuration file. The fallbacks live in
/// [`LibraryConfig`](crate::config::LibraryConfig), which in turn falls back to
/// msspeculator's defaults.
#[derive(clap::Args, Debug)]
pub struct BuildLibraryArgs {
    /// Sequence database to digest.
    #[arg(long, value_name = "PATH")]
    pub fasta: PathBuf,

    /// Where to write the library. The suffix picks the format:
    /// `.mzspeclib.txt` writes mzSpecLib, anything else writes DIA-NN TSV, and
    /// a trailing `.gz` compresses either. `.mzspeclib.txt.gz` is the
    /// recommended spelling.
    #[arg(long, short = 'o', value_name = "PATH")]
    pub out: PathBuf,

    /// Path to the configuration file, TOML or JSON. Its `[library]` section
    /// supplies anything not given as a flag.
    #[arg(long, short = 'c', value_name = "PATH")]
    pub config: Option<PathBuf>,

    /// A `.safetensors` artifact, or `builtin:NAME` for one compiled in.
    #[arg(long, value_name = "MODEL")]
    pub model: Option<String>,

    #[arg(long, value_name = "N")]
    pub missed_cleavages: Option<usize>,
    #[arg(long, value_name = "N")]
    pub min_length: Option<usize>,
    #[arg(long, value_name = "N")]
    pub max_length: Option<usize>,
    #[arg(long, value_name = "Z")]
    pub min_charge: Option<i64>,
    #[arg(long, value_name = "Z")]
    pub max_charge: Option<i64>,

    /// Fixed modification rule, repeatable.
    #[arg(long = "fixed-mod", value_name = "TARGETS[MOD]")]
    fixed_mods: Option<Vec<String>>,
    /// Predict with no fixed modifications at all.
    ///
    /// Needed because the default is carbamidomethyl and an empty
    /// `--fixed-mod` list cannot be spelled: repeating a flag zero times is
    /// indistinguishable from not passing it.
    #[arg(long, conflicts_with = "fixed_mods")]
    pub no_fixed_mods: bool,
    /// Variable modification rule, repeatable.
    #[arg(long = "variable-mod", value_name = "TARGETS[MOD]")]
    pub variable_mods: Option<Vec<String>>,
    #[arg(long, value_name = "N")]
    pub max_variable_mods: Option<usize>,

    // No acquisition or chromatography context. The model artifact has its own
    // default for both, and selecting a different one is a decision about the
    // model rather than about the library this command writes -- so it belongs
    // wherever the model is chosen, not on a second set of flags here.
    /// Drop fragments below this base-peak-relative intensity.
    #[arg(long, value_name = "FRACTION")]
    pub min_intensity: Option<f64>,
    /// Keep at most this many of the strongest fragments per precursor, applied
    /// after `--min-intensity`.
    #[arg(long, value_name = "N")]
    pub max_fragments: Option<usize>,

    /// Add pseudo-reversed decoy precursors. `--no-decoys` turns them off.
    ///
    /// Two flags on one field, so "not given" stays distinct from "off": a bare
    /// `bool` could not disable a configured `decoys = true`, which left the
    /// configuration file with the last word on a flag.
    #[arg(long, overrides_with = "no_decoys")]
    decoys: bool,
    #[arg(long, overrides_with = "decoys")]
    no_decoys: bool,

    /// Replace the library if it already exists.
    ///
    /// Predicting a proteome takes minutes, so overwriting one by accident is
    /// expensive in both directions: the old library is gone and the new one
    /// costs the same wait again.
    #[arg(short = 'O', long)]
    pub overwrite: bool,

    /// Where to write the resolved-configuration sidecar. Defaults to
    /// `<out>.config.json`.
    #[arg(long, value_name = "PATH", conflicts_with = "no_config_out")]
    pub config_out: Option<PathBuf>,
    /// Skip the resolved-configuration sidecar.
    #[arg(long)]
    pub no_config_out: bool,
}

impl BuildLibraryArgs {
    /// `Some(true)`/`Some(false)` when either decoy flag was given, `None` when
    /// neither was -- which is what lets the configuration file supply the value
    /// only in the absence of a flag.
    ///
    /// `overrides_with` makes the later flag win rather than erroring, so
    /// `--decoys --no-decoys` is `Some(false)` and both being set at once is
    /// unreachable.
    pub fn decoys(&self) -> Option<bool> {
        match (self.decoys, self.no_decoys) {
            (false, false) => None,
            (decoys, _) => Some(decoys),
        }
    }

    /// The fixed-modification rules asked for, where an empty list is a real
    /// answer and `None` means "not given".
    pub fn fixed_mods(&self) -> Option<Vec<String>> {
        if self.no_fixed_mods {
            return Some(Vec::new());
        }
        self.fixed_mods.clone()
    }
}
