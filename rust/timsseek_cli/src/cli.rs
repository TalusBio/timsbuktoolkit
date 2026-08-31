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
/// Nothing is swallowed: every argument is a flag, so no positional can absorb
/// a stray token. An unknown flag is an error and so is an unrecognised first
/// word. `nothing_is_swallowed` pins that, because it holds only while that
/// stays true.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Command>,

    #[command(flatten)]
    pub search: SearchArgs,
}

#[derive(clap::Subcommand, Debug)]
pub enum Command {
    /// Search raw files against a library. The same thing a bare invocation
    /// does, named.
    Search(Box<SearchArgs>),
    /// Predict a spectral library from a FASTA and write it.
    BuildLibrary(Box<BuildLibraryArgs>),
}

impl Cli {
    /// The search arguments this invocation resolved to, whether they were
    /// given bare or under `search`.
    pub fn search_args(&self) -> &SearchArgs {
        match &self.command {
            Some(Command::Search(args)) => args,
            _ => &self.search,
        }
    }
}

#[derive(clap::Args, Debug, Default)]
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
#[derive(clap::Args, Debug, Default)]
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
    pub fixed_mods: Option<Vec<String>>,
    /// Variable modification rule, repeatable.
    #[arg(long = "variable-mod", value_name = "TARGETS[MOD]")]
    pub variable_mods: Option<Vec<String>>,
    #[arg(long, value_name = "N")]
    pub max_variable_mods: Option<usize>,

    /// Acquisition context: a named setup, or the factors spelled out as
    /// "INSTRUMENT::DETECTOR::FRAGMENTATION::ENERGY".
    #[arg(long, value_name = "CONTEXT", conflicts_with = "nce")]
    pub ms_context: Option<String>,
    /// Collision energy only, for an otherwise unknown setup.
    #[arg(long, value_name = "NCE", conflicts_with = "ms_context")]
    pub nce: Option<f32>,
    /// Named chromatography context for a real retention time. Absent means the
    /// context-free index.
    #[arg(long, value_name = "DATASET")]
    pub chrom_context: Option<String>,

    /// Drop fragments below this base-peak-relative intensity.
    #[arg(long, value_name = "FRACTION")]
    pub min_intensity: Option<f64>,
    /// Keep at most this many of the strongest fragments per precursor, applied
    /// after `--min-intensity`.
    #[arg(long, value_name = "N")]
    pub max_fragments: Option<usize>,

    /// Add pseudo-reversed decoy precursors.
    #[arg(long)]
    pub decoys: bool,

    /// Where to write the resolved-configuration sidecar. Defaults to
    /// `<out>.config.json`.
    #[arg(long, value_name = "PATH", conflicts_with = "no_config_out")]
    pub config_out: Option<PathBuf>,
    /// Skip the resolved-configuration sidecar.
    #[arg(long)]
    pub no_config_out: bool,
}
