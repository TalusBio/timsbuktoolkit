use clap::{
    Parser,
    ValueEnum,
};
use std::path::PathBuf;
use timsseek::DecoyPolicy;
use timsseek::ml::RescoreModel;

/// Clap mirror of [`timsseek::ml::RescoreModel`].
///
/// The model enum itself lives in the library, next to the rescorers it selects
/// between. `clap::ValueEnum` is a foreign trait, so it cannot be implemented
/// for the lib-owned type from this crate — hence this mirror, which exists
/// ONLY to render `[possible values: gbm, lda, hybrid]` and parse the flag. The
/// `From` impl below is exhaustive on purpose: a new variant in the library
/// fails to compile here instead of silently going unreachable from the CLI.
///
/// The TOML `rescore_model` field deserializes straight into the lib type; this
/// mirror is not in that path.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
#[clap(rename_all = "lowercase")]
pub enum CliRescoreModel {
    Gbm,
    Lda,
    Hybrid,
}

impl From<CliRescoreModel> for RescoreModel {
    fn from(v: CliRescoreModel) -> Self {
        match v {
            CliRescoreModel::Gbm => RescoreModel::Gbm,
            CliRescoreModel::Lda => RescoreModel::Lda,
            CliRescoreModel::Hybrid => RescoreModel::Hybrid,
        }
    }
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
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
    pub raw_inputs: Vec<String>,

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

    /// Skip writing results.feature_stats.json sidecar after rescoring.
    #[arg(long)]
    pub no_feature_stats: bool,
}
