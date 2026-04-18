use clap::Parser;
use std::path::PathBuf;
use timsseek::DecoyStrategy;

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

    /// Path to the JSON configuration file (optional, uses defaults if not provided)
    #[arg(short, long)]
    pub config: Option<PathBuf>,

    /// Path to the .d file(s) (will over-write the config file)
    /// Can be specified multiple times for batch processing
    #[arg(short, long, value_name = "FILE")]
    pub dotd_files: Vec<PathBuf>,

    /// Path to the speclib file (will over-write the config file)
    #[arg(short, long)]
    pub speclib_file: Option<PathBuf>,

    /// Path to a calibration library (optional).
    /// If provided, Phase 1 prescore uses this library for calibrant selection,
    /// while Phase 3 scoring uses the main speclib.
    /// If not provided, the main speclib is used for both phases.
    #[arg(long)]
    pub calib_lib: Option<PathBuf>,

    /// Path to the output directory
    #[arg(short, long)]
    pub output_dir: Option<PathBuf>,

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
    pub decoy_strategy: Option<DecoyStrategy>,

    /// Print the default TOML configuration to stdout and exit.
    #[arg(long, conflicts_with = "write_default_config")]
    pub print_default_config: bool,

    /// Write the default TOML configuration to the given path and exit.
    #[arg(long, value_name = "PATH")]
    pub write_default_config: Option<PathBuf>,
}
