use clap::Parser;
use std::path::PathBuf;
use timsseek::DecoyStrategy;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    /// Increase logging verbosity (can be repeated: -v for debug, -vv for trace)
    #[arg(short, long, action = clap::ArgAction::Count, global = true)]
    pub verbose: u8,

    /// Decrease logging verbosity (can be repeated: -q for warn, -qq for error)
    #[arg(short, long, action = clap::ArgAction::Count, global = true)]
    pub quiet: u8,

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

    /// Path to the output directory
    #[arg(short, long)]
    pub output_dir: Option<PathBuf>,

    /// Overwrite existing output directory if it exists
    #[arg(short = 'O', long)]
    pub overwrite: bool,

    /// Decoy generation strategy
    /// Options: if-missing (default), force, never
    /// - if-missing: Generate mass-shift decoys only if library has none
    /// - force: Drop library decoys and regenerate mass-shift decoys
    /// - never: Use library as-is without decoy generation
    #[arg(long, value_name = "STRATEGY")]
    pub decoy_strategy: Option<DecoyStrategy>,
}
