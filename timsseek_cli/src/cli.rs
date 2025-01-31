use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    /// Path to the JSON configuration file
    #[arg(short, long)]
    pub config: PathBuf,

    /// Path to the .d file (will over-write the config file)
    #[arg(short, long)]
    pub dotd_file: Option<PathBuf>,

    /// Path to the speclib file (will over-write the config file)
    #[arg(short, long)]
    pub speclib_file: Option<PathBuf>,

    /// Path to the output directory
    #[arg(short, long)]
    pub output_dir: Option<PathBuf>,

    /// Whether to enable the full output mode
    #[arg(short, long)]
    pub full_output: bool,
}
