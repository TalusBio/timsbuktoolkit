use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    /// Increase logging verbosity (can be repeated: -v for debug, -vv for trace)
    #[arg(short, long, action = clap::ArgAction::Count, global = true)]
    pub verbose: u8,

    /// Decrease logging verbosity (can be repeated: -q for warn, -qq for error)
    #[arg(short, long, action = clap::ArgAction::Count, global = true)]
    pub quiet: u8,

    #[arg(
        long,
        value_name = "FILE",
        help = "Path to raw MS data or a cached index",
        short = 'r'
    )]
    pub raw_data_path: Option<PathBuf>,
    #[arg(
        long,
        value_name = "FILE",
        help = "Path to a spectral library with ion-annotated fragments and reference intensities",
        short = 'e'
    )]
    pub elution_groups_path: Option<PathBuf>,
}
