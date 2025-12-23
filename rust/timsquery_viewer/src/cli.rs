use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    #[arg(
        long,
        value_name = "FILE",
        help = "Path to raw .d directory",
        short = 'r'
    )]
    pub raw_data_path: Option<PathBuf>,
    #[arg(
        long,
        value_name = "FILE",
        help = "Path to elution groups file (.json or .txt)",
        short = 'e'
    )]
    pub elution_groups_path: Option<PathBuf>,
}
