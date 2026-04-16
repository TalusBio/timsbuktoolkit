mod cli;
mod config;
mod dedup;
mod decoys;
mod koina;
mod mods;

use clap::Parser;
use cli::Cli;
use config::SpeclibBuildConfig;

fn main() {
    let cli = Cli::parse();
    let config = match SpeclibBuildConfig::from_cli(&cli) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Error loading config: {e}");
            std::process::exit(1);
        }
    };
    if let Err(e) = config.validate() {
        eprintln!("Invalid config: {e}");
        std::process::exit(1);
    }
    eprintln!("Config loaded. Pipeline not yet implemented.");
}
