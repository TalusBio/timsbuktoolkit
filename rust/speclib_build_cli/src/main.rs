mod cli;
mod config;
mod dedup;
mod decoys;
mod entry;
mod koina;
mod mods;
mod pipeline;

use clap::Parser;
use cli::Cli;
use config::SpeclibBuildConfig;

fn main() {
    let cli = Cli::parse();

    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| "info".into()),
        )
        .with_writer(std::io::stderr)
        .init();

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

    let rt = tokio::runtime::Runtime::new().unwrap();
    if let Err(e) = rt.block_on(pipeline::run(&config)) {
        eprintln!("Pipeline error: {e}");
        std::process::exit(1);
    }
}
