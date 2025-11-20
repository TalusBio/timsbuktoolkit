mod cli;
mod commands;
mod error;
mod processing;

use clap::Parser;
use tracing::subscriber::set_global_default;
use tracing_subscriber::EnvFilter;
use tracing_subscriber::fmt::format::FmtSpan;
use tracing_subscriber::prelude::*;
use tracing_subscriber::registry::Registry;

use crate::cli::{
    Args,
    Commands,
};
use crate::commands::{
    main_query_index,
    main_write_template,
};
use crate::error::CliError;

// mimalloc seems to work better for windows
// ... more accurately ... not using it causes everything to
// be extremely slow on windows...
#[cfg(target_os = "windows")]
use mimalloc::MiMalloc;

#[cfg(target_os = "windows")]
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

fn main() -> Result<(), CliError> {
    let env_filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info"));
    let subscriber = Registry::default()
        .with(env_filter)
        .with(tracing_subscriber::fmt::layer().with_span_events(FmtSpan::CLOSE));

    set_global_default(subscriber).expect("Setting default subscriber failed");
    let args = Args::parse();

    match args.command {
        Some(Commands::QueryIndex(args)) => main_query_index(args)?,
        Some(Commands::WriteTemplate(args)) => main_write_template(args)?,
        None => {
            println!("No command provided");
        }
    }
    Ok(())
}
