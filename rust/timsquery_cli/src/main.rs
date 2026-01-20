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

/// Converts verbosity flags to a log level string.
/// Returns the log level based on verbose/quiet counts.
/// If RUST_LOG is set, it takes precedence.
fn get_log_level(verbose: u8, quiet: u8) -> String {
    // RUST_LOG environment variable takes precedence
    if std::env::var("RUST_LOG").is_ok() {
        return std::env::var("RUST_LOG").unwrap();
    }

    // Calculate effective verbosity: positive = more verbose, negative = more quiet
    let effective = verbose as i8 - quiet as i8;

    match effective {
        2.. => "trace".to_string(),
        1 => "debug".to_string(),
        0 => "info".to_string(),
        -1 => "warn".to_string(),
        _ => "error".to_string(),
    }
}

fn main() -> Result<(), CliError> {
    let args = Args::parse();

    let log_level = get_log_level(args.verbose, args.quiet);
    let env_filter = EnvFilter::builder()
        .with_default_directive(log_level.parse().unwrap())
        .from_env_lossy();
    let subscriber = Registry::default()
        .with(env_filter)
        .with(tracing_subscriber::fmt::layer().with_span_events(FmtSpan::CLOSE));

    set_global_default(subscriber).expect("Setting default subscriber failed");

    match args.command {
        Some(Commands::QueryIndex(args)) => main_query_index(args)?,
        Some(Commands::WriteTemplate(args)) => main_write_template(args)?,
        None => {
            println!("No command provided");
        }
    }
    Ok(())
}
