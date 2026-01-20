/// This file is mainly used for testing and demonstration purposes.
///
/// In particular is used to check that the code to generate
/// spectral libs from python works correctly in Rust.
use clap::{
    Parser,
    Subcommand,
};
// use serde_json::to_string_pretty;
use timsseek::Speclib;
use tracing_subscriber::EnvFilter;
use tracing_subscriber::fmt::format::FmtSpan;
use tracing_subscriber::prelude::*;
use tracing_subscriber::registry::Registry;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Increase logging verbosity (can be repeated: -v for debug, -vv for trace)
    #[arg(short, long, action = clap::ArgAction::Count, global = true)]
    verbose: u8,

    /// Decrease logging verbosity (can be repeated: -q for warn, -qq for error)
    #[arg(short, long, action = clap::ArgAction::Count, global = true)]
    quiet: u8,

    #[command(subcommand)]
    command: SubCommands,
}

#[derive(Subcommand)]
enum SubCommands {
    // Sample,
    Parse {
        #[arg(short, long)]
        speclib_file: String,
    },
}

// fn generate_sample_speclib() {
//     // This function generates a sample speclib and prints it to stdout.
//     let speclib = Speclib::sample();
//     println!("{}", to_string_pretty(&speclib).unwrap());
// }

fn parse_speclib(speclib_file: &str) -> Result<(), timsseek::errors::LibraryReadingError> {
    // This function parses a speclib file and prints the parsed content to stdout.
    let speclib = Speclib::from_file(
        std::path::Path::new(speclib_file),
        timsseek::DecoyStrategy::default(),
    )?;
    println!("{:#?}", speclib);
    Ok(())
}

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

fn setup_logger(verbose: u8, quiet: u8) {
    let log_level = get_log_level(verbose, quiet);
    let env_filter = EnvFilter::builder()
        .with_default_directive(log_level.parse().unwrap())
        .from_env_lossy();
    let subscriber = Registry::default()
        .with(env_filter)
        .with(tracing_subscriber::fmt::layer().with_span_events(FmtSpan::CLOSE));
    subscriber.init();
}

fn main() -> Result<(), timsseek::errors::LibraryReadingError> {
    let cli = Cli::parse();
    setup_logger(cli.verbose, cli.quiet);

    match &cli.command {
        SubCommands::Parse { speclib_file } => {
            parse_speclib(speclib_file)?;
            Ok(())
        }
    }
}
