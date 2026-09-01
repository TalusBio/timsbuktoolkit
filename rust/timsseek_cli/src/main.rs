mod artifacts;
mod build_library;
mod build_progress;
mod cli;
mod config;
#[cfg(feature = "dashboard")]
mod dashboard;
mod errors;
mod logging;
mod output_sink;
mod predicted_library;
mod processing;
mod run_inputs;
mod search;

use clap::Parser;

use cli::{
    Cli,
    Command,
};

#[cfg(all(
    any(target_os = "windows", target_env = "musl"),
    not(feature = "track-alloc")
))]
use mimalloc::MiMalloc;

#[cfg(all(
    any(target_os = "windows", target_env = "musl"),
    not(feature = "track-alloc")
))]
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

#[cfg(feature = "track-alloc")]
#[global_allocator]
static GLOBAL: alloc_track::TrackingAllocator = alloc_track::TrackingAllocator::new();

fn main() {
    if let Err(e) = run() {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

fn run() -> std::result::Result<(), errors::CliError> {
    // Before parsing: clap resolves a deprecated spelling to its canonical name
    // and does not report which one was typed.
    cli::warn_deprecated_spellings(std::env::args().skip(1));

    let cli = Cli::parse();
    if let Some(Command::BuildLibrary(build)) = &cli.command {
        // Ahead of any work, so prediction is not silent. Search installs its
        // own richer subscriber, which needs the resolved output directory to
        // place the log file and so cannot be hoisted here.
        logging::init_stderr_logging();
        let config = match build.config.as_deref() {
            Some(path) => config::parse_config_file(path)?,
            None => config::BuildConfig::default(),
        };
        return build_library::run(&build_library::resolve_build(build, &config));
    }
    search::search(cli.search_args())
}
