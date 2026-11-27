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
    if let Err(e) = build_global_pool() {
        eprintln!("failed to build the global Rayon pool: {e}");
        std::process::exit(1);
    }
    if let Err(e) = run() {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

/// Builds the process-wide Rayon pool before anything can use the default one.
/// `RAYON_NUM_THREADS` keeps Rayon's meaning: unset, unparsable or 0 selects
/// the available parallelism.
fn build_global_pool() -> Result<(), rayon::ThreadPoolBuildError> {
    let requested = std::env::var("RAYON_NUM_THREADS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .filter(|&count| count > 0);
    let thread_count = match requested {
        Some(count) => count,
        None => std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get),
    };
    rayon::ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .thread_name(|index| format!("timsseek-w{index}"))
        .build_global()?;
    assert_eq!(rayon::current_num_threads(), thread_count);
    Ok(())
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
