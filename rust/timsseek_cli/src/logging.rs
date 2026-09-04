//! Subscriber installation, and the log-file placement rules that go with it.

use tims_stage::is_remote_uri;
use tracing_subscriber::filter::EnvFilter;
use tracing_subscriber::fmt::format::FmtSpan;
use tracing_subscriber::fmt::{
    self,
};
use tracing_subscriber::prelude::*;

#[cfg(feature = "instrumentation")]
use tracing_profile::{
    PrintTreeConfig,
    PrintTreeLayer,
};

use crate::cli::SearchArgs;
use crate::run_inputs::ResolvedInputs;

/// Holds the `instrumentation` flush guard, which flushes aggregated spans on
/// drop and so must outlive every traced span. Type-erased via `Box<dyn Any>`
/// to avoid naming the `tracing_profile` type outside the feature gate.
pub(crate) struct TracingHandle {
    #[cfg(feature = "instrumentation")]
    _tree_guard: Box<dyn std::any::Any>,
    #[cfg(not(feature = "instrumentation"))]
    _empty: (),
}

/// The filter both subscribers use. `EnvFilter` is not `Clone`, so callers
/// needing one per layer call this repeatedly.
fn env_filter(level: &str) -> EnvFilter {
    EnvFilter::builder()
        .with_default_directive(level.parse().unwrap_or_else(|_| "info".parse().unwrap()))
        .from_env_lossy()
        .add_directive("forust_ml=warn".parse().unwrap())
        .add_directive("timscentroid::storage=warn".parse().unwrap())
}

/// A stderr-only subscriber, for subcommands with no output directory to put a
/// log file in.
///
/// `build-library` predicts for minutes and has nothing else to say meanwhile,
/// so its progress has to reach the terminal. Without a subscriber installed
/// every `info!` on that path is discarded and the command looks hung; the
/// binary this replaced installed one for the same reason.
pub(crate) fn init_stderr_logging() {
    tracing_subscriber::registry()
        .with(
            fmt::layer()
                .without_time()
                .with_writer(std::io::stderr)
                .with_filter(env_filter("info")),
        )
        .init();
}

/// Install the tracing subscriber. Log file resolution order:
///   1. `--log-path -`            → stderr-only, no file
///   2. `--log-path <p>`          → that path verbatim
///   3. default, local output     → `<output_dir>/timsseek-<ts>.log`
///   4. default, no/remote output → `./timsseek-<ts>.log` in CWD
///
/// The `YYYYMMDDTHHMMSS` local-time suffix avoids clobbering previous runs
/// sharing the directory. Logs never reach the terminal unless `--log-path -`
/// opts in.
pub(crate) fn init_tracing(args: &SearchArgs, resolved: &ResolvedInputs) -> TracingHandle {
    let log_file_path: Option<std::path::PathBuf> = match args.log_path {
        Some(ref p) if p.to_str() == Some("-") => None,
        Some(ref p) => Some(p.clone()),
        None => {
            let base: std::path::PathBuf = Some(&resolved.output_uri)
                .filter(|d| !is_remote_uri(d.as_str()))
                .map(std::path::PathBuf::from)
                .unwrap_or_else(|| {
                    std::env::current_dir().unwrap_or_else(|_| std::path::PathBuf::from("."))
                });
            let ts = chrono::Local::now().format("%Y%m%dT%H%M%S");
            Some(base.join(format!("timsseek-{ts}.log")))
        }
    };

    let build_env_filter = || env_filter(&args.log_level);

    let (file_layer, stderr_warn_layer, stderr_all_layer) =
        if let Some(ref log_path) = log_file_path {
            if let Some(parent) = log_path.parent() {
                let _ = std::fs::create_dir_all(parent);
            }
            let log_file = std::fs::File::create(log_path).expect("Failed to create log file");
            let fl = fmt::layer()
                .with_writer(std::sync::Mutex::new(log_file))
                .with_filter(build_env_filter());
            // Warnings, plus the index build's one-line tolerance report: a
            // build parameter the user should see next to "Loading index".
            let sl = fmt::layer()
                .with_writer(std::io::stderr)
                .without_time()
                .with_filter(EnvFilter::new("warn,timscentroid::indexing=info"));
            (Some(fl), Some(sl), None)
        } else {
            let sl = fmt::layer()
                .with_writer(std::io::stderr)
                .with_filter(build_env_filter());
            (None, None, Some(sl))
        };

    let spans_layer = log_file_path.as_ref().map(|log_path| {
        let mut spans_path = log_path.clone();
        let fname = spans_path
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("timsseek.log")
            .to_string();
        spans_path.set_file_name(format!("{fname}.spans.jsonl"));
        let spans_file = std::fs::File::create(&spans_path).expect("Failed to create spans file");
        fmt::layer()
            .json()
            .with_span_events(FmtSpan::NEW | FmtSpan::CLOSE)
            .with_writer(std::sync::Mutex::new(spans_file))
            .with_filter(build_env_filter())
    });

    #[cfg(feature = "instrumentation")]
    let perf_filter = EnvFilter::builder()
        .with_default_directive("trace".parse().unwrap())
        .with_env_var("RUST_PERF_LOG")
        .from_env_lossy()
        .add_directive("forust_ml::gradientbooster=warn".parse().unwrap());

    #[cfg(feature = "instrumentation")]
    let events_filter = tracing_subscriber::filter::filter_fn(|metadata| !metadata.is_event());

    #[cfg(feature = "instrumentation")]
    let (tree_layer, _tree_guard) = PrintTreeLayer::new(PrintTreeConfig {
        attention_above_percent: 25.0,
        relevant_above_percent: 2.5,
        hide_below_percent: 0.0,
        display_unaccounted: true,
        no_color: false,
        accumulate_spans_count: false,
        accumulate_events: false,
        aggregate_similar_siblings: true,
    });
    #[cfg(feature = "instrumentation")]
    let tree_layer = tree_layer
        .with_filter(perf_filter)
        .with_filter(events_filter);

    let reg = tracing_subscriber::registry()
        .with(file_layer)
        .with(spans_layer)
        .with(stderr_warn_layer)
        .with(stderr_all_layer);

    #[cfg(feature = "instrumentation")]
    let reg = reg.with(tree_layer);

    reg.init();

    println!("timsseek v{}", env!("CARGO_PKG_VERSION"));
    match log_file_path {
        Some(ref log_path) => println!("Log: {}", log_path.display()),
        None => println!("Log: <stderr> (--log-path -)"),
    }
    println!();

    TracingHandle {
        #[cfg(feature = "instrumentation")]
        _tree_guard: Box::new(_tree_guard),
        #[cfg(not(feature = "instrumentation"))]
        _empty: (),
    }
}
