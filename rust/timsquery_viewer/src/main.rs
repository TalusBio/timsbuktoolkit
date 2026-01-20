mod app;
mod chromatogram_processor;
mod cli;
mod computed_state;
mod domain;
mod error;
mod file_loader;
mod plot_renderer;
mod ui;

use eframe::egui;
use tracing_subscriber::EnvFilter;
use tracing_subscriber::fmt::format::FmtSpan;
use tracing_subscriber::registry::Registry;

use clap::Parser;
use std::fmt::Write as FMTWrite;
use tracing_subscriber::layer::SubscriberExt;
use tracing_subscriber::util::SubscriberInitExt;

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

fn setup_logger(verbose: u8, quiet: u8) {
    let log_level = get_log_level(verbose, quiet);
    let env_filter = match EnvFilter::builder().parse(&log_level) {
        Ok(filter) => filter,
        Err(_) => {
            let mut warning_msg = String::new();
            let _ = writeln!(
                &mut warning_msg,
                "Warning: Invalid log level: {}. Falling back to 'info'.",
                log_level
            );
            eprintln!("{}", warning_msg);
            EnvFilter::new("info")
        }
    };

    // Initialize Subscriber
    let subscriber = Registry::default()
        .with(env_filter)
        .with(tracing_subscriber::fmt::layer().with_span_events(FmtSpan::CLOSE));

    subscriber.init();
}

fn main() -> eframe::Result {
    let args = cli::Cli::parse();

    setup_logger(args.verbose, args.quiet);
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1400.0, 800.0])
            .with_min_inner_size([800.0, 600.0])
            .with_app_id("timsquery_viewer"),
        ..Default::default()
    };

    eframe::run_native(
        "TimsQuery Viewer",
        options,
        Box::new(|cc| Ok(Box::new(app::ViewerApp::new(cc, &args)))),
    )
}
