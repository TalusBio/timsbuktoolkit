mod app;
mod chromatogram_processor;
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
use std::fs::File;
use std::io::Write as IOWrite;
use tracing_subscriber::layer::SubscriberExt;
use tracing_subscriber::util::SubscriberInitExt;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Path to record session events to (or "-" for stdout)
    #[arg(long)]
    record_session: Option<String>,
}

#[cfg(target_os = "windows")]
use mimalloc::MiMalloc;

#[cfg(target_os = "windows")]
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

fn setup_logger() {
    let app_level = std::env::var("RUST_LOG").unwrap_or_else(|_| "info".to_string());
    let env_filter = match EnvFilter::builder()
        .parse(&app_level) {
            Ok(filter) => filter,
            Err(_) => {
                let mut warning_msg = String::new();
                let _ = writeln!(
                    &mut warning_msg,
                    "Warning: Invalid RUST_LOG value: {}. Falling back to 'info'.",
                    app_level
                );
                eprintln!("{}", warning_msg);
                EnvFilter::new("info")
            }
        };

    // 5. Initialize Subscriber
    let subscriber = Registry::default()
        .with(env_filter)
        .with(tracing_subscriber::fmt::layer().with_span_events(FmtSpan::CLOSE));

    subscriber.init(); // simpler than set_global_default + expect
}

fn main() -> eframe::Result {
    let args = Cli::parse();

    setup_logger();

    let session_log_writer: Option<Box<dyn IOWrite + Send + Sync>> =
        match args.record_session.as_deref() {
            Some("-") => Some(Box::new(std::io::stdout())),
            Some(path) => {
                let file = File::create(path).expect("Failed to create session log file");
                Some(Box::new(file))
            }
            None => None,
        };

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
        Box::new(|cc| Ok(Box::new(app::ViewerApp::new(cc, session_log_writer)))),
    )
}
