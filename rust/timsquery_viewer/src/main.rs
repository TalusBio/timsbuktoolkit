mod app;
mod chromatogram_processor;
mod error;
mod file_loader;
mod plot_renderer;
mod ui;

use eframe::egui;
use tracing::subscriber::set_global_default;
use tracing_subscriber::EnvFilter;
use tracing_subscriber::fmt::format::FmtSpan;
use tracing_subscriber::prelude::*;
use tracing_subscriber::registry::Registry;

#[cfg(target_os = "windows")]
use mimalloc::MiMalloc;

#[cfg(target_os = "windows")]
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

fn main() -> eframe::Result {
    let env_filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info"));
    let subscriber = Registry::default()
        .with(env_filter)
        .with(tracing_subscriber::fmt::layer().with_span_events(FmtSpan::CLOSE));

    set_global_default(subscriber).expect("Setting default subscriber failed");

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1400.0, 800.0])
            .with_min_inner_size([800.0, 600.0]),
        ..Default::default()
    };

    eframe::run_native(
        "TimsQuery Viewer",
        options,
        Box::new(|cc| Ok(Box::new(app::ViewerApp::new(cc)))),
    )
}
