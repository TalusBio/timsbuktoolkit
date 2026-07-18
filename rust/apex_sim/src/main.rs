//! apex_sim -- interactive simulator for the timsseek apex finder.
//!
//! Synthesizes extraction products (real co-eluting fragments, contaminants,
//! noise) with live-tunable knobs, feeds them to the REAL `timsseek` apex
//! finder, and displays the raw chromatograms, a heatmap, and every per-cycle
//! intermediate trace with both scoring passes' apexes marked.

use apex_sim::app;
use eframe::egui;

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1400.0, 900.0])
            .with_min_inner_size([900.0, 600.0])
            .with_app_id("apex_sim"),
        ..Default::default()
    };

    eframe::run_native(
        "Apex Finder Simulator",
        options,
        Box::new(|_cc| Ok(Box::new(app::SimApp::default()))),
    )
}
