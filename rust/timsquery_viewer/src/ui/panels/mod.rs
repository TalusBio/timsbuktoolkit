pub mod config_panel;
pub mod mobility_panel;
pub mod precursor_table_panel;
pub mod spectrum_display_panel;

pub use config_panel::{ConfigPanel, ScreenshotAction, ScreenshotState};
pub use mobility_panel::MobilityPanel;
pub use precursor_table_panel::TablePanel;
pub use spectrum_display_panel::SpectrumPanel;

use egui::Color32;

/// Shared color mapping for ion labels.
pub(crate) fn ion_color(label: &str) -> Color32 {
    match label.chars().next() {
        Some('b') | Some('B') => Color32::from_rgb(100, 149, 237), // Blue (Cornflower)
        Some('y') | Some('Y') => Color32::from_rgb(220, 80, 80),   // Red
        Some('P') => Color32::from_rgb(255, 200, 50),              // Yellow (Precursor)
        _ => Color32::from_rgb(50, 205, 50),                       // Green (Lime)
    }
}
