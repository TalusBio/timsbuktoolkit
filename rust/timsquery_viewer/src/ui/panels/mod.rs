pub mod config_panel;
pub mod mobility_panel;
pub mod precursor_table_panel;
pub mod spectrum_display_panel;

pub use config_panel::{
    ConfigPanel,
    ScreenshotAction,
    ScreenshotState,
};
pub use mobility_panel::MobilityPanel;
pub use precursor_table_panel::TablePanel;
pub use spectrum_display_panel::SpectrumPanel;

use egui::Color32;

/// Shared color mapping for ion labels.
pub(crate) fn ion_color(label: &str, alpha: Option<u8>) -> Color32 {
    let alpha = alpha.unwrap_or(255);
    match label.chars().next() {
        Some('b') | Some('B') => Color32::from_rgba_unmultiplied(100, 149, 237, alpha), /* Blue (Cornflower) */
        Some('y') | Some('Y') => Color32::from_rgba_unmultiplied(220, 80, 80, alpha),   // Red
        Some('P') | Some('p') => Color32::from_rgba_unmultiplied(50, 205, 50, alpha), /* Green (Lime) */
        _ => Color32::from_rgba_unmultiplied(255, 200, 50, alpha), // Yellow (Gold)
    }
}
