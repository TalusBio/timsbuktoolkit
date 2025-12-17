use eframe::egui::{
    self,
    Color32,
};
use egui_plot::{
    Line,
    Plot,
    PlotPoint,
    PlotPoints,
    Text,
};

use crate::ui::{
    Panel,
    PanelContext,
};

/// Panel for displaying MS2 spectrum
pub struct SpectrumPanel;

impl SpectrumPanel {
    pub fn new() -> Self {
        Self
    }

    /// Get color based on fragment label prefix
    fn get_fragment_color(label: &str) -> Color32 {
        match label.chars().next() {
            Some('b') | Some('B') => Color32::from_rgb(100, 149, 237), // Blue (Cornflower)
            Some('y') | Some('Y') => Color32::from_rgb(220, 20, 60),   // Red (Crimson)
            Some('p') | Some('P') => Color32::from_rgb(255, 200, 0),   // Yellow
            _ => Color32::from_rgb(50, 205, 50),                       // Green (Lime)
        }
    }
}

impl Panel for SpectrumPanel {
    fn render(&mut self, ui: &mut egui::Ui, ctx: &mut PanelContext) {
        if let Some(spec) = &ctx.computed.ms2_spectrum {
            ui.label(format!("RT: {:.2} seconds", spec.rt_seconds));
            ui.separator();

            Plot::new("ms2_spectrum")
                .height(ui.available_height())
                .show_axes([true, true])
                .allow_zoom(true)
                .allow_drag(true)
                .x_axis_label("m/z")
                .y_axis_label("Intensity")
                .include_y(0.0)
                .show(ui, |plot_ui| {
                    // Calculate label offset based on max intensity
                    let max_intensity = spec.intensities.iter().cloned().fold(0.0f32, f32::max);
                    let label_offset = (max_intensity * 0.03) as f64; // 3% of max intensity

                    // Draw each peak as a vertical line from 0 to intensity
                    for (idx, (&mz, &intensity)) in
                        spec.mz_values.iter().zip(&spec.intensities).enumerate()
                    {
                        let label_str = &spec.fragment_labels[idx];
                        let color = Self::get_fragment_color(label_str);

                        let points = PlotPoints::new(vec![[mz, 0.0], [mz, intensity as f64]]);
                        let line = Line::new(label_str, points).color(color);
                        plot_ui.line(line);

                        // Add label above the peak with offset
                        let label = Text::new(
                            label_str,
                            PlotPoint::new(mz, intensity as f64 + label_offset),
                            label_str,
                        )
                        .color(color);
                        plot_ui.text(label);
                    }
                });
        } else {
            ui.centered_and_justified(|ui| {
                ui.label("Click on XIC plot to view MS2 spectrum");
            });
        }
    }

    fn title(&self) -> &str {
        "MS2"
    }
}

impl Default for SpectrumPanel {
    fn default() -> Self {
        Self::new()
    }
}
