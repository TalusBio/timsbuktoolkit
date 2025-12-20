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

use crate::plot_renderer::MS2Spectrum;

/// Panel for displaying MS2 spectrum
pub struct SpectrumPanel;

impl SpectrumPanel {
    pub fn new() -> Self {
        Self
    }

    pub fn title(&self) -> &str {
        "MS2"
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

    pub fn render(
        &mut self,
        ui: &mut egui::Ui,
        ms2_spectrum: &Option<MS2Spectrum>,
        expected_intensities: &Option<std::collections::HashMap<String, f32>>,
    ) {
        if let Some(spec) = ms2_spectrum {
            ui.label(format!("RT: {:.2} seconds", spec.rt_seconds));
            ui.separator();

            // Calculate label offset based on max intensity
            let max_intensity = spec.intensities.iter().cloned().fold(0.0f32, f32::max);
            let norm_factor = max_intensity.max(1.0);
            let label_offset = 0.03f64; // 3% of max intensity

            let expected_intensities = expected_intensities.as_ref();
            let expected_norm_factor = expected_intensities
                .map(|exp| exp.values().cloned().fold(0.0f32, f32::max).max(1.0));

            Plot::new("ms2_spectrum")
                .height(ui.available_height())
                .show_axes([true, true])
                .allow_zoom(true)
                .allow_drag(true)
                .x_axis_label("m/z")
                .y_axis_label("Intensity")
                .include_y(0.0)
                .show(ui, |plot_ui| {
                    // Draw each peak as a vertical line from 0 to intensity
                    for (idx, (&mz, &intensity)) in
                        spec.mz_values.iter().zip(&spec.intensities).enumerate()
                    {
                        let label_str = &spec.fragment_labels[idx];
                        let color = Self::get_fragment_color(label_str);
                        let y_value = (intensity / norm_factor) as f64;

                        let points = PlotPoints::new(vec![[mz, 0.0], [mz, y_value]]);
                        let line = Line::new(label_str, points).color(color);
                        plot_ui.line(line);

                        // Add label above the peak with offset
                        let label = Text::new(
                            label_str,
                            PlotPoint::new(mz, y_value + label_offset),
                            label_str,
                        )
                        .color(color);
                        plot_ui.text(label);

                        // If expected intensities are provided, draw them as dashed lines
                        // in the negative direction
                        if let Some(expected) = expected_intensities
                            && let Some(&expected_intensity) = expected.get(label_str) {
                                let y_ref_value =
                                    (expected_intensity / expected_norm_factor.expect("Expected norm factor calculated if we have expected intensities"))
                                        as f64;
                                let expected_points = PlotPoints::new(vec![
                                    [mz, 0.0],
                                    [mz, -y_ref_value],
                                ]);
                                let expected_line =
                                    Line::new(format!("expected_{}", label_str), expected_points)
                                        .color(color);
                                plot_ui.line(expected_line);
                            }
                    }
                });
        } else {
            ui.centered_and_justified(|ui| {
                ui.label("Click on XIC plot to view MS2 spectrum");
            });
        }
    }
}

impl Default for SpectrumPanel {
    fn default() -> Self {
        Self::new()
    }
}
