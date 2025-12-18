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
                    // Check if we have library intensities for mirror mode
                    let has_library = spec.library_fragment_intensities.is_some();

                    // Max-normalize observed intensities when in mirror mode
                    let observed_max = spec.intensities.iter().cloned().fold(0.0f32, f32::max);
                    let normalized_observed: Vec<f64> = if has_library && observed_max > 0.0 {
                        spec.intensities
                            .iter()
                            .map(|&v| (v / observed_max) as f64)
                            .collect()
                    } else {
                        spec.intensities.iter().map(|&v| v as f64).collect()
                    };

                    // Calculate label offset (3% of 1.0 for normalized, or 3% of max for raw)
                    let obs_max_for_offset = normalized_observed.iter().cloned().fold(0.0f64, f64::max);
                    let label_offset = obs_max_for_offset * 0.03;

                    // Draw library mirror (if present) as negative intensities
                    // Normalize only if max library intensity > 1, otherwise use raw values
                    if let Some(lib) = &spec.library_fragment_intensities {
                        let lib_max = lib.iter().cloned().fold(0.0f32, f32::max);
                        let should_normalize_lib = lib_max > 1.0;
                        let lib_display_max = if should_normalize_lib { 1.0 } else { lib_max as f64 };
                        let lib_label_offset = lib_display_max * 0.03;

                        for (idx, (&mz, &lib_int)) in spec.mz_values.iter().zip(lib).enumerate() {
                            // Normalize library intensity only if max > 1
                            let display_lib = if should_normalize_lib && lib_max > 0.0 {
                                (lib_int / lib_max) as f64
                            } else {
                                lib_int as f64
                            };

                            // Use same fragment-type color as observed spectrum
                            let label_str = &spec.fragment_labels[idx];
                            let color = Self::get_fragment_color(label_str);
                            let points = PlotPoints::new(vec![[mz, 0.0], [mz, -display_lib]]);
                            let line = Line::new(format!("lib_{}", idx), points)
                                .stroke(egui::Stroke::new(2.0, color));
                            plot_ui.line(line);

                            // Label below mirrored peak for all fragments with intensity > 0
                            if display_lib > 0.0 {
                                let label = Text::new(
                                    label_str,
                                    PlotPoint::new(mz, -display_lib - lib_label_offset),
                                    label_str,
                                )
                                .color(color);
                                plot_ui.text(label);
                            }
                        }
                    }

                    // Draw observed (normalized in mirror mode) peaks on top
                    for (idx, (&mz, &intensity)) in
                        spec.mz_values.iter().zip(&normalized_observed).enumerate()
                    {
                        let label_str = &spec.fragment_labels[idx];
                        let color = Self::get_fragment_color(label_str);

                        let points = PlotPoints::new(vec![[mz, 0.0], [mz, intensity]]);
                        let line = Line::new(label_str, points)
                            .stroke(egui::Stroke::new(2.0, color));
                        plot_ui.line(line);

                        // Label above peak for all fragments with intensity > 0
                        if intensity > 0.0 {
                            let label = Text::new(
                                label_str,
                                PlotPoint::new(mz, intensity + label_offset),
                                label_str,
                            )
                            .color(color);
                            plot_ui.text(label);
                        }
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
