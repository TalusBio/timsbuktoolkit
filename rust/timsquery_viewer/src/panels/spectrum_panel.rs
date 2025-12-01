use eframe::egui;
use egui_plot::{Bar, BarChart, Plot};

use crate::plot_renderer::MS2Spectrum;

/// Panel for displaying MS2 spectrum
pub struct SpectrumPanel;

impl SpectrumPanel {
    pub fn new() -> Self {
        Self
    }

    pub fn render(&self, ui: &mut egui::Ui, spectrum: Option<&MS2Spectrum>) {
        ui.heading("MS2 Spectrum");

        if let Some(spec) = spectrum {
            ui.label(format!("RT: {:.2} seconds", spec.rt_seconds));
            ui.separator();

            Plot::new("ms2_spectrum")
                .height(ui.available_height())
                .show_axes([true, true])
                .allow_zoom(true)
                .allow_drag(true)
                .show(ui, |plot_ui| {
                    let bars: Vec<Bar> = spec
                        .mz_values
                        .iter()
                        .zip(&spec.intensities)
                        .enumerate()
                        .map(|(idx, (&mz, &intensity))| {
                            Bar::new(mz, intensity as f64)
                                .width(0.5)
                                .name(&spec.fragment_labels[idx])
                        })
                        .collect();

                    plot_ui.bar_chart(BarChart::new("MS2 Spectrum", bars));
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
