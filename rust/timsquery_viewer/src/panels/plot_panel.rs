use eframe::egui;

use crate::app::ComputedState;
use crate::plot_renderer::{self, PlotMode};

/// Panel for displaying chromatogram plots
pub struct PlotPanel;

impl PlotPanel {
    pub fn new() -> Self {
        Self
    }

    /// Render the plot panel
    /// Returns the clicked RT (in seconds) if the plot was clicked
    pub fn render(
        &self,
        ui: &mut egui::Ui,
        computed: &ComputedState,
        selected_index: Option<usize>,
        show_split: bool,
    ) -> Option<f64> {
        self.render_header(ui);
        ui.separator();
        self.render_content(ui, computed, selected_index, show_split)
    }

    fn render_header(&self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.heading("Chromatogram");
        });
    }

    fn render_content(
        &self,
        ui: &mut egui::Ui,
        computed: &ComputedState,
        selected_index: Option<usize>,
        show_split: bool,
    ) -> Option<f64> {
        if let Some(chromatogram) = &computed.chromatogram {
            if show_split {
                let mut clicked_rt = None;
                let link_id = "split_xic_x_axis";

                ui.label(format!("Elution Group ID: {}", chromatogram.reference_id));
                ui.label(format!(
                    "RT: {:.2} s, Mobility: {:.4}",
                    chromatogram.reference_rt_seconds, chromatogram.reference_ook0
                ));
                ui.separator();

                let remaining_height = ui.available_height();
                let separator_space = ui.spacing().item_spacing.y * 3.0;
                let plot_height = (remaining_height - separator_space) * 0.5;

                ui.allocate_ui_with_layout(
                    egui::vec2(ui.available_width(), plot_height),
                    egui::Layout::top_down(egui::Align::LEFT),
                    |ui| {
                        ui.heading("Precursor Traces");
                        if let Some(rt) = plot_renderer::render_chromatogram_plot(
                            ui,
                            chromatogram,
                            PlotMode::PrecursorsOnly,
                            Some(link_id),
                            false,
                        ) {
                            clicked_rt = Some(rt);
                        }
                    },
                );

                ui.separator();

                ui.heading("Fragment Traces");
                if let Some(rt) = plot_renderer::render_chromatogram_plot(
                    ui,
                    chromatogram,
                    PlotMode::FragmentsOnly,
                    Some(link_id),
                    false,
                ) {
                    clicked_rt = Some(rt);
                }

                clicked_rt
            } else {
                plot_renderer::render_chromatogram_plot(
                    ui,
                    chromatogram,
                    PlotMode::All,
                    None,
                    true,
                )
            }
        } else {
            if selected_index.is_some() {
                ui.label("Generating chromatogram...");
            } else {
                ui.label("Select a precursor from the table to view chromatogram");
            }
            None
        }
    }
}

impl Default for PlotPanel {
    fn default() -> Self {
        Self::new()
    }
}
