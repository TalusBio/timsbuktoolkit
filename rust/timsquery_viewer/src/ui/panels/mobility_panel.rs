use eframe::egui;
use egui::Color32;
use egui_plot::{HLine, Plot, PlotPoints, Points};

use super::ion_color;
use crate::computed_state::MobilityData;

/// Consolidated mobility visualization.
/// x = intensity-weighted mean m/z, y = intensity-weighted mean mobility.
/// Point size proportional to relative intensity, color = ion type.
/// Dashed horizontal lines = 1x (integration) mobility range.
/// Solid horizontal lines = 2x (query) mobility range.
#[derive(Default)]
pub struct MobilityPanel;

impl MobilityPanel {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn title(&self) -> &str {
        "Mobility"
    }

    pub fn render(&self, ui: &mut egui::Ui, data: Option<&MobilityData>) {
        let Some(data) = data else {
            ui.centered_and_justified(|ui| {
                ui.label("Click on XIC plot to view mobility");
            });
            return;
        };

        const MIN_RADIUS: f32 = 3.0;
        const MAX_RADIUS: f32 = 14.0;

        let max_intensity = data
            .ions
            .iter()
            .map(|ion| ion.total_intensity)
            .fold(f64::NEG_INFINITY, f64::max);

        // Include a generation counter in the plot ID so egui_plot
        // resets its persisted zoom bounds whenever the data changes.
        let plot_id = format!("mobility_summary_{}", data.generation);

        let plot = Plot::new(plot_id)
            .height(ui.available_height())
            .x_axis_label("m/z")
            .y_axis_label("1/K0 (V·s/cm²)")
            .allow_zoom(true)
            .allow_drag(true);

        plot.show(ui, |plot_ui| {
            // 2x (wide query) range — solid lines
            let (wide_lo, wide_hi) = data.wide_mobility_range;
            plot_ui.hline(HLine::new("2x range lo", wide_lo).color(Color32::from_rgb(120, 120, 120)));
            plot_ui.hline(HLine::new("2x range hi", wide_hi).color(Color32::from_rgb(120, 120, 120)));

            // 1x (integration) range — dashed lines
            let (tol_lo, tol_hi) = data.mobility_range;
            plot_ui.hline(
                HLine::new("1x range lo", tol_lo)
                    .color(Color32::from_rgb(200, 80, 80))
                    .style(egui_plot::LineStyle::dashed_dense()),
            );
            plot_ui.hline(
                HLine::new("1x range hi", tol_hi)
                    .color(Color32::from_rgb(200, 80, 80))
                    .style(egui_plot::LineStyle::dashed_dense()),
            );

            // Reference mobility — thin dashed
            plot_ui.hline(
                HLine::new("ref mobility", data.ref_mobility)
                    .color(Color32::from_rgb(255, 60, 60))
                    .style(egui_plot::LineStyle::dashed_loose()),
            );

            for ion in &data.ions {
                let frac = if max_intensity > 0.0 {
                    (ion.total_intensity / max_intensity) as f32
                } else {
                    0.0
                };
                let radius = MIN_RADIUS + frac.sqrt() * (MAX_RADIUS - MIN_RADIUS);
                let color = ion_color(&ion.label);

                plot_ui.points(
                    Points::new(&ion.label, PlotPoints::new(vec![[ion.mean_mz, ion.mean_mobility]]))
                        .color(color)
                        .radius(radius),
                );
                plot_ui.text(
                    egui_plot::Text::new(
                        format!("{}_label", ion.label),
                        egui_plot::PlotPoint::new(ion.mean_mz, ion.mean_mobility),
                        &ion.label,
                    )
                    .color(color),
                );
            }
        });
    }
}
