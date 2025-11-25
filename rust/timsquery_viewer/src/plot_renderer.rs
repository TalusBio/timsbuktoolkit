use eframe::egui;
use egui_plot::{
    Legend,
    Line,
    Plot,
    PlotPoints,
    PlotPoint,
};

use tracing::instrument;
use crate::chromatogram_processor::ChromatogramOutput;

struct LineData {
    points: Vec<PlotPoint>,
    name: String,
    stroke: egui::Stroke,
}

impl LineData {
    fn to_plot_line<'a>(&'a self) -> Line<'a> {
        Line::new(&self.name, self.points.as_slice())
            .stroke(self.stroke)
    }
}

pub struct ChromatogramLines {
    precursor_lines: Vec<LineData>,
    fragment_lines: Vec<LineData>,
    reference_id: u64,
    reference_ook0: f64,
    reference_rt_seconds: f64,
    intensity_range: (f64, f64),
    rt_seconds_range: (f64, f64),
}

impl ChromatogramLines {
    #[instrument(skip(chromatogram))]
    pub(crate) fn from_chromatogram(chromatogram: &ChromatogramOutput) -> Self {
        let precursor_lines = chromatogram
            .precursor_mzs
            .iter()
            .zip(chromatogram.precursor_intensities.iter())
            .enumerate()
            .map(|(i, (mz, intensities))| {
                let points: Vec<PlotPoint> = chromatogram
                    .retention_time_results_seconds
                    .iter()
                    .zip(intensities.iter())
                    .map(|(&rt, &intensity)| PlotPoint::new(rt as f64, intensity as f64))
                    .collect();

                let color = get_precursor_color(i);
                LineData {
                    points,
                    name: format!("Precursor m/z {:.4}", mz),
                    stroke: egui::Stroke::new(2.0, color),
                }
            })
            .collect();

        let fragment_lines = chromatogram
            .fragment_mzs
            .iter()
            .zip(chromatogram.fragment_intensities.iter())
            .enumerate()
            .map(|(i, (mz, intensities))| {
                let points: Vec<PlotPoint> = chromatogram
                    .retention_time_results_seconds
                    .iter()
                    .zip(intensities.iter())
                    .map(|(&rt, &intensity)| PlotPoint::new(rt as f64, intensity as f64))
                    .collect();

                let color = get_fragment_color(i);
                LineData {
                    points,
                    name: format!("Fragment m/z {:.4}", mz),
                    stroke: egui::Stroke::new(1.5, color),
                }
            })
            .collect();

        let rt_seconds_range = (
            chromatogram
                .retention_time_results_seconds
                .iter()
                .cloned()
                .fold(f32::INFINITY, f32::min) as f64,
            chromatogram
                .retention_time_results_seconds
                .iter()
                .cloned()
                .fold(f32::NEG_INFINITY, f32::max) as f64,
        );

        let intensity_range = {
            let mut min_intensity = f64::INFINITY;
            let mut max_intensity = f64::NEG_INFINITY;

            for intensities in chromatogram
                .precursor_intensities
                .iter()
                .chain(chromatogram.fragment_intensities.iter())
            {
                for &intensity in intensities {
                    if intensity > 0.0 {
                        min_intensity = min_intensity.min(intensity as f64);
                        max_intensity = max_intensity.max(intensity as f64);
                    }
                }
            }

            (min_intensity, max_intensity)
        };

        Self {
            precursor_lines,
            fragment_lines,
            reference_id: chromatogram.id,
            reference_ook0: chromatogram.mobility_ook0 as f64,
            reference_rt_seconds: chromatogram.rt_seconds as f64,
            intensity_range,
            rt_seconds_range,
        }
    }
}

/// Renders a chromatogram plot using egui_plot
pub fn render_chromatogram_plot(
    ui: &mut egui::Ui,
    chromatogram: &ChromatogramLines,
) {
    ui.label(format!("Elution Group ID: {}", chromatogram.reference_id));
    ui.label(format!(
        "RT: {:.2} s, Mobility: {:.4}",
        chromatogram.reference_rt_seconds, chromatogram.reference_ook0
    ));
    ui.separator();

    // Create the plot with box select zoom enabled
    let plot = Plot::new("chromatogram_plot")
        .legend(Legend::default())
        .show_axes([true, true])
        .x_axis_label("Retention Time (s)")
        .y_axis_label("Intensity");

    plot.show(ui, |plot_ui| {
        for x in chromatogram.precursor_lines.iter() {
            plot_ui.line(x.to_plot_line());
        }

        for y in chromatogram.fragment_lines.iter() {
            plot_ui.line(y.to_plot_line());
        }

    });
}

/// Get a color for precursor traces (blue-ish tones)
fn get_precursor_color(index: usize) -> egui::Color32 {
    let colors = [
        egui::Color32::from_rgb(0, 114, 178),   // Blue
        egui::Color32::from_rgb(86, 180, 233),  // Sky Blue
        egui::Color32::from_rgb(0, 158, 115),   // Bluish Green
        egui::Color32::from_rgb(0, 0, 255),     // Bright Blue
        egui::Color32::from_rgb(100, 149, 237), // Cornflower Blue
    ];
    colors[index % colors.len()]
}

/// Get a color for fragment traces (orange/red tones)
fn get_fragment_color(index: usize) -> egui::Color32 {
    let colors = [
        egui::Color32::from_rgb(230, 159, 0),   // Orange
        egui::Color32::from_rgb(213, 94, 0),    // Vermillion
        egui::Color32::from_rgb(204, 121, 167), // Reddish Purple
        egui::Color32::from_rgb(240, 228, 66),  // Yellow
        egui::Color32::from_rgb(255, 165, 0),   // Bright Orange
        egui::Color32::from_rgb(220, 50, 47),   // Red
        egui::Color32::from_rgb(255, 99, 71),   // Tomato
        egui::Color32::from_rgb(255, 140, 0),   // Dark Orange
    ];
    colors[index % colors.len()]
}
