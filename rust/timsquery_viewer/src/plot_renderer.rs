use eframe::egui;
use egui_plot::{
    Legend,
    Line,
    Plot,
    PlotPoint,
    PlotPoints,
    Polygon,
};

use crate::chromatogram_processor::ChromatogramOutput;
use tracing::{
    info,
    instrument,
};

const REFERENCE_RT_BAND_WIDTH_SECONDS: f64 = 10.0;

/// Specifies which traces to render in the chromatogram plot
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlotMode {
    /// Show all traces (precursors + fragments)
    All,
    /// Show only precursor traces
    PrecursorsOnly,
    /// Show only fragment traces
    FragmentsOnly,
}

#[derive(Debug)]
pub struct ChromatogramLines {
    precursor_lines: Vec<ChromatogramLine>,
    fragment_lines: Vec<ChromatogramLine>,
    pub reference_id: u64,
    pub reference_ook0: f64,
    pub reference_rt_seconds: f64,
    intensity_max: f64,
    pub rt_seconds_range: (f64, f64),
}

impl ChromatogramLines {
    #[instrument(skip(chromatogram))]
    pub(crate) fn from_chromatogram(chromatogram: &ChromatogramOutput) -> Self {
        let mut global_max_intensity = f32::NEG_INFINITY;

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

                let intensity_max = intensities
                    .iter()
                    .cloned()
                    .fold(f32::NEG_INFINITY, f32::max) as f64;
                global_max_intensity = global_max_intensity.max(intensity_max as f32);

                let color = get_precursor_color(i);
                ChromatogramLine {
                    data: LineData {
                        points,
                        name: format!("Precursor m/z {:.4}", mz),
                        stroke: egui::Stroke::new(2.0, color),
                    },
                    intensity_max,
                }
            })
            .collect();

        let fragment_lines = chromatogram
            .fragment_mzs
            .iter()
            .zip(chromatogram.fragment_intensities.iter())
            .zip(chromatogram.fragment_labels.iter())
            .enumerate()
            .map(|(i, ((mz, intensities), label))| {
                let points: Vec<PlotPoint> = chromatogram
                    .retention_time_results_seconds
                    .iter()
                    .zip(intensities.iter())
                    .map(|(&rt, &intensity)| PlotPoint::new(rt as f64, intensity as f64))
                    .collect();

                let intensity_max = intensities
                    .iter()
                    .cloned()
                    .fold(f32::NEG_INFINITY, f32::max) as f64;
                global_max_intensity = global_max_intensity.max(intensity_max as f32);

                let color = get_fragment_color(i);
                ChromatogramLine {
                    data: LineData {
                        points,
                        name: format!("{} mz={:.4}", label, mz),
                        stroke: egui::Stroke::new(1.5, color),
                    },
                    intensity_max,
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

        Self {
            precursor_lines,
            fragment_lines,
            reference_id: chromatogram.id,
            reference_ook0: chromatogram.mobility_ook0 as f64,
            reference_rt_seconds: chromatogram.rt_seconds as f64,
            intensity_max: global_max_intensity as f64,
            rt_seconds_range,
        }
    }

    fn get_fragment_intensity_max(&self) -> f64 {
        self.fragment_lines
            .iter()
            .map(|line| line.intensity_max)
            .fold(f64::NEG_INFINITY, f64::max)
    }

    fn get_precursor_intensity_max(&self) -> f64 {
        self.precursor_lines
            .iter()
            .map(|line| line.intensity_max)
            .fold(f64::NEG_INFINITY, f64::max)
    }
}

/// MS2 spectrum data at a specific retention time
#[derive(Debug, Clone)]
pub struct MS2Spectrum {
    pub mz_values: Vec<f64>,
    pub intensities: Vec<f32>,
    pub rt_seconds: f64,
    pub fragment_labels: Vec<String>,
}

/// Renders a chromatogram plot using egui_plot with custom zoom/pan controls
/// Returns the clicked RT (in seconds) if the plot was clicked
///
/// If `link_group_id` is provided, the X-axis will be linked to other plots with the same ID
/// If `show_header` is false, the elution group ID and reference RT/mobility labels are not shown
/// If `reset_bounds_applied` is false, the plot bounds will be reset to show the full data range, and the flag will be set to true
pub fn render_chromatogram_plot(
    ui: &mut egui::Ui,
    chromatogram: &ChromatogramLines,
    mode: PlotMode,
    link_group_id: Option<&str>,
    show_header: bool,
    reset_bounds_applied: &mut bool,
) -> Option<f64> {
    let mut clicked_rt = None;

    // Optionally show header information
    if show_header {
        ui.label(format!("Elution Group ID: {}", chromatogram.reference_id));
        ui.label(format!(
            "RT: {:.2} s, Mobility: {:.4}",
            chromatogram.reference_rt_seconds, chromatogram.reference_ook0
        ));
    }

    let (scroll_delta, _shift_pressed) = ui.input(|i| (i.smooth_scroll_delta, i.modifiers.shift));

    let plot_id = match mode {
        PlotMode::All => "chromatogram_plot",
        PlotMode::PrecursorsOnly => "chromatogram_plot_precursors",
        PlotMode::FragmentsOnly => "chromatogram_plot_fragments",
    };
    let mut plot = Plot::new(plot_id)
        .legend(Legend::default())
        .show_axes([true, true])
        .x_axis_label("Retention Time (s)")
        .y_axis_label("Intensity")
        .allow_zoom(false)
        .allow_drag(false)
        .allow_scroll(false);

    if let Some(link_id) = link_group_id {
        plot = plot.link_axis(link_id.to_string(), [true, false]);
    }

    plot.show(ui, |plot_ui| {
        let rt_band_half_width = REFERENCE_RT_BAND_WIDTH_SECONDS / 2.0;
        let max_polygon_height = match mode {
            PlotMode::All => chromatogram.intensity_max,
            PlotMode::PrecursorsOnly => chromatogram.get_precursor_intensity_max(),
            PlotMode::FragmentsOnly => chromatogram.get_fragment_intensity_max(),
        };

        let reference_band = Polygon::new(
            "Reference RT",
            PlotPoints::new(vec![
                [chromatogram.reference_rt_seconds - rt_band_half_width, 0.0],
                [chromatogram.reference_rt_seconds + rt_band_half_width, 0.0],
                [
                    chromatogram.reference_rt_seconds + rt_band_half_width,
                    max_polygon_height,
                ],
                [
                    chromatogram.reference_rt_seconds - rt_band_half_width,
                    max_polygon_height,
                ],
            ]),
        )
        .fill_color(egui::Color32::from_rgba_premultiplied(128, 128, 128, 26))
        .stroke(egui::Stroke::NONE);

        plot_ui.polygon(reference_band);

        match mode {
            PlotMode::All | PlotMode::PrecursorsOnly => {
                for line in chromatogram.precursor_lines.iter() {
                    plot_ui.line(line.data.to_plot_line());
                }
            }
            _ => {}
        }

        match mode {
            PlotMode::All | PlotMode::FragmentsOnly => {
                for line in chromatogram.fragment_lines.iter() {
                    plot_ui.line(line.data.to_plot_line());
                }
            }
            _ => {}
        }

        let plot_hovered = plot_ui.response().hovered();
        if plot_hovered && scroll_delta.length_sq() > 0.0 {
            let zoom_speed = 0.05;
            let scroll_y = scroll_delta.y;
            let scroll_x = scroll_delta.x;

            let zoom_amount_y = (scroll_y * zoom_speed / 10.0).exp();
            let zoom_amount_x = (scroll_x * zoom_speed / 10.0).exp();

            let zoom_factor = egui::Vec2::new(zoom_amount_x, zoom_amount_y);

            plot_ui.zoom_bounds_around_hovered(zoom_factor);
        }

        let pointer_drag_delta = plot_ui.pointer_coordinate_drag_delta();
        if pointer_drag_delta.x != 0.0 || pointer_drag_delta.y != 0.0 {
            let pan_delta = egui::Vec2::new(-pointer_drag_delta.x, -pointer_drag_delta.y);
            plot_ui.translate_bounds(pan_delta);
        }

        if !*reset_bounds_applied {
            plot_ui.set_plot_bounds_x(
                chromatogram.rt_seconds_range.0..=chromatogram.rt_seconds_range.1,
            );
            plot_ui.set_plot_bounds_y(0.0..=max_polygon_height);
            *reset_bounds_applied = true;
        } else {
            let bounds = plot_ui.plot_bounds();

            let y_min = bounds.min()[1];
            let y_max = bounds.max()[1];
            let clamped_y_min = 0.0;
            let clamped_y_max = y_max.min(max_polygon_height);

            if y_min != clamped_y_min || y_max != clamped_y_max {
                plot_ui.set_plot_bounds_y(clamped_y_min..=clamped_y_max);
            }

            let x_min = bounds.min()[0];
            let x_max = bounds.max()[0];
            let clamped_x_min = x_min.max(chromatogram.rt_seconds_range.0);
            let clamped_x_max = x_max.min(chromatogram.rt_seconds_range.1);

            if x_min != clamped_x_min || x_max != clamped_x_max {
                plot_ui.set_plot_bounds_x(clamped_x_min..=clamped_x_max);
            }
        }

        if plot_ui.response().clicked()
            && let Some(pointer_pos) = plot_ui.pointer_coordinate()
        {
            clicked_rt = Some(pointer_pos.x);
            info!("Plot clicked at RT: {:.2}s", pointer_pos.x);
        }
    });

    clicked_rt
}

#[derive(Debug)]
struct LineData {
    points: Vec<PlotPoint>,
    name: String,
    stroke: egui::Stroke,
}

impl LineData {
    fn to_plot_line<'a>(&'a self) -> Line<'a> {
        Line::new(&self.name, self.points.as_slice()).stroke(self.stroke)
    }
}

#[derive(Debug)]
pub struct ChromatogramLine {
    data: LineData,
    intensity_max: f64,
}

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
