use eframe::egui;
use egui_plot::{
    Legend,
    Line,
    Plot,
    PlotPoints,
    PlotPoint,
    Polygon,
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
            let mut min_intensity = 0.0f64;
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

            for intensities in chromatogram
                .fragment_intensities
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
            println!("Intensity range: min {}, max {}", min_intensity, max_intensity);

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

/// Renders a chromatogram plot using egui_plot with custom zoom/pan controls
pub fn render_chromatogram_plot(
    ui: &mut egui::Ui,
    chromatogram: &ChromatogramLines,
) {
    ui.label(format!("Elution Group ID: {}", chromatogram.reference_id));
    ui.label(format!(
        "RT: {:.2} s, Mobility: {:.4}",
        chromatogram.reference_rt_seconds, chromatogram.reference_ook0
    ));
    ui.label("💡 Scroll: zoom X-axis | Shift+Scroll: zoom Y-axis | Drag: pan");
    ui.separator();

    // Get input state before entering plot closure
    let scroll_delta = ui.input(|i| i.smooth_scroll_delta);
    let shift_pressed = ui.input(|i| i.modifiers.shift);

    // Create the plot with custom controls (disable defaults)
    // Manually set bounds to match data range
    let plot = Plot::new("chromatogram_plot")
        .legend(Legend::default())
        .show_axes([true, true])
        .x_axis_label("Retention Time (s)")
        .y_axis_label("Intensity")
        .allow_zoom(false)
        .allow_drag(false)
        .allow_scroll(false)
        .include_x(chromatogram.rt_seconds_range.0)
        .include_x(chromatogram.rt_seconds_range.1)
        .include_y(chromatogram.intensity_range.0)
        .include_y(chromatogram.intensity_range.1);

    plot.show(ui, |plot_ui| {
        // Add reference RT band (10 seconds wide, subtle gray)
        let rt_band_half_width = 5.0;
        let reference_band = Polygon::new(
            "Reference RT",
            PlotPoints::new(vec![
                [chromatogram.reference_rt_seconds - rt_band_half_width, 0.0],
                [chromatogram.reference_rt_seconds + rt_band_half_width, 0.0],
                [chromatogram.reference_rt_seconds + rt_band_half_width, chromatogram.intensity_range.1],
                [chromatogram.reference_rt_seconds - rt_band_half_width, chromatogram.intensity_range.1],
            ])
        )
        .fill_color(egui::Color32::from_rgba_premultiplied(128, 128, 128, 26)) // 10% opacity gray
        .stroke(egui::Stroke::NONE);

        plot_ui.polygon(reference_band);

        // Draw precursor lines
        for line in chromatogram.precursor_lines.iter() {
            plot_ui.line(line.to_plot_line());
        }

        // Draw fragment lines
        for line in chromatogram.fragment_lines.iter() {
            plot_ui.line(line.to_plot_line());
        }

        // Custom zoom handling (using input state captured before closure)
        if scroll_delta.length_sq() > 0.0 {
            // Zoom speed factor
            let zoom_speed = 0.1;
            let scroll_y = scroll_delta.y;

            // Calculate zoom factor based on scroll
            let zoom_amount = (scroll_y * zoom_speed / 10.0).exp();

            // Apply zoom only to x-axis by default, y-axis if Shift is pressed
            let zoom_factor = if shift_pressed {
                egui::Vec2::new(1.0, zoom_amount) // Zoom y-axis only
            } else {
                egui::Vec2::new(zoom_amount, 1.0) // Zoom x-axis only
            };

            plot_ui.zoom_bounds_around_hovered(zoom_factor);
        }

        // Custom pan handling
        let pointer_drag_delta = plot_ui.pointer_coordinate_drag_delta();
        if pointer_drag_delta.x != 0.0 || pointer_drag_delta.y != 0.0 {
            // Invert the drag delta for natural panning
            let pan_delta = egui::Vec2::new(-pointer_drag_delta.x, -pointer_drag_delta.y);
            plot_ui.translate_bounds(pan_delta);
        }

        // Clamp bounds to valid ranges
        let bounds = plot_ui.plot_bounds();

        // Clamp y-axis to [0, max_intensity]
        let y_min = bounds.min()[1];
        let y_max = bounds.max()[1];
        let clamped_y_min = y_min.max(0.0);
        let clamped_y_max = y_max.min(chromatogram.intensity_range.1);

        if y_min != clamped_y_min || y_max != clamped_y_max {
            plot_ui.set_plot_bounds_y(clamped_y_min..=clamped_y_max);
        }

        // Optional: clamp x-axis to data range
        let x_min = bounds.min()[0];
        let x_max = bounds.max()[0];
        let clamped_x_min = x_min.max(chromatogram.rt_seconds_range.0);
        let clamped_x_max = x_max.min(chromatogram.rt_seconds_range.1);

        if x_min != clamped_x_min || x_max != clamped_x_max {
            plot_ui.set_plot_bounds_x(clamped_x_min..=clamped_x_max);
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
