use eframe::egui;
use egui_plot::{
    Legend,
    Line,
    Plot,
    PlotPoints,
};

use crate::chromatogram_processor::ChromatogramOutput;

/// Renders a chromatogram plot using egui_plot
pub fn render_chromatogram_plot(
    ui: &mut egui::Ui,
    chromatogram: &ChromatogramOutput,
    reset_bounds: &mut bool,
    reset_x: &mut bool,
    reset_y: &mut bool,
) {
    ui.label(format!("Elution Group ID: {}", chromatogram.id));
    ui.label(format!(
        "RT: {:.2} s, Mobility: {:.4}",
        chromatogram.rt_seconds, chromatogram.mobility_ook0
    ));
    ui.separator();

    // Calculate data bounds for reset functionality
    let mut min_rt = f64::MAX;
    let mut max_rt = f64::MIN;
    let mut min_intensity = f64::MAX;
    let mut max_intensity = f64::MIN;

    for &rt in &chromatogram.retention_time_results_seconds {
        min_rt = min_rt.min(rt as f64);
        max_rt = max_rt.max(rt as f64);
    }

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

    // Add some padding
    let rt_padding = (max_rt - min_rt) * 0.05;
    let intensity_padding = (max_intensity - min_intensity) * 0.05;

    // Create the plot with box select zoom enabled
    let mut plot = Plot::new("chromatogram_plot")
        .legend(Legend::default())
        .allow_boxed_zoom(true)
        .allow_drag(egui::Vec2b::new(true, true))
        .show_axes([true, true])
        .x_axis_label("Retention Time (s)")
        .y_axis_label("Intensity");

    // Apply reset if requested
    if *reset_bounds {
        plot = plot.include_x(min_rt - rt_padding);
        plot = plot.include_x(max_rt + rt_padding);
        plot = plot.include_y(min_intensity - intensity_padding);
        plot = plot.include_y(max_intensity + intensity_padding);
        plot = plot.auto_bounds(egui::Vec2b::new(true, true));
        *reset_bounds = false;
        *reset_x = false;
        *reset_y = false;
    } else if *reset_x {
        plot = plot.include_x(min_rt - rt_padding);
        plot = plot.include_x(max_rt + rt_padding);
        plot = plot.auto_bounds(egui::Vec2b::new(true, false));
        *reset_x = false;
    } else if *reset_y {
        plot = plot.include_y(min_intensity - intensity_padding);
        plot = plot.include_y(max_intensity + intensity_padding);
        plot = plot.auto_bounds(egui::Vec2b::new(false, true));
        *reset_y = false;
    }

    plot.show(ui, |plot_ui| {
        // Plot precursor intensities
        for (i, (mz, intensities)) in chromatogram
            .precursor_mzs
            .iter()
            .zip(chromatogram.precursor_intensities.iter())
            .enumerate()
        {
            let points: PlotPoints = chromatogram
                .retention_time_results_seconds
                .iter()
                .zip(intensities.iter())
                .map(|(&rt, &intensity)| [rt as f64, intensity as f64])
                .collect();

            let color = get_precursor_color(i);
            let line = Line::new(points)
                .name(format!("Precursor m/z {:.4}", mz))
                .stroke(egui::Stroke::new(2.0, color));

            plot_ui.line(line);
        }

        // Plot fragment intensities
        for (i, (mz, intensities)) in chromatogram
            .fragment_mzs
            .iter()
            .zip(chromatogram.fragment_intensities.iter())
            .enumerate()
        {
            let points: PlotPoints = chromatogram
                .retention_time_results_seconds
                .iter()
                .zip(intensities.iter())
                .map(|(&rt, &intensity)| [rt as f64, intensity as f64])
                .collect();

            let color = get_fragment_color(i);
            let line = Line::new(points)
                .name(format!("Fragment m/z {:.4}", mz))
                .stroke(egui::Stroke::new(1.5, color));

            plot_ui.line(line);
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
