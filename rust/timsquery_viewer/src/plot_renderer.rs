use eframe::egui;
use egui_plot::{
    Legend,
    Line,
    Plot,
    PlotPoint,
    PlotPoints,
    Polygon,
};
use timscentroid::rt_mapping::{
    CycleToRTMapping,
    MS1CycleIndex,
    RTIndex,
};
use timsseek::scoring::apex_finding::{
    ApexScore,
    ScoreTraces,
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

pub enum AutoZoomMode {
    Disabled,
    PeakApex,
    Range,
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

#[derive(Debug)]
pub struct ScoreLines {
    main_score_line: LineData,
    lines: Vec<LineData>,
    apex_score: ApexScore,
    rt_seconds_range: (f64, f64),
}

impl ScoreLines {
    #[instrument(skip_all)]
    pub(crate) fn from_scores(
        apex: ApexScore,
        scores: &ScoreTraces,
        mapper: &CycleToRTMapping<MS1CycleIndex>,
        cycle_offset: usize,
    ) -> Self {
        let mut lines: Vec<_> = scores
            .iter_scores()
            .map(|(name, trace)| {
                let max_val = trace.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
                info!("Max score for {}: {}", name, max_val);
                let norm_factor = max_val.max(1e-6);
                let inv_norm_factor = (1.0 / norm_factor) as f64;
                let inv_norm_factor = if name == "main_score" {
                    info!("Main score trace length: {}", trace.len());
                    1.0
                } else {
                    inv_norm_factor
                };
                let points: Vec<PlotPoint> = trace
                    .iter()
                    .enumerate()
                    .map(|(cycle_idx, score)| {
                        let global_idx = cycle_idx + cycle_offset;
                        mapper
                            .rt_milis_for_index(&MS1CycleIndex::new(global_idx as u32))
                            .map(|rt| {
                                let rt_seconds = rt as f64 / 1000.0;
                                PlotPoint::new(rt_seconds, *score as f64 * inv_norm_factor)
                            })
                            .unwrap()
                    })
                    .collect();

                LineData {
                    points,
                    name: name.into(),
                    stroke: egui::Stroke::new(1.5, egui::Color32::LIGHT_BLUE),
                }
            })
            .collect();

        // Check that  all lines are the same length
        let first_line_len = lines.first().map(|line| line.points.len()).unwrap_or(0);
        for line in &lines {
            assert_eq!(
                line.points.len(),
                first_line_len,
                "All score lines should have the same number of points"
            );
        }

        let main_score_line = lines
            .pop_if(|x| x.name == "main_score")
            .expect("There should be a main_score line");

        let rt_seconds_range = (
            lines
                .first()
                .and_then(|line| line.points.first())
                .map(|pt| pt.x)
                .unwrap_or(0.0),
            lines
                .last()
                .and_then(|line| line.points.last())
                .map(|pt| pt.x)
                .unwrap_or(0.0),
        );
        Self {
            main_score_line,
            lines,
            apex_score: apex,
            rt_seconds_range,
        }
    }

    pub fn render(
        &self,
        ui: &mut egui::Ui,
        link_group_id: Option<&str>,
        auto_zoom_frame_counter: &mut u8,
        auto_zoom_mode: &AutoZoomMode,
        // Lines to add vertical markers for
        label_lines: &[(String, f64)],
    ) -> Option<f64> {
        let (plot_id_top, plot_id_bot) = match link_group_id {
            Some(id) => (
                format!("score_traces_top_{}", id),
                format!("score_traces_bot_{}", id),
            ),
            None => (
                "score_traces_top".to_string(),
                "score_traces_bot".to_string(),
            ),
        };
        let (scroll_delta, _shift_pressed) =
            ui.input(|i| (i.smooth_scroll_delta, i.modifiers.shift));
        let half_height = ui.available_height() * 0.5;
        let top_plot = Plot::new(plot_id_top)
            .legend(Legend::default())
            .height(half_height)
            .show_axes([true, true])
            .allow_zoom(false)
            .allow_drag(false)
            .allow_scroll(false)
            .x_axis_label("Retention Time (s)")
            .y_axis_label("Normalized Score");

        let bot_plot = Plot::new(plot_id_bot)
            .height(half_height)
            .show_axes([true, true])
            .allow_zoom(false)
            .allow_drag(false)
            .allow_scroll(false)
            .x_axis_label("Retention Time (s)")
            .y_axis_label("Main Score");

        let (top_plot, bot_plot) = if let Some(link_id) = link_group_id {
            const ONLY_X_AXIS: [bool; 2] = [true, false];
            (
                top_plot.link_axis(link_id.to_string(), ONLY_X_AXIS),
                bot_plot.link_axis(link_id.to_string(), ONLY_X_AXIS),
            )
        } else {
            (top_plot, bot_plot)
        };

        let mut clicked_rt = None;

        top_plot.include_y(0.0).show(ui, |plot_ui| {
            for line in &self.lines {
                plot_ui.line(line.to_plot_line());
            }
            plot_reflines(label_lines, plot_ui, 0.0, 1.0);
            zoom_behavior(plot_ui, &scroll_delta);
            if *auto_zoom_frame_counter <= 0 {
                // Also ... we are letting all auto zoom modes apply on the
                // companion plot.
                // Since they are normalized the max will be 1.0
                clamp_bounds(plot_ui, 1.0, self.rt_seconds_range);
            }
            if plot_ui.response().clicked()
                && let Some(pointer_pos) = plot_ui.pointer_coordinate()
            {
                clicked_rt = Some(pointer_pos.x);
                info!("Plot clicked at RT: {:.2}s", pointer_pos.x);
            }
        });

        bot_plot.include_y(0.0).show(ui, |plot_ui| {
            let max_y = self
                .main_score_line
                .points
                .iter()
                .map(|pt| pt.y)
                .fold(0.0, f64::max);
            plot_ui.line(self.main_score_line.to_plot_line());
            plot_reflines(label_lines, plot_ui, 0.0, max_y);

            let apex_rt_seconds = self.apex_score.retention_time_ms as f64 / 1000.0;
            let apex_line = Polygon::new(
                "Apex RT",
                PlotPoints::new(vec![
                    [apex_rt_seconds - 0.1, 0.0],
                    [apex_rt_seconds + 0.1, 0.0],
                    [apex_rt_seconds + 0.1, max_y],
                    [apex_rt_seconds - 0.1, max_y],
                ]),
            )
            .fill_color(egui::Color32::from_rgba_premultiplied(255, 0, 0, 26))
            .stroke(egui::Stroke::NONE);

            plot_ui.polygon(apex_line);

            zoom_behavior(plot_ui, &scroll_delta);
            if *auto_zoom_frame_counter > 0 {
                // plot_ui.set_plot_bounds_x(self.rt_seconds_range.0..=self.rt_seconds_range.1);
                match auto_zoom_mode {
                    AutoZoomMode::Range => {
                        plot_ui
                            .set_plot_bounds_x(self.rt_seconds_range.0..=self.rt_seconds_range.1);
                        plot_ui.set_plot_bounds_y(0.0..=max_y);
                    }
                    AutoZoomMode::PeakApex => {
                        plot_ui
                            .set_plot_bounds_x((apex_rt_seconds - 30.0)..=(apex_rt_seconds + 30.0));
                        plot_ui.set_plot_bounds_y(0.0..=max_y);
                    }
                    AutoZoomMode::Disabled => {}
                }
                *auto_zoom_frame_counter -= 1;
            } else {
                clamp_bounds(plot_ui, max_y, self.rt_seconds_range);
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
/// If `auto_zoom_frame_counter` is greater than 0, the plot bounds will be reset to show the full data range, and the counter will be decremented
pub fn render_chromatogram_plot(
    ui: &mut egui::Ui,
    chromatogram: &ChromatogramLines,
    mode: PlotMode,
    link_group_id: Option<&str>,
    show_header: bool,
    auto_zoom_frame_counter: &mut u8,
    auto_zoom_mode: &AutoZoomMode,
    // Lines to add vertical markers for
    label_lines: &[(String, f64)],
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
        const ONLY_X_AXIS: [bool; 2] = [true, false];
        plot = plot.link_axis(link_id.to_string(), ONLY_X_AXIS);
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
        plot_reflines(label_lines, plot_ui, 0.0, max_polygon_height);

        zoom_behavior(plot_ui, &scroll_delta);
        if *auto_zoom_frame_counter > 0 {
            match auto_zoom_mode {
                AutoZoomMode::Range => {
                    plot_ui.set_plot_bounds_x(
                        chromatogram.rt_seconds_range.0..=chromatogram.rt_seconds_range.1,
                    );
                    plot_ui.set_plot_bounds_y(0.0..=max_polygon_height);
                }
                AutoZoomMode::PeakApex => {
                    plot_ui.set_plot_bounds_y(0.0..=max_polygon_height);
                }
                AutoZoomMode::Disabled => {}
            }
            *auto_zoom_frame_counter -= 1;
        } else {
            clamp_bounds(plot_ui, max_polygon_height, chromatogram.rt_seconds_range);
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

fn zoom_behavior(plot_ui: &mut egui_plot::PlotUi, scroll_delta: &egui::Vec2) {
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
}

fn clamp_bounds(plot_ui: &mut egui_plot::PlotUi, y_max_clamp: f64, x: (f64, f64)) {
    let bounds = plot_ui.plot_bounds();

    let y_min = bounds.min()[1];
    let y_max = bounds.max()[1];
    let clamped_y_min = 0.0;
    let clamped_y_max = y_max.min(y_max_clamp);

    if y_min != clamped_y_min || y_max != clamped_y_max {
        plot_ui.set_plot_bounds_y(clamped_y_min..=clamped_y_max);
    }

    let x_min = bounds.min()[0];
    let x_max = bounds.max()[0];
    let clamped_x_min = x_min.max(x.0);
    let clamped_x_max = x_max.min(x.1);

    if x_min != clamped_x_min || x_max != clamped_x_max {
        plot_ui.set_plot_bounds_x(clamped_x_min..=clamped_x_max);
    }
}

fn plot_reflines(
    label_lines: &[(String, f64)],
    plot_ui: &mut egui_plot::PlotUi,
    min_y: f64,
    max_y: f64,
) {
    for (label, rt_seconds) in label_lines {
        let vertical_line = Line::new(
            label.as_str(),
            PlotPoints::new(vec![[*rt_seconds, min_y], [*rt_seconds, max_y]]),
        )
        .color(egui::Color32::from_rgb(254, 0, 0))
        .stroke(egui::Stroke::new(1.0, egui::Color32::from_rgb(254, 0, 0)))
        .style(egui_plot::LineStyle::dashed_loose());
        plot_ui.line(vertical_line);
    }
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
