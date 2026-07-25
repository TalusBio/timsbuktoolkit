//! egui rendering helpers: chromatogram overlay, heatmap, intermediate traces.

use eframe::egui::{
    self,
    Color32,
    ColorImage,
    TextureHandle,
    TextureOptions,
};
use egui_plot::{
    Bar,
    BarChart,
    Legend,
    Line,
    LineStyle,
    Plot,
    PlotPoints,
};

use crate::scorer::ScoreResult;
use crate::sim::SimData;

/// Color a transition by its ion series (b/y/a/c/x/z) or unknown (`?`).
/// Series inferred from the leading char of the label.
pub fn fragment_type_color(label: &str) -> Color32 {
    match label.chars().next() {
        Some('b') => Color32::from_rgb(31, 119, 180),  // blue
        Some('y') => Color32::from_rgb(214, 39, 40),   // red
        Some('a') => Color32::from_rgb(23, 190, 207),  // teal
        Some('c') => Color32::from_rgb(44, 160, 44),   // green
        Some('x') => Color32::from_rgb(255, 127, 14),  // orange
        Some('z') => Color32::from_rgb(148, 103, 189), // purple
        Some('?') => Color32::from_rgb(227, 119, 194), // pink (unknown ion)
        _ => Color32::from_rgb(160, 160, 160),         // other
    }
}

/// Gold, reserved for precursor (M+n) traces.
const PRECURSOR_COLOR: Color32 = Color32::from_rgb(212, 175, 55);

/// Distinct-ish colors for score traces (bottom panel).
pub fn palette(idx: usize) -> Color32 {
    const COLORS: [Color32; 8] = [
        Color32::from_rgb(31, 119, 180),
        Color32::from_rgb(255, 127, 14),
        Color32::from_rgb(44, 160, 44),
        Color32::from_rgb(214, 39, 40),
        Color32::from_rgb(148, 103, 189),
        Color32::from_rgb(140, 86, 75),
        Color32::from_rgb(227, 119, 194),
        Color32::from_rgb(23, 190, 207),
    ];
    COLORS[idx % COLORS.len()]
}

fn line_from(name: &str, ys: &[f32]) -> Line<'static> {
    let pts: PlotPoints = ys
        .iter()
        .enumerate()
        .map(|(x, y)| [x as f64, *y as f64])
        .collect();
    Line::new(name.to_owned(), pts)
}

/// Bright green, reserved for the ground-truth apex.
const TRUE_APEX_COLOR: Color32 = Color32::from_rgb(0, 230, 90);
const PASS1_COLOR: Color32 = Color32::GRAY;
const PASS2_COLOR: Color32 = Color32::from_rgb(255, 80, 80);

/// Draw a vertical marker as a two-point line spanning [0, y_max].
///
/// An empty `name` drops the marker from the legend (egui_plot filters empty
/// names) -- used for the detected-apex markers so the legend stays clean; the
/// caption above the plot documents their colors.
fn vmarker(
    name: &str,
    x: f64,
    y_max: f64,
    color: Color32,
    style: LineStyle,
    width: f32,
) -> Line<'static> {
    let pts = PlotPoints::new(vec![[x, 0.0], [x, y_max]]);
    Line::new(name.to_owned(), pts)
        .color(color)
        .style(style)
        .width(width)
}

/// One-line legend for the vertical apex markers (colored inline).
pub fn apex_caption(ui: &mut egui::Ui) {
    ui.horizontal(|ui| {
        ui.label("apex markers:");
        ui.colored_label(TRUE_APEX_COLOR, "TRUE (realized, off-grid)");
        ui.label("·");
        ui.colored_label(PASS1_COLOR, "pass1 (quick)");
        ui.label("·");
        ui.colored_label(PASS2_COLOR, "pass2 (full)");
    });
}

/// Top panel: overlaid chromatograms. Real solid, contaminants dashed,
/// precursors thicker. Vertical markers: TRUE (green) + pass1/pass2 apexes.
pub fn chromatograms(
    ui: &mut egui::Ui,
    data: &SimData,
    score: Option<&ScoreResult>,
    true_apex: f32,
) {
    let mut y_max = 1.0_f64;
    for r in data.fragment_rows.iter().chain(data.precursor_rows.iter()) {
        for &v in &r.intensities {
            y_max = y_max.max(v as f64);
        }
    }

    Plot::new("chromatograms")
        .legend(Legend::default())
        .height(230.0)
        .include_y(0.0)
        .show(ui, |pui| {
            for r in data.fragment_rows.iter() {
                let mut line =
                    line_from(&r.label, &r.intensities).color(fragment_type_color(&r.label));
                // Absent (expected-but-unobserved) transitions drawn dashed.
                if r.is_absent {
                    line = line.style(LineStyle::dashed_loose());
                }
                pui.line(line);
            }
            for r in data.precursor_rows.iter() {
                let line = line_from(&r.label, &r.intensities)
                    .color(PRECURSOR_COLOR)
                    .width(2.5);
                pui.line(line);
            }
            // TRUE apex kept named (green) so it is easy to locate the real
            // peak; detected apexes get empty names (dropped from legend).
            pui.line(vmarker(
                "TRUE apex",
                true_apex as f64,
                y_max,
                TRUE_APEX_COLOR,
                LineStyle::Solid,
                2.5,
            ));
            if let Some(s) = score {
                pui.line(vmarker(
                    "",
                    s.pass1.apex_cycle as f64,
                    y_max,
                    PASS1_COLOR,
                    LineStyle::dashed_dense(),
                    2.0,
                ));
                pui.line(vmarker(
                    "",
                    s.pass2.joint_apex_cycle as f64,
                    y_max,
                    PASS2_COLOR,
                    LineStyle::Solid,
                    2.0,
                ));
            }
        });
}

/// Map t in [0,1] to an RGB triple (dark-blue -> cyan -> yellow -> red).
fn colormap(t: f32) -> [u8; 3] {
    let t = t.clamp(0.0, 1.0);
    let (r, g, b) = if t < 0.33 {
        let u = t / 0.33;
        (0.0, u * 0.7, 0.3 + u * 0.7)
    } else if t < 0.66 {
        let u = (t - 0.33) / 0.33;
        (u, 0.7 + u * 0.3, 1.0 - u)
    } else {
        let u = (t - 0.66) / 0.34;
        (1.0, 1.0 - u * 0.8, 0.0)
    };
    [(r * 255.0) as u8, (g * 255.0) as u8, (b * 255.0) as u8]
}

/// Build a heatmap texture: rows = transitions (fragments), cols = cycles.
pub fn build_heatmap(ctx: &egui::Context, data: &SimData) -> (TextureHandle, usize) {
    let rows = &data.fragment_rows;
    let n_rows = rows.len().max(1);
    let n_cols = rows
        .first()
        .map(|r| r.intensities.len())
        .unwrap_or(1)
        .max(1);

    let global_max = rows
        .iter()
        .flat_map(|r| r.intensities.iter())
        .cloned()
        .fold(0.0_f32, f32::max)
        .max(1e-6);

    let mut rgb = vec![0u8; n_cols * n_rows * 3];
    for (y, r) in rows.iter().enumerate() {
        for (x, &v) in r.intensities.iter().enumerate() {
            // log-scale for dynamic range
            let t = ((v / global_max).max(0.0)).powf(0.4);
            let c = colormap(t);
            let idx = (y * n_cols + x) * 3;
            rgb[idx..idx + 3].copy_from_slice(&c);
        }
    }

    let image = ColorImage::from_rgb([n_cols, n_rows], &rgb);
    let tex = ctx.load_texture("heatmap", image, TextureOptions::NEAREST);
    (tex, n_rows)
}

/// Middle panel: heatmap image + apex column overlays (TRUE + pass1/pass2).
pub fn heatmap(
    ui: &mut egui::Ui,
    tex: &TextureHandle,
    n_cols: usize,
    score: Option<&ScoreResult>,
    true_apex: f32,
) {
    let avail_w = ui.available_width();
    let height = 150.0;
    let size = egui::vec2(avail_w, height);
    let resp = ui.add(
        egui::Image::new(tex)
            .fit_to_exact_size(size)
            .texture_options(TextureOptions::NEAREST),
    );
    let rect = resp.rect;

    let painter = ui.painter_at(rect);
    let col_w = rect.width() / n_cols.max(1) as f32;
    let x_of = |cycle: f32| rect.left() + (cycle + 0.5) * col_w;
    painter.vline(
        x_of(true_apex),
        rect.y_range(),
        egui::Stroke::new(2.0, TRUE_APEX_COLOR),
    );
    if let Some(s) = score {
        painter.vline(
            x_of(s.pass1.apex_cycle as f32),
            rect.y_range(),
            egui::Stroke::new(1.5, PASS1_COLOR),
        );
        painter.vline(
            x_of(s.pass2.joint_apex_cycle as f32),
            rect.y_range(),
            egui::Stroke::new(2.0, PASS2_COLOR),
        );
    }
}

/// Which intermediate traces to display.
#[derive(Debug, Clone)]
pub struct TraceToggles {
    pub cosine: bool,
    pub scribe: bool,
    pub lazyscore: bool,
    pub log_intensity: bool,
    pub ms1_precursor: bool,
    pub apex_profile: bool,
}

impl Default for TraceToggles {
    fn default() -> Self {
        Self {
            cosine: true,
            scribe: true,
            lazyscore: false,
            log_intensity: false,
            ms1_precursor: false,
            apex_profile: true,
        }
    }
}

/// Bottom panel: the apex-finder intermediate traces, min-max normalized to
/// [0,1] so heterogeneous scales overlay legibly. Apex markers overlaid.
pub fn traces(ui: &mut egui::Ui, score: &ScoreResult, tog: &TraceToggles, true_apex: f32) {
    let t = &score.traces;
    let mut series: Vec<(&str, &[f32])> = Vec::new();
    if tog.cosine {
        series.push(("cosine", &t.cosine_trace));
    }
    if tog.scribe {
        series.push(("scribe", &t.ms2_scribe));
    }
    if tog.lazyscore {
        series.push(("lazyscore", &t.ms2_lazyscore));
    }
    if tog.log_intensity {
        series.push(("log_intensity", &t.ms2_log_intensity));
    }
    if tog.ms1_precursor {
        series.push(("ms1_precursor", &t.ms1_precursor_trace));
    }
    if tog.apex_profile {
        series.push(("apex_profile", &t.apex_profile));
    }

    Plot::new("traces")
        .legend(Legend::default())
        .height(230.0)
        .include_y(0.0)
        .include_y(1.05)
        .show(ui, |pui| {
            for (i, (name, ys)) in series.iter().enumerate() {
                pui.line(line_from(name, &normalize(ys)).color(palette(i)));
            }
            pui.line(vmarker(
                "TRUE apex",
                true_apex as f64,
                1.05,
                TRUE_APEX_COLOR,
                LineStyle::Solid,
                2.5,
            ));
            pui.line(vmarker(
                "",
                score.pass1.apex_cycle as f64,
                1.05,
                PASS1_COLOR,
                LineStyle::dashed_dense(),
                2.0,
            ));
            pui.line(vmarker(
                "",
                score.pass2.joint_apex_cycle as f64,
                1.05,
                PASS2_COLOR,
                LineStyle::Solid,
                2.0,
            ));
        });
}

/// Blue: the real peptide is present.
pub const PRESENT_COLOR: Color32 = Color32::from_rgb(31, 119, 180);
/// Red: the pure-noise twin of a present peptide.
pub const ABSENT_COLOR: Color32 = Color32::from_rgb(214, 39, 40);

/// Number of histogram bins over the log10 score range.
const N_BINS: usize = 48;

/// log10 of the finite, strictly-positive entries. A score of 0 (common for the
/// pure-noise twin) has no log and is dropped; callers report how many.
fn log10_positive(vals: &[f32]) -> Vec<f64> {
    vals.iter()
        .filter(|v| v.is_finite() && **v > 0.0)
        .map(|v| (*v as f64).log10())
        .collect()
}

/// Fraction-per-bin bars over `[lo, hi]`: heights sum to 1, so populations of
/// different size stay comparable. Returns the bars and the tallest height.
///
/// `hi` must exceed `lo` — an empty range would give zero-width bins, so
/// callers widen degenerate (all-equal) populations first.
fn hist_bars(logs: &[f64], lo: f64, hi: f64) -> (Vec<Bar>, f64) {
    debug_assert!(hi > lo, "empty bin range [{lo}, {hi}]");
    let bin_width = (hi - lo) / N_BINS as f64;
    let mut counts = vec![0usize; N_BINS];
    for &v in logs {
        let bin = (((v - lo) / bin_width).floor().max(0.0) as usize).min(N_BINS - 1);
        counts[bin] += 1;
    }
    let n_values = logs.len().max(1) as f64;
    let peak = counts.iter().max().copied().unwrap_or(0) as f64 / n_values;
    let bars = counts
        .iter()
        .enumerate()
        .map(|(bin, &count)| {
            Bar::new(lo + (bin as f64 + 0.5) * bin_width, count as f64 / n_values).width(bin_width)
        })
        .collect();
    (bars, peak)
}

/// One overlaid population in a score histogram.
pub struct HistSeries<'a> {
    pub name: &'a str,
    pub values: &'a [f32],
    pub color: Color32,
}

/// Overlaid log10 histograms of `main_score` populations (the score spans
/// orders of magnitude, so raw-scale bins are useless). Bars are fractions of
/// their own population, so series of different size stay comparable.
/// `marker_score` draws a vertical reference line. `id` must be unique per plot.
pub fn score_histogram(
    ui: &mut egui::Ui,
    id: &str,
    series: &[HistSeries<'_>],
    marker_score: Option<f32>,
) {
    let logs: Vec<Vec<f64>> = series.iter().map(|s| log10_positive(s.values)).collect();

    let dropped: usize = series
        .iter()
        .zip(&logs)
        .map(|(s, l)| s.values.len() - l.len())
        .sum();
    if dropped > 0 {
        ui.label(format!(
            "{dropped} value(s) scored <= 0 (not plottable in log10)"
        ));
    }
    if logs.iter().all(|l| l.is_empty()) {
        ui.label("(no positive scores to plot)");
        return;
    }

    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for &v in logs.iter().flatten() {
        lo = lo.min(v);
        hi = hi.max(v);
    }
    // Degenerate (all-equal) populations still need a non-zero bin width.
    if hi - lo < 1e-9 {
        lo -= 0.5;
        hi += 0.5;
    }

    let binned: Vec<(Vec<Bar>, f64)> = logs.iter().map(|l| hist_bars(l, lo, hi)).collect();
    let y_max = binned
        .iter()
        .map(|(_, peak)| *peak)
        .fold(1e-6_f64, f64::max);

    Plot::new(id)
        .legend(Legend::default())
        .height(200.0)
        .include_y(0.0)
        .include_y(y_max * 1.05)
        .show(ui, |pui| {
            // Reverse order: the first series is drawn last, so it stays on top.
            for (s, (bars, _)) in series.iter().zip(binned).rev() {
                pui.bar_chart(BarChart::new(s.name, bars).color(s.color));
            }
            if let Some(m) = marker_score.filter(|m| m.is_finite() && *m > 0.0) {
                pui.line(vmarker(
                    "current run",
                    (m as f64).log10(),
                    y_max * 1.05,
                    TRUE_APEX_COLOR,
                    LineStyle::Solid,
                    2.0,
                ));
            }
        });
}

/// Min-max normalize to [0,1]; constant series -> zeros.
fn normalize(ys: &[f32]) -> Vec<f32> {
    let mut lo = f32::INFINITY;
    let mut hi = f32::NEG_INFINITY;
    for &v in ys {
        if v.is_finite() {
            lo = lo.min(v);
            hi = hi.max(v);
        }
    }
    let range = hi - lo;
    if !range.is_finite() || range <= 0.0 {
        return vec![0.0; ys.len()];
    }
    ys.iter()
        .map(|&v| ((v - lo) / range).clamp(0.0, 1.0))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log10_positive_drops_nonpositive_and_nonfinite() {
        let got = log10_positive(&[100.0, 0.0, -1.0, f32::NAN, f32::INFINITY, 10.0]);
        assert_eq!(got, vec![2.0, 1.0]);
    }

    #[test]
    fn hist_bars_bins_everything_as_fractions() {
        let logs = vec![0.0, 0.5, 1.0]; // lowest, middle, and the top edge
        let (bars, peak) = hist_bars(&logs, 0.0, 1.0);
        let total: f64 = bars.iter().map(|b| b.value).sum();
        assert!(
            (total - 1.0).abs() < 1e-9,
            "bars must sum to 1, got {total}"
        );
        // The top-edge value lands in the last bin, not out of bounds.
        assert!(bars.last().unwrap().value > 0.0);
        assert!((peak - 1.0 / 3.0).abs() < 1e-9);
    }
}
