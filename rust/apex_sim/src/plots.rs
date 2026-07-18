//! egui rendering helpers: chromatogram overlay, heatmap, intermediate traces.

use eframe::egui::{
    self,
    Color32,
    ColorImage,
    TextureHandle,
    TextureOptions,
};
use egui_plot::{
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
        ui.colored_label(TRUE_APEX_COLOR, "TRUE (configured)");
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
