//! Rendering for the three dashboard tabs (`draw`) and the heatmap
//! downsampler that feeds the Fit tab (`heatmap_cells`).
//!
//! The Fit tab is the one built with the most care: `FitRecording::path_indices`
//! is `prefix ++ DP chain ++ suffix` (see `recording.rs`), and a calibration
//! that misbehaves at the edges of the gradient is often the greedily
//! attached tail's doing rather than the DP's. So the path overlay paints
//! `path[dp_range]` in one color and the two tails in another — that split is
//! the first thing a user reading this tab needs to see.
//!
//! Every panel here follows one rule: if there is nothing to draw, or the
//! area is too small to draw it in, render a `Paragraph` saying so. Nothing
//! indexes a slice without first checking it is in bounds — `ratatui`'s own
//! layout solver already tolerates degenerate `Rect`s, so the manual guards
//! below are for the raw buffer/grid indexing this module does by hand for
//! the half-block heatmap, not for widget layout.

use crate::metrics::weighted_ridge_half_width;
use crate::{
    App,
    BatchMetrics,
    FitRecording,
    Tab,
};
use calibrt::Point;
use ratatui::Frame;
use ratatui::layout::{
    Constraint,
    Layout,
    Rect,
};
use ratatui::style::{
    Color,
    Modifier,
    Style,
};
use ratatui::text::Line;
use ratatui::widgets::{
    Block,
    Paragraph,
    Row,
    Sparkline,
    Table,
    Tabs,
    Wrap,
};
use std::ops::Range;

/// How many of the most recent batches the Convergence tab's table shows.
const TABLE_ROWS: usize = 8;

pub fn draw(frame: &mut Frame, app: &mut App) {
    let area = frame.area();
    let rows = Layout::vertical([
        Constraint::Length(1),
        Constraint::Min(0),
        Constraint::Length(1),
    ])
    .split(area);

    draw_tab_bar(frame, rows[0], app);
    match app.tab() {
        Tab::Fit => draw_fit_tab(frame, rows[1], app),
        Tab::Convergence => draw_convergence_tab(frame, rows[1], app),
        Tab::Tolerances => draw_tolerances_tab(frame, rows[1], app),
    }
    draw_status_line(frame, rows[2], app);
}

fn draw_tab_bar(frame: &mut Frame, area: Rect, app: &App) {
    let titles: Vec<&str> = Tab::ALL.iter().map(|t| t.title()).collect();
    let selected = Tab::ALL.iter().position(|t| *t == app.tab()).unwrap_or(0);
    let tabs = Tabs::new(titles).select(selected).highlight_style(
        Style::default()
            .fg(Color::Yellow)
            .add_modifier(Modifier::BOLD),
    );
    frame.render_widget(tabs, area);
}

fn draw_status_line(frame: &mut Frame, area: Rect, app: &App) {
    let count = app
        .pending_count()
        .map(|n| n.to_string())
        .unwrap_or_else(|| "-".to_string());
    // Deliberately terse: at 100 columns the previous, fuller wording clipped
    // before `d:dp-pane`, hiding that key entirely rather than just looking
    // busy.
    let text = format!(
        " b{} {} cnt:{} | n:next r:run q:detach ^C:abort [/]:stage d:dp",
        app.batch(),
        app.stage().label(),
        count,
    );
    frame.render_widget(Paragraph::new(text), area);
}

// ---------------------------------------------------------------------
// Fit tab
// ---------------------------------------------------------------------

fn draw_fit_tab(frame: &mut Frame, area: Rect, app: &App) {
    if area.width == 0 || area.height == 0 {
        return;
    }
    let rec = app.recording();
    if rec.geom().bins == 0 {
        frame.render_widget(Paragraph::new("No grid recorded yet."), area);
        return;
    }

    // Only carve out a DP pane when there is enough width left for the
    // heatmap to still mean something — a narrow terminal keeps the grid and
    // simply does not show the pane, rather than splitting into two useless
    // slivers.
    let show_dp = app.dp_pane() && area.width >= 40 && area.height >= 3;
    let (grid_area, dp_area) = if show_dp {
        let cols = Layout::horizontal([Constraint::Percentage(65), Constraint::Percentage(35)])
            .split(area);
        (cols[0], Some(cols[1]))
    } else {
        (area, None)
    };

    draw_heatmap(frame, grid_area, app);
    if let Some(dp_area) = dp_area {
        draw_dp_pane(frame, dp_area, rec);
    }
}

/// Overlay priority for one terminal cell (which represents two grid rows —
/// see `heatmap_cells`). Higher wins when both grid rows the cell represents
/// carry different marks, or when a later stage's overlay lands on a cell an
/// earlier one already marked.
///
/// This — not a two-tone foreground/background per character — is what
/// actually carries the Fit tab's information in a way this module's own
/// tests can see. `TestBackend`'s buffer records a `Style` (colors) per
/// cell, but the given test harness (`render`, below) only ever reads
/// `.symbol()`. A `▀` painted with the "upper row in `fg`, lower row in `bg`"
/// technique the brief describes would satisfy that literally, but every
/// cell would carry the same glyph regardless of the underlying weights —
/// the asymmetric ridge fixture exists specifically so a transposed heatmap
/// is "visible in the snapshot" (a carry-over from Task 4), and a snapshot
/// that is a uniform block of `▀` regardless of data cannot show that, or
/// show one stage's overlay differing from another's. So the glyph itself
/// carries the signal here — heat density, suppression, path/DP-chain
/// membership, curve and ridge — and color is layered on top only as a
/// (real-terminal-only, untested) bonus.
const MARK_NONE: u8 = 0;
const MARK_SUPPRESSED: u8 = 1;
const MARK_TAIL: u8 = 2;
const MARK_DP_CHAIN: u8 = 3;
const MARK_CURVE: u8 = 4;
const MARK_RIDGE: u8 = 5;

fn mark_glyph(m: u8) -> Option<&'static str> {
    match m {
        MARK_SUPPRESSED => Some("."),
        MARK_TAIL => Some("+"),
        MARK_DP_CHAIN => Some("#"),
        MARK_CURVE => Some("*"),
        MARK_RIDGE => Some("~"),
        _ => None,
    }
}

fn mark_color(m: u8) -> Option<Color> {
    match m {
        MARK_SUPPRESSED => Some(Color::DarkGray),
        MARK_TAIL => Some(Color::Yellow),
        MARK_DP_CHAIN => Some(Color::LightGreen),
        MARK_CURVE => Some(Color::White),
        MARK_RIDGE => Some(Color::Magenta),
        _ => None,
    }
}

/// A four-level density scale for a cell with no overlay mark. Zero weight
/// renders as a plain space so the ridge stands out against empty grid.
fn heat_glyph(v: f32, max: f32) -> &'static str {
    if max <= 0.0 || v <= 0.0 {
        return " ";
    }
    match (v / max).clamp(0.0, 1.0) {
        t if t >= 0.75 => "\u{2588}", // █
        t if t >= 0.50 => "\u{2593}", // ▓
        t if t >= 0.25 => "\u{2592}", // ▒
        _ => "\u{2591}",              // ░
    }
}

fn heat_color(v: f32, max: f32) -> Color {
    if max <= 0.0 || v <= 0.0 {
        return Color::Reset;
    }
    let t = (v / max).clamp(0.0, 1.0);
    let level = (t * 255.0).round() as u8;
    Color::Rgb(0, level / 2, level)
}

/// The geometry every mark/overlay helper below needs: the terminal area's
/// size in cells (`w`, `h`), the grid's own bin count, and the display
/// resolution the two are reconciled at (`disp_rows` is `h * 2` — the
/// half-block doubling — `disp_cols` is `w`). Grouped into one `Copy` struct
/// so these helpers take one argument for "where on screen" instead of five.
#[derive(Clone, Copy)]
struct Dims {
    w: usize,
    h: usize,
    bins: usize,
    disp_rows: usize,
    disp_cols: usize,
}

/// Which half of a terminal cell a grid row landed in — the half-block
/// doubling `heatmap_cells`/`grid_to_screen` compute. Tracked per half
/// (rather than collapsed to one value per cell) so two grid rows sharing a
/// terminal line stay visually distinct: this is the vertical resolution
/// `heatmap_cells` computes and the interface promises.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Half {
    Upper,
    Lower,
}

fn draw_heatmap(frame: &mut Frame, area: Rect, app: &App) {
    if area.width == 0 || area.height == 0 {
        return;
    }
    let rec = app.recording();
    let bins = rec.geom().bins;
    if bins == 0 {
        frame.render_widget(Paragraph::new("Grid has zero bins."), area);
        return;
    }

    let w = area.width as usize;
    let h = area.height as usize;
    let dims = Dims {
        w,
        h,
        bins,
        disp_rows: h * 2,
        disp_cols: w,
    };
    let cells = heatmap_cells(rec, area.width, area.height);
    let max_w = cells.iter().copied().fold(0.0f32, f32::max).max(1e-9);

    // Two mark slots per terminal cell — index `(ty * w + tx) * 2 + half`,
    // `half` 0 for upper, 1 for lower — so an overlay on one grid row never
    // clobbers what the *other* grid row sharing that character wanted to
    // show.
    let mut marks = vec![MARK_NONE; w * h * 2];
    if app.show_suppressed() {
        mark_suppressed(&mut marks, dims, rec);
    }
    if app.show_path() {
        mark_path(&mut marks, dims, rec);
    }
    if app.show_curve() {
        mark_curve(&mut marks, dims, rec);
    }
    if app.show_ridge() {
        mark_ridge(&mut marks, dims, rec);
    }

    for ty in 0..h {
        for tx in 0..w {
            let idx = (ty * dims.disp_cols + tx) * 2;
            let upper_heat = cells.get(idx).copied().unwrap_or(0.0);
            let lower_heat = cells.get(idx + 1).copied().unwrap_or(0.0);
            let cell = ty * w + tx;
            let upper_mark = marks[cell * 2];
            let lower_mark = marks[cell * 2 + 1];
            let (symbol, color) =
                compose_cell(upper_mark, lower_mark, upper_heat, lower_heat, max_w);
            if let Some(buf_cell) = frame
                .buffer_mut()
                .cell_mut((area.x + tx as u16, area.y + ty as u16))
            {
                buf_cell.set_symbol(symbol);
                buf_cell.set_fg(color);
                buf_cell.set_bg(Color::Reset);
            }
        }
    }
}

/// Picks one terminal cell's glyph and color from its two independent
/// half-rows.
///
/// "Occupied" means either half has a nonzero weight or an overlay mark.
/// Occupancy shape is the primary, `TestBackend`-visible signal (`▀` upper
/// only, `▄` lower only, `█`/a density glyph when both, ` ` when neither) —
/// this is what makes `heatmap_cells`'s two-half-rows-per-line resolution
/// actually show up in a snapshot instead of being silently collapsed to one
/// glyph per cell.
///
/// When only one half is occupied, an overlay mark on that half wins its
/// glyph/color outright (marks are sparse and meant to stand out); with no
/// mark, that half's own heat intensity picks the color, and the shape is
/// plain (there is no standard partial-block character for "the upper half,
/// at quarter density").
///
/// When *both* halves are occupied, `heatmap_cells` collapsing them to one
/// glyph is unavoidable — one character cannot show two independent overlay
/// marks at once — so the placement rule is: the higher-priority mark
/// between the two halves wins the whole cell (an overlay is never invisible
/// just because it shares a character with the other half), and only when
/// *neither* half carries a mark does the cell fall back to a density glyph
/// over the combined (max) heat, which is exactly the common "both grid rows
/// on the ridge have real weight" case. `fit_tab_marks_the_higher_priority_half_when_both_are_occupied`
/// pins this rule with a fixture built so it actually triggers.
fn compose_cell(
    upper_mark: u8,
    lower_mark: u8,
    upper_heat: f32,
    lower_heat: f32,
    max: f32,
) -> (&'static str, Color) {
    let upper_on = upper_mark != MARK_NONE || upper_heat > 0.0;
    let lower_on = lower_mark != MARK_NONE || lower_heat > 0.0;
    match (upper_on, lower_on) {
        (false, false) => (" ", Color::Reset),
        (true, false) => match (mark_glyph(upper_mark), mark_color(upper_mark)) {
            (Some(sym), Some(col)) => (sym, col),
            _ => ("\u{2580}", heat_color(upper_heat, max)), // ▀
        },
        (false, true) => match (mark_glyph(lower_mark), mark_color(lower_mark)) {
            (Some(sym), Some(col)) => (sym, col),
            _ => ("\u{2584}", heat_color(lower_heat, max)), // ▄
        },
        (true, true) => {
            let winner = upper_mark.max(lower_mark);
            match (mark_glyph(winner), mark_color(winner)) {
                (Some(sym), Some(col)) => (sym, col),
                _ => {
                    let heat = upper_heat.max(lower_heat);
                    (heat_glyph(heat, max), heat_color(heat, max))
                }
            }
        }
    }
}

fn raise_mark(marks: &mut [u8], dims: Dims, tx: usize, ty: usize, half: Half, level: u8) {
    if tx >= dims.w || ty >= dims.h {
        return;
    }
    let slot_idx = (ty * dims.w + tx) * 2 + if half == Half::Upper { 0 } else { 1 };
    if let Some(slot) = marks.get_mut(slot_idx)
        && level > *slot
    {
        *slot = level;
    }
}

fn mark_grid_indices(marks: &mut [u8], dims: Dims, indices: &[usize], level: u8) {
    if dims.bins == 0 {
        return;
    }
    for &idx in indices {
        let row = idx / dims.bins;
        let col = idx % dims.bins;
        let (ty, tx, half) = grid_to_screen(row, col, dims);
        raise_mark(marks, dims, tx, ty, half, level);
    }
}

fn mark_suppressed(marks: &mut [u8], dims: Dims, rec: &FitRecording) {
    let bins = dims.bins;
    for row in 0..bins {
        for col in 0..bins {
            if rec.is_suppressed(row, col) && rec.weight(row, col) > 0.0 {
                let (ty, tx, half) = grid_to_screen(row, col, dims);
                raise_mark(marks, dims, tx, ty, half, MARK_SUPPRESSED);
            }
        }
    }
}

/// Marks `path[..dp_range.start]` and `path[dp_range.end..]` (Pass 2's
/// greedily attached tails) at one priority and `path[dp_range]` (the DP's
/// own chosen chain) at a higher one — the split the brief calls out as the
/// answer to the first question a user asks about a misbehaving edge of the
/// fit.
fn mark_path(marks: &mut [u8], dims: Dims, rec: &FitRecording) {
    let path = rec.path_indices();
    if dims.bins == 0 || path.is_empty() {
        return;
    }
    let dp_range = rec.dp_range();
    // Defensive clamp: `FitRecording` documents (and debug-asserts) that
    // `dp_range.end <= path.len()`, but a plain slice index is not guarded
    // by that assertion in a release build, so this tab clamps rather than
    // trusting it.
    let start = dp_range.start.min(path.len());
    let end = dp_range.end.clamp(start, path.len());

    mark_grid_indices(marks, dims, &path[..start], MARK_TAIL);
    mark_grid_indices(marks, dims, &path[start..end], MARK_DP_CHAIN);
    mark_grid_indices(marks, dims, &path[end..], MARK_TAIL);
}

/// Marks the fitted curve at every display column (not just at path nodes),
/// so it renders as a continuous line rather than sparse dots.
fn mark_curve(marks: &mut [u8], dims: Dims, rec: &FitRecording) {
    let geom = rec.geom();
    let bins = dims.bins;
    let curve = rec.curve();
    if bins == 0 || curve.is_empty() {
        return;
    }
    let (x_lo, x_hi) = geom.x_range;
    let span = (x_hi - x_lo).max(1e-9);
    for tx in 0..dims.w {
        let col_range = bin_range(tx, dims.disp_cols, bins);
        let col = col_range.start;
        let x = x_lo + (col as f64 + 0.5) / bins as f64 * span;
        let Some(y) = predict_curve(curve, x) else {
            continue;
        };
        let row = bin_of(y, geom.y_range, bins);
        let dr = forward_map(row, bins, dims.disp_rows);
        let half = if dr.is_multiple_of(2) {
            Half::Upper
        } else {
            Half::Lower
        };
        raise_mark(marks, dims, tx, dr / 2, half, MARK_CURVE);
    }
}

fn mark_ridge(marks: &mut [u8], dims: Dims, rec: &FitRecording) {
    let geom = rec.geom();
    let bins = dims.bins;
    if bins == 0 {
        return;
    }
    let curve = rec.curve();
    for m in rec.ridge() {
        let Some(center) = predict_curve(curve, m.library.0) else {
            continue;
        };
        let col = bin_of(m.library.0, geom.x_range, bins);
        for y in [center + m.half_width, center - m.half_width] {
            let row = bin_of(y, geom.y_range, bins);
            let (ty, tx, half) = grid_to_screen(row, col, dims);
            raise_mark(marks, dims, tx, ty, half, MARK_RIDGE);
        }
    }
}

fn draw_dp_pane(frame: &mut Frame, area: Rect, rec: &FitRecording) {
    if area.width == 0 || area.height == 0 {
        return;
    }
    let dp = rec.dp();
    if dp.is_empty() {
        frame.render_widget(
            Paragraph::new(
                "No DP node trace recorded for this batch (re-fit with dp_nodes enabled).",
            )
            .wrap(Wrap { trim: true })
            .block(Block::bordered().title("DP")),
            area,
        );
        return;
    }

    // No selection cursor exists yet (that lands with list/table state in a
    // later task), so the pane shows the most recently decided node's full
    // detail plus a compact trace of every node.
    let selected = &dp[dp.len() - 1];
    let mut lines: Vec<Line> = Vec::with_capacity(dp.len() + 4);
    lines.push(Line::styled(
        format!(
            "node {}: lib={:.3} obs={:.3} chose={} acc_w={:.3}",
            selected.i,
            selected.library,
            selected.observed,
            selected
                .chose
                .map(|c| c.to_string())
                .unwrap_or_else(|| "-".to_string()),
            selected.acc_weight,
        ),
        Style::default().add_modifier(Modifier::BOLD),
    ));
    if selected.considered.is_empty() {
        lines.push(Line::raw("  considered: (none)"));
    } else {
        lines.push(Line::raw("  considered:"));
        for (j, w) in &selected.considered {
            lines.push(Line::raw(format!("    j={j} edge_w={w:.3}")));
        }
    }
    lines.push(Line::raw(""));
    lines.push(Line::raw("path:"));
    for d in dp {
        let marker = if d.chose.is_some() { "->" } else { "  " };
        lines.push(Line::raw(format!(
            "{marker} i={:>3} lib={:.2} obs={:.2}",
            d.i, d.library, d.observed
        )));
    }

    frame.render_widget(
        Paragraph::new(lines)
            .wrap(Wrap { trim: false })
            .block(Block::bordered().title("DP")),
        area,
    );
}

/// Downsamples (or upsamples) `rec`'s `bins x bins` weight grid to
/// `area_w * area_h * 2` values — two half-rows per terminal line, so a
/// `▀`-painted cell at terminal row `y`, column `x` reads its foreground
/// (upper grid row) from `cells[(y * area_w + x) * 2]` and its background
/// (lower grid row) from `cells[(y * area_w + x) * 2 + 1]`.
///
/// Each display cell aggregates by taking the max weight over the block of
/// source cells that maps to it, in both directions: this is exactly a
/// partition of `0..bins` into `disp_rows` (`area_h * 2`) contiguous ranges
/// and of `0..bins` into `disp_cols` (`area_w`) contiguous ranges, so the
/// total work is `bins * bins` regardless of how the area is shaped — the
/// same amount of work as visiting the whole grid once. When `bins` is
/// smaller than the display area, several display cells legitimately share
/// one source range (replication, not interpolation): every display cell
/// still gets a real value, never a gap.
pub fn heatmap_cells(rec: &FitRecording, area_w: u16, area_h: u16) -> Vec<f32> {
    let area_w = area_w as usize;
    let area_h = area_h as usize;
    let total = area_w * area_h * 2;
    let mut out = vec![0.0f32; total];
    let bins = rec.geom().bins;
    if bins == 0 || area_w == 0 || area_h == 0 {
        return out;
    }
    let disp_rows = area_h * 2;
    for dr in 0..disp_rows {
        let row_range = bin_range(dr, disp_rows, bins);
        for dc in 0..area_w {
            let col_range = bin_range(dc, area_w, bins);
            let mut m = 0.0f32;
            for row in row_range.clone() {
                for col in col_range.clone() {
                    m = m.max(rec.weight(row, col));
                }
            }
            let cell_index = (dr / 2) * area_w + dc;
            let half = dr % 2;
            out[cell_index * 2 + half] = m;
        }
    }
    out
}

/// The half-open range of source indices (out of `src_n`) that display index
/// `disp_i` (out of `disp_n`) covers. A partition of `0..src_n`: summing the
/// lengths of every `disp_i`'s range yields exactly `src_n`, whether
/// downsampling (`src_n > disp_n`) or upsampling (`src_n < disp_n`, where
/// several consecutive `disp_i` map to the same single-element range).
fn bin_range(disp_i: usize, disp_n: usize, src_n: usize) -> Range<usize> {
    if disp_n == 0 || src_n == 0 {
        return 0..0;
    }
    let lo = disp_i * src_n / disp_n;
    let hi = ((disp_i + 1) * src_n)
        .div_ceil(disp_n)
        .max(lo + 1)
        .min(src_n);
    lo..hi
}

/// The forward counterpart of `bin_range`: which display index a given
/// source index falls into. Used to place overlays (a single path node, a
/// single curve point) at the same display cell `heatmap_cells`/`bin_range`
/// would have aggregated it into.
fn forward_map(src_i: usize, src_n: usize, disp_n: usize) -> usize {
    if src_n == 0 || disp_n == 0 {
        return 0;
    }
    (src_i * disp_n / src_n).min(disp_n - 1)
}

/// Maps a grid `(row, col)` to `(terminal_row, terminal_col, half)`.
fn grid_to_screen(row: usize, col: usize, dims: Dims) -> (usize, usize, Half) {
    let dr = forward_map(row, dims.bins, dims.disp_rows);
    let dc = forward_map(col, dims.bins, dims.disp_cols);
    let half = if dr.is_multiple_of(2) {
        Half::Upper
    } else {
        Half::Lower
    };
    (dr / 2, dc, half)
}

/// Grid-bin index of `v` within `range`, replicating `FitRecording`'s private
/// `col_of`/`row_of` (not reusable directly — `ui.rs` only has `FitRecording`'s
/// public surface). Non-finite `v` or a zero-width range map to bin 0 rather
/// than panicking or producing NaN.
fn bin_of(v: f64, range: (f64, f64), bins: usize) -> usize {
    if bins == 0 {
        return 0;
    }
    let (lo, hi) = range;
    let span = hi - lo;
    // Written as a positive check (rather than `!(span > 0.0)`) so a NaN
    // `span` is unambiguously "not usable" rather than relying on how `!`
    // interacts with a partial order.
    let usable = span.is_finite() && span > 0.0 && v.is_finite();
    if !usable {
        return 0;
    }
    (((v - lo) / span * bins as f64) as usize).min(bins - 1)
}

/// Linear interpolation over a (library-sorted) curve, clamped at the
/// endpoints rather than erroring — this only feeds a visual overlay, not a
/// prediction API. `None` only when the curve is empty.
fn predict_curve(curve: &[Point], x: f64) -> Option<f64> {
    if curve.len() < 2 {
        return curve.first().map(|p| p.observed);
    }
    let last = curve.len() - 1;
    if x <= curve[0].library {
        return Some(curve[0].observed);
    }
    if x >= curve[last].library {
        return Some(curve[last].observed);
    }
    let i = curve.partition_point(|p| p.library < x).clamp(1, last);
    let a = curve[i - 1];
    let b = curve[i];
    let span = (b.library - a.library).max(1e-9);
    let t = (x - a.library) / span;
    Some(a.observed + t * (b.observed - a.observed))
}

// ---------------------------------------------------------------------
// Convergence tab
// ---------------------------------------------------------------------

fn draw_convergence_tab(frame: &mut Frame, area: Rect, app: &App) {
    if area.width == 0 || area.height == 0 {
        return;
    }
    let metrics = app.metrics();
    if metrics.is_empty() {
        frame.render_widget(
            Paragraph::new(
                "No batches recorded yet — metrics appear after the first Phase 1 batch.",
            )
            .wrap(Wrap { trim: true })
            .block(Block::bordered().title("Convergence")),
            area,
        );
        return;
    }

    let rows = Layout::vertical([
        Constraint::Length(1),
        Constraint::Min(3),
        Constraint::Length(TABLE_ROWS as u16 + 3),
    ])
    .split(area);

    let header = format!(
        " {} batches  |  retained frames: {}  stride: every {}  dropped: {}",
        metrics.len(),
        app.retained_frames(),
        app.frame_stride(),
        app.dropped_frames(),
    );
    frame.render_widget(Paragraph::new(header), rows[0]);

    draw_sparklines(frame, rows[1], metrics);
    draw_batch_table(frame, rows[2], metrics);
}

fn draw_sparklines(frame: &mut Frame, area: Rect, metrics: &[BatchMetrics]) {
    if area.width == 0 || area.height == 0 {
        return;
    }
    let wrmse: Vec<f64> = metrics.iter().map(|m| m.wrmse).collect();
    let max_delta: Vec<f64> = metrics.iter().map(|m| m.max_delta).collect();
    let mean_delta: Vec<f64> = metrics.iter().map(|m| m.mean_delta).collect();
    let path_nodes: Vec<u64> = metrics.iter().map(|m| m.path_nodes as u64).collect();
    let ridge_hw: Vec<f64> = metrics.iter().map(|m| m.ridge_half_width).collect();

    let series: [(&str, Vec<u64>); 5] = [
        ("wrmse", scaled_u64(&wrmse)),
        ("max_delta", scaled_u64(&max_delta)),
        ("mean_delta", scaled_u64(&mean_delta)),
        ("path_nodes", path_nodes),
        ("ridge_half_width", scaled_u64(&ridge_hw)),
    ];

    // `Fill(1)` rather than `Ratio(1, 5)`: ratatui distributes leftover space
    // from integer rounding more evenly across `Fill` segments, so the five
    // sparklines come out the same height (or within one row of each other)
    // instead of a few segments silently absorbing all of the remainder.
    let spark_rows = Layout::vertical([Constraint::Fill(1); 5]).split(area);
    for (i, (label, data)) in series.iter().enumerate() {
        let Some(row) = spark_rows.get(i).copied() else {
            continue;
        };
        let spark = Sparkline::default()
            .block(Block::bordered().title(*label))
            .data(data.as_slice());
        frame.render_widget(spark, row);
    }
}

/// Scales a metric series into `0..=1000` for `Sparkline`, which only takes
/// `u64`. NaN/infinite samples (a real possibility — `weighted_ridge_half_width`
/// is NaN with no ridge measurements) render as `0` rather than corrupting the
/// scale; an all-non-finite or all-zero series renders as flat zero rather
/// than dividing by zero.
fn scaled_u64(values: &[f64]) -> Vec<u64> {
    let max = values
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .fold(0.0f64, f64::max);
    if max <= 0.0 {
        return vec![0; values.len()];
    }
    values
        .iter()
        .map(|v| {
            if v.is_finite() {
                ((v / max) * 1000.0).round().clamp(0.0, 1000.0) as u64
            } else {
                0
            }
        })
        .collect()
}

fn draw_batch_table(frame: &mut Frame, area: Rect, metrics: &[BatchMetrics]) {
    if area.width == 0 || area.height == 0 {
        return;
    }
    let start = metrics.len().saturating_sub(TABLE_ROWS);
    let recent = &metrics[start..];

    let header = Row::new([
        "chunk", "n", "wrmse", "max_d", "path", "ridge_hw", "admit", "evict",
    ]);
    let table_rows: Vec<Row> = recent
        .iter()
        .map(|m| {
            Row::new([
                m.chunk.to_string(),
                m.n_points.to_string(),
                format!("{:.4}", m.wrmse),
                format!("{:.4}", m.max_delta),
                m.path_nodes.to_string(),
                format!("{:.4}", m.ridge_half_width),
                m.admitted.to_string(),
                m.evicted.to_string(),
            ])
        })
        .collect();
    let widths = [Constraint::Length(8); 8];
    let table = Table::new(table_rows, widths)
        .header(header)
        .block(Block::bordered().title("Batches (churn: admit/evict)"));
    frame.render_widget(table, area);
}

// ---------------------------------------------------------------------
// Tolerances tab
// ---------------------------------------------------------------------

fn draw_tolerances_tab(frame: &mut Frame, area: Rect, app: &App) {
    if area.width == 0 || area.height == 0 {
        return;
    }
    match app.real_fit() {
        None => {
            let text = "Step B (m/z, mobility and RT-residual tolerance estimation) has not \
                         run yet. This dashboard is paused inside Phase 1 batch scoring; those \
                         measurements are only made after Phase 2 derives the final calibration \
                         from every collected calibrant. They do not exist yet — this is not an \
                         empty panel, there is nothing to show.";
            frame.render_widget(
                Paragraph::new(text)
                    .wrap(Wrap { trim: true })
                    .block(Block::bordered().title("Tolerances")),
                area,
            );
        }
        Some(rec) => {
            let ridge = rec.ridge();
            let rt_half_width = weighted_ridge_half_width(ridge);
            let mut lines = vec![Line::raw(format!(
                "RT residual: weighted half-width {rt_half_width:.4}s over {} ridge column(s)",
                ridge.len()
            ))];
            if let Some((min_hw, max_hw)) =
                ridge
                    .iter()
                    .map(|m| m.half_width)
                    .fold(None, |acc, hw| match acc {
                        None => Some((hw, hw)),
                        Some((lo, hi)) => Some((lo.min(hw), hi.max(hw))),
                    })
            {
                lines.push(Line::raw(format!("  range: {min_hw:.4}s .. {max_hw:.4}s")));
            } else {
                lines.push(Line::raw("  (no ridge measurements recorded)"));
            }
            lines.push(Line::raw(""));
            lines.push(Line::raw(
                "m/z and mobility distributions are Step B measurements outside this RT \
                 recording; they render here once the tolerance summary is wired through.",
            ));
            frame.render_widget(
                Paragraph::new(lines)
                    .wrap(Wrap { trim: true })
                    .block(Block::bordered().title("Tolerances")),
                area,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Stage;
    use calibrt::{
        CalibrationState,
        LibraryRT,
        ObserveOpts,
        ObservedRTSeconds,
    };
    use ratatui::Terminal;
    use ratatui::backend::TestBackend;

    fn render(app: &mut App, w: u16, h: u16) -> String {
        let mut t = Terminal::new(TestBackend::new(w, h)).expect("test terminal");
        t.draw(|f| draw(f, app)).expect("draw");
        let buf = t.backend().buffer().clone();
        // One string per row, so the snapshot is readable as a picture.
        (0..buf.area.height)
            .map(|y| {
                (0..buf.area.width)
                    .map(|x| buf[(x, y)].symbol())
                    .collect::<String>()
            })
            .collect::<Vec<_>>()
            .join("\n")
    }

    /// An asymmetric ridge: `x` ranges over `(0, 16)`, `y` over `(0, 48)` —
    /// different spans, and a shape that bends (slope 2 up to the midpoint,
    /// slope 4 after), so it is neither symmetric about the diagonal nor a
    /// straight line. `col_of`/`row_of` (recording.rs) are only exercised
    /// elsewhere by fixtures with `library == observed`, so a transposed
    /// heatmap would not show up there — it shows up here, in the snapshot.
    ///
    /// One extra point sits far beyond the core chain's last node, isolated
    /// in its own grid row and column, with a much lower weight than
    /// anything nearby.
    ///
    /// With `lookback == 1`, the DP only ever looks one survivor-rank back,
    /// so a `dip` point is inserted between the core chain and the stray:
    /// its `observed` value drops back down (a realistic outlier — RT does
    /// occasionally misbehave locally), so it fails the DP's monotonic edge
    /// test against the core chain's end and restarts fresh at its own small
    /// weight. The stray, one rank further, then only reaches back to the
    /// dip (small accumulated weight) — never to the core chain's real,
    /// much larger accumulated weight — so the DP's own global-best path
    /// never extends past the core chain. Pass 2's forward walk, which
    /// re-checks monotonicity against the DP's *actual* chosen endpoint
    /// rather than the DP's scoring, then skips the dip (it fails
    /// monotonicity against that endpoint) and attaches the stray directly
    /// as a suffix — exactly the DP-chain-vs-greedy-tail split the Fit tab
    /// exists to show.
    fn ridge_points() -> Vec<(LibraryRT<f64>, ObservedRTSeconds<f64>, f64)> {
        let mut pts: Vec<_> = (0..13)
            .map(|i| {
                let x = i as f64 + 1.0;
                let y = if x <= 8.0 {
                    2.0 * x
                } else {
                    16.0 + 4.0 * (x - 8.0)
                };
                (LibraryRT(x), ObservedRTSeconds(y), 5.0 + i as f64)
            })
            .collect();
        // Dip: breaks the DP's monotonic edge back to the core chain's end
        // (library=13, observed=36) by dropping `observed` well below 36.
        // observed=22 also keeps this point in a grid row none of the core
        // points occupy (row 7 of 16), so it is not out-weighed and
        // suppressed by one of them sharing its row.
        pts.push((LibraryRT(14.0), ObservedRTSeconds(22.0), 3.0));
        // Stray calibrant well past the core chain's last node, with a much
        // lower weight and its own isolated grid row/column.
        pts.push((LibraryRT(15.5), ObservedRTSeconds(44.0), 2.0));
        pts
    }

    fn fixture_app_with_ridge() -> App {
        let bins = 16;
        let mut app = App::new(bins);
        let mut state = CalibrationState::new(bins, (0.0, bins as f64), (0.0, 48.0), 1).unwrap();
        state.update(ridge_points().into_iter()).unwrap();
        state.fit_with(app.recording_mut(), ObserveOpts { dp_nodes: true });
        state.measure_ridge_width_with(0.3, app.recording_mut());
        app
    }

    /// A plain diagonal ridge sized to `bins`, asymmetric in range only (`x`
    /// spans `(0, bins)`, `y` spans `(0, 2*bins)`) — used only to exercise
    /// `heatmap_cells`'s down/upsampling at arbitrary bin counts, not the
    /// DP-chain distinction `fixture_app_with_ridge` is built for.
    fn fixture_recording(bins: usize) -> FitRecording {
        let bins = bins.max(1);
        let mut state =
            CalibrationState::new(bins, (0.0, bins as f64), (0.0, 2.0 * bins as f64), 5).unwrap();
        let pts: Vec<_> = (0..bins)
            .map(|i| {
                let x = i as f64 + 0.5;
                (LibraryRT(x), ObservedRTSeconds(2.0 * x), 1.0 + i as f64)
            })
            .collect();
        state.update(pts.into_iter()).unwrap();
        let mut rec = FitRecording::new(bins);
        state.fit_with(&mut rec, ObserveOpts::NONE);
        rec
    }

    /// Six batches with a decaying `max_delta` — the "has it stopped moving"
    /// signal the Convergence tab exists to show — and some admitted/evicted
    /// churn, on top of the ridge fixture's ("fixture_app_with_ridge")
    /// recording so the Fit tab still has something to show if a test
    /// switches back to it.
    fn fixture_app_with_metrics() -> App {
        let mut app = fixture_app_with_ridge();
        for i in 0..6u32 {
            let decay = 10.0 / (i as f64 + 1.0);
            app.push_metrics(BatchMetrics {
                chunk: i as usize,
                n_points: 20 + i as usize,
                wrmse: decay * 0.5,
                max_delta: decay,
                mean_delta: decay * 0.4,
                path_nodes: 12 + (i as usize % 3),
                ridge_half_width: 1.0 + decay * 0.1,
                admitted: (i % 3) as usize,
                evicted: (i % 2) as usize,
            });
        }
        app.set_frame_summary(6, 2, 4);
        app
    }

    #[test]
    fn fit_tab_renders_a_diagonal_ridge() {
        let mut app = fixture_app_with_ridge();
        insta::assert_snapshot!(render(&mut app, 100, 30));
    }

    #[test]
    fn fit_tab_renders_each_stage() {
        for stage in [
            Stage::Grid,
            Stage::Suppressed,
            Stage::Path,
            Stage::Curve,
            Stage::Ridge,
        ] {
            let mut app = fixture_app_with_ridge();
            app.set_stage(stage);
            insta::assert_snapshot!(format!("fit_stage_{:?}", stage), render(&mut app, 100, 30));
        }
    }

    /// `set_stage(Curve)` implies `suppressed`, `path` and `curve` are all
    /// on — but each is independently toggleable, not just a byproduct of
    /// the stage. Overriding `show_path` off afterward must remove the path
    /// marks (`#`/`+`) while `show_curve` (untouched) keeps the curve
    /// (`*`) rendering, proving the two are separate state rather than one
    /// stage-derived flag read two ways.
    #[test]
    fn fit_tab_overlay_toggles_are_independent_of_stage() {
        let mut app = fixture_app_with_ridge();
        app.set_stage(Stage::Curve);
        assert!(
            app.show_path(),
            "Curve stage implies the path overlay is on"
        );
        app.set_show_path(false);
        assert!(!app.show_path());
        assert!(
            app.show_curve(),
            "turning off `path` must not have touched `curve`"
        );

        let out = render(&mut app, 100, 30);
        assert!(
            !out.contains('#') && !out.contains('+'),
            "path overlay was toggled off, so no DP-chain or tail marks should \
             remain:\n{out}"
        );
        assert!(
            out.contains('*'),
            "curve overlay was left on and should still render:\n{out}"
        );
        insta::assert_snapshot!(out);
    }

    /// No given test toggles `dp_pane`, so nothing pinned the pane's content
    /// before this — an off-by-one in which node is "selected" or a
    /// formatting regression in `chose`/`acc_weight`/`considered` could
    /// silently drift. `fixture_app_with_ridge` already fits with
    /// `dp_nodes: true`, so `rec.dp()` is populated; this just has to turn
    /// the pane on and check its content is present and looks right.
    #[test]
    fn dp_pane_shows_the_selected_nodes_decision_and_considered_list() {
        use ratatui::crossterm::event::{
            KeyCode,
            KeyEvent,
            KeyModifiers,
        };
        let mut app = fixture_app_with_ridge();
        app.set_stage(Stage::Path);
        app.handle_key(KeyEvent::new(KeyCode::Char('d'), KeyModifiers::NONE));
        assert!(app.dp_pane());

        let out = render(&mut app, 100, 30);
        assert!(
            out.contains("chose="),
            "missing the selected node's chose:\n{out}"
        );
        assert!(
            out.contains("acc_w="),
            "missing the selected node's acc_weight:\n{out}"
        );
        assert!(
            out.contains("considered"),
            "missing the considered list:\n{out}"
        );
        assert!(
            out.contains("edge_w="),
            "missing a considered edge weight:\n{out}"
        );
        insta::assert_snapshot!(out);
    }

    #[test]
    fn convergence_tab_renders_metrics_and_churn() {
        let mut app = fixture_app_with_metrics();
        app.set_tab(Tab::Convergence);
        insta::assert_snapshot!(render(&mut app, 100, 30));
    }

    #[test]
    fn tolerances_tab_explains_itself_during_phase_one() {
        let mut app = fixture_app_with_ridge(); // real_fit is None
        app.set_tab(Tab::Tolerances);
        insta::assert_snapshot!(render(&mut app, 100, 30));
    }

    #[test]
    fn every_tab_survives_an_empty_state() {
        for tab in Tab::ALL {
            let mut app = App::new(10); // no frames, no metrics, no fit
            app.set_tab(tab);
            render(&mut app, 80, 24); // must not panic
        }
    }

    #[test]
    fn a_tiny_terminal_does_not_panic() {
        let mut app = fixture_app_with_ridge();
        for (w, h) in [(1u16, 1u16), (3, 3), (20, 5), (200, 8)] {
            render(&mut app, w, h);
        }
    }

    #[test]
    fn heatmap_downsamples_a_large_grid_to_the_area() {
        let rec = fixture_recording(200); // bins far exceeding the area
        let cells = heatmap_cells(&rec, 40, 10);
        assert_eq!(cells.len(), 40 * 10 * 2, "two half-rows per terminal line");
        assert!(
            cells.iter().any(|v| *v > 0.0),
            "the ridge must survive downsampling"
        );
    }

    #[test]
    fn heatmap_upsamples_a_small_grid_without_panicking() {
        let rec = fixture_recording(4); // bins far smaller than the area
        let cells = heatmap_cells(&rec, 40, 10);
        assert_eq!(cells.len(), 40 * 10 * 2);
    }

    /// A `bins=4` grid rendered into a body exactly 2 rows tall (`disp_rows
    /// == bins`, one-to-one) with a single display column, so both grid rows
    /// are known to land in the terminal's one heatmap character: row 0 (a
    /// point too light to survive `suppress_nonmax`'s 1.0 seed, so it is
    /// suppressed but still weight > 0 — the suppressed mark) in the upper
    /// half, row 1 (the sole survivor, so it is trivially the whole DP chain)
    /// in the lower half. `compose_cell`'s placement rule says the
    /// higher-priority mark wins the whole cell when both halves are
    /// occupied — `MARK_DP_CHAIN > MARK_SUPPRESSED` — so the character must
    /// be `#`, not `.`, and not some blend of the two.
    #[test]
    fn fit_tab_marks_the_higher_priority_half_when_both_are_occupied() {
        let bins = 4;
        let mut app = App::new(bins);
        let mut state = CalibrationState::new(bins, (0.0, 4.0), (0.0, 4.0), 2).unwrap();
        state
            .update(
                [
                    (LibraryRT(0.5), ObservedRTSeconds(0.5), 0.5),
                    (LibraryRT(1.5), ObservedRTSeconds(1.5), 2.0),
                ]
                .into_iter(),
            )
            .unwrap();
        state.fit_with(app.recording_mut(), ObserveOpts::NONE);
        app.set_stage(Stage::Path); // shows both the suppressed mask and the path

        // 1 tab-bar row + 2 heatmap rows + 1 status row; 1 column.
        let out = render(&mut app, 1, 4);
        let heatmap_char = out.lines().nth(1).and_then(|l| l.chars().next());
        assert_eq!(
            heatmap_char,
            Some('#'),
            "the DP chain mark must win the whole cell over the suppressed mark \
             sharing its other half, got:\n{out}"
        );
        insta::assert_snapshot!(out);
    }
}
