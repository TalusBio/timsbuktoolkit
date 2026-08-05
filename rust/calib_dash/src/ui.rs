//! Rendering for the three dashboard tabs (`draw`) and the heatmap
//! downsampler that feeds the Fit tab (`heatmap_cells`).
//!
//! `FitRecording::path_indices` is `prefix ++ DP chain ++ suffix` (see
//! `recording.rs`), and a misbehaving edge of the gradient is often the greedily
//! attached tail's doing, not the DP's. So the path overlay paints
//! `path[dp_range]` with one glyph and the two tails with another.
//!
//! # Inversion means mode, hue means mark kind
//!
//! `Modifier::REVERSED` says *what mode the screen is in*: the selected tab, a
//! scrubbed (not live) frame, an accumulating count, a marked heatmap cell.
//! Color says *what a thing is*: the density ramp and the mark kinds. Apart,
//! neither channel has to mean three things at once, and inversion — a swap of
//! whatever colors are already there — survives a colorblind reader, a
//! monochrome screenshot and any terminal theme. `Mark` is the only place that
//! adds to this.
//!
//! # Exactly one mark layer paints per frame
//!
//! `paint_heatmap` runs exactly one `mark_*` arm per draw (`crate::Layer`, cycled
//! by `m`/`M`), and each layer is internally overlap-free: the DP chain and greedy
//! tails partition the path by construction, and the other three layers only emit
//! `Mark::Region`. So no cross-layer priority exists to resolve below — `Mark`'s
//! `Ord` only arbitrates a DP node against a tail node rounded into one cell.
//!
//! Nothing to draw, or no room to draw it in, renders a `Paragraph` saying so
//! (`empty_panel`). The manual size guards are for the raw buffer indexing the
//! half-block heatmap does by hand; widget layout needs none of them.

use crate::ToleranceSummary;
use crate::app::{
    App,
    Layer,
    Tab,
};
use crate::metrics::BatchMetrics;
use crate::recording::{
    DpDecision,
    FitRecording,
};
use calibrt::{
    LibraryRT,
    RidgeSummary,
};
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
use ratatui::text::{
    Line,
    Span,
};
use ratatui::widgets::{
    Block,
    Clear,
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

/// Guards a division by a span that may be zero (a degenerate axis range, an
/// all-zero heat field, two curve points at the same library RT).
const EPS: f64 = 1e-9;

/// The `?` overlay's width: the help rows are `{keys:<11}{help}` plus two
/// border columns, and the longest help text lands well inside this.
const KEYS_OVERLAY_WIDTH: u16 = 52;

const DP_TITLE: &str = "DP";
const TOLERANCES_TITLE: &str = "Tolerances";
const FIT_TITLE: &str = "Fit";

pub(crate) fn draw(frame: &mut Frame, app: &mut App) {
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

    if app.show_keys() {
        draw_keys_overlay(frame, area);
    }
}

fn draw_tab_bar(frame: &mut Frame, area: Rect, app: &App) {
    let titles: Vec<&str> = Tab::ALL.iter().map(|t| t.title()).collect();
    let selected = Tab::ALL
        .iter()
        .position(|t| *t == app.tab())
        .expect("ALL lists every variant");
    let tabs = Tabs::new(titles)
        .select(selected)
        .highlight_style(Style::default().add_modifier(Modifier::REVERSED | Modifier::BOLD));
    frame.render_widget(tabs, area);
}

/// One binding, spelled once for both places it is advertised: `hint` is the
/// status line's terse label (it has to fit beside four others on a narrow
/// terminal), `help` the `?` overlay's spelled-out one. These tables keep the two
/// *renderers* in step; `App::handle_key` does its own routing.
#[derive(Clone, Copy)]
struct Binding {
    keys: &'static str,
    hint: &'static str,
    help: &'static str,
}

impl Binding {
    const fn new(keys: &'static str, hint: &'static str, help: &'static str) -> Self {
        Self { keys, hint, help }
    }

    /// One `?`-overlay row, keys in a fixed-width column so descriptions line up.
    fn help_line(&self) -> Line<'static> {
        Line::raw(format!("  {:<11}{}", self.keys, self.help))
    }
}

/// Bindings every tab answers to, in the order the `?` overlay lists them.
const GLOBAL_KEYS: &[Binding] = &[
    Binding::new("n", "next", "next batch"),
    Binding::new("r", "run", "run to end"),
    Binding::new("h l", "tab", "switch tab"),
    Binding::new("q", "detach", "detach"),
    Binding::new("^C", "abort", "abort"),
];

/// Where a tab's own bindings slot into `GLOBAL_KEYS` on the status line: after
/// the run controls, ahead of the tab bar and the two exits.
const TAB_KEYS_AT: usize = 2;

/// `?` itself, listed apart from `GLOBAL_KEYS` because the status line pins it
/// to its own right-hand column — see `fit_status_hints`.
const KEYS_OVERLAY_KEY: Binding = Binding::new("?", "keys", "this screen");

const FIT_KEYS: &[Binding] = &[
    Binding::new("< >", "frame", "scrub retained frames"),
    Binding::new("m M", "layer", "cycle mark layer"),
    Binding::new("d", "dp", "toggle DP pane"),
];

/// The tab-local bindings to advertise: the Fit tab's keys do nothing elsewhere.
fn tab_keys(tab: Tab) -> &'static [Binding] {
    if tab == Tab::Fit { FIT_KEYS } else { &[] }
}

/// One `key`/`action` hint per binding: the key `BOLD`, the action `DarkGray` —
/// weight is what makes the bound letter scannable. `actions` false renders the
/// keys alone, which is how the degrade below sheds action words but not keys.
fn binding_spans(bindings: &[Binding], actions: bool) -> Vec<Span<'static>> {
    let mut spans = Vec::new();
    for (i, b) in bindings.iter().enumerate() {
        if i > 0 {
            spans.push(Span::styled(
                " \u{b7} ",
                Style::default().fg(Color::DarkGray),
            ));
        }
        spans.push(Span::styled(
            b.keys.to_string(),
            Style::default().add_modifier(Modifier::BOLD),
        ));
        if actions && !b.hint.is_empty() {
            spans.push(Span::raw(" "));
            spans.push(Span::styled(
                b.hint.to_string(),
                Style::default().fg(Color::DarkGray),
            ));
        }
    }
    spans
}

/// The most detailed hint line (tab-local plus global bindings, ordered as
/// `TAB_KEYS_AT` describes) that still fits `width`, degrading in stages: full
/// text, keys only, global keys only, nothing. `? keys` is never part of this
/// line — it has its own pinned-right column, so it survives every stage.
fn fit_status_hints(tab_local: &[Binding], width: usize) -> Line<'static> {
    let (head, tail) = GLOBAL_KEYS.split_at(TAB_KEYS_AT);
    let full: Vec<Binding> = head.iter().chain(tab_local).chain(tail).copied().collect();
    for (bindings, actions) in [
        (full.as_slice(), true),
        (full.as_slice(), false),
        (GLOBAL_KEYS, false),
    ] {
        let line = Line::from(binding_spans(bindings, actions));
        if line.width() <= width {
            return line;
        }
    }
    Line::default()
}

/// The status line: batch/pending-count on the left, key hints in the middle
/// (degrading as `fit_status_hints` describes), `? keys` pinned to the right —
/// one `Paragraph` per column, so the right one stays right-aligned.
fn draw_status_line(frame: &mut Frame, area: Rect, app: &App) {
    let mut state_spans = vec![Span::raw(format!(" b{} ", app.batch()))];
    // `REVERSED` because the pending count is a mode.
    if let Some(n) = app.pending_count() {
        state_spans.push(Span::styled(
            n.to_string(),
            Style::default().add_modifier(Modifier::REVERSED),
        ));
        state_spans.push(Span::raw(" "));
    }
    let state = Line::from(state_spans);
    let state_w = (state.width() as u16).min(area.width);

    // Width measured from the binding, so renaming `keys` cannot silently clip the
    // one hint meant to survive every degrade. The extra column separates it.
    let keys_hint = Line::from(binding_spans(&[KEYS_OVERLAY_KEY], true));
    let right_w = (keys_hint.width() as u16 + 1).min(area.width);
    let mid_w = area.width.saturating_sub(state_w).saturating_sub(right_w) as usize;
    let mid_line = fit_status_hints(tab_keys(app.tab()), mid_w);

    let cols = Layout::horizontal([
        Constraint::Length(state_w),
        Constraint::Min(0),
        Constraint::Length(right_w),
    ])
    .split(area);

    frame.render_widget(Paragraph::new(state), cols[0]);
    frame.render_widget(Paragraph::new(mid_line), cols[1]);
    frame.render_widget(Paragraph::new(keys_hint), cols[2]);
}

/// The `?` overlay: the full key map, grouped by heading, drawn over whatever the
/// tab underneath was showing (`Clear` first, or the cells behind show through).
/// Dismissed on any key, so it needs no "press any key" footer.
fn draw_keys_overlay(frame: &mut Frame, area: Rect) {
    let heading = |s: &'static str| Line::styled(s, Style::default().add_modifier(Modifier::BOLD));
    let mut lines = vec![heading("Every tab")];
    lines.extend(GLOBAL_KEYS.iter().map(Binding::help_line));
    lines.push(KEYS_OVERLAY_KEY.help_line());
    lines.push(Line::raw(""));
    lines.push(heading("Fit tab"));
    lines.extend(FIT_KEYS.iter().map(Binding::help_line));

    // Height from the rows themselves plus two border rows: a literal height
    // silently clips the last group whenever a binding is added.
    let w = KEYS_OVERLAY_WIDTH.min(area.width);
    let h = (lines.len() as u16 + 2).min(area.height);
    if w == 0 || h == 0 {
        return;
    }
    let popup = Rect {
        x: area.x + (area.width - w) / 2,
        y: area.y + (area.height - h) / 2,
        width: w,
        height: h,
    };
    frame.render_widget(Clear, popup);
    frame.render_widget(
        Paragraph::new(lines).block(Block::bordered().title(" Keys ")),
        popup,
    );
}

// ---------------------------------------------------------------------
// Fit tab
// ---------------------------------------------------------------------

/// Minimum body size for the DP pane to be worth carving out: below this the
/// split leaves two useless slivers, so a narrow terminal keeps the whole grid.
const DP_PANE_MIN: (u16, u16) = (40, 3);

/// The DP pane's share of the Fit tab's width; the heatmap keeps the rest.
const DP_PANE_PERCENT: u16 = 35;

fn draw_fit_tab(frame: &mut Frame, area: Rect, app: &App) {
    if area.is_empty() {
        return;
    }
    let rec = app.active_recording();
    if rec.geom().bins == 0 {
        empty_panel(frame, area, FIT_TITLE, "No grid recorded yet.");
        return;
    }

    // The banner always wins the one row it needs, even when that leaves the
    // heatmap zero rows: `<`/`>` only work if it is never in doubt which frame is
    // on screen. An unreadably short heatmap beats one mistakable for live.
    let area = match app.scrub_frame() {
        Some(i) => {
            let rows = Layout::vertical([Constraint::Length(1), Constraint::Min(0)]).split(area);
            draw_scrub_banner(frame, rows[0], app, i);
            rows[1]
        }
        None => area,
    };

    let show_dp = app.dp_pane() && area.width >= DP_PANE_MIN.0 && area.height >= DP_PANE_MIN.1;
    let (grid_area, dp_area) = if show_dp {
        let cols = Layout::horizontal([
            Constraint::Percentage(100 - DP_PANE_PERCENT),
            Constraint::Percentage(DP_PANE_PERCENT),
        ])
        .split(area);
        (cols[0], Some(cols[1]))
    } else {
        (area, None)
    };

    draw_heatmap(frame, grid_area, app);
    if let Some(dp_area) = dp_area {
        // The pane needs a recording fit with `dp_nodes: true`. A scrubbed `rec`
        // always is one; the live batch fits with `ObserveOpts::NONE` to keep that
        // cost off the hot path, so it reads `sync_dp`'s on-demand recording,
        // falling back to `rec` (which renders "no DP trace").
        let dp_rec = if app.scrub_frame().is_some() {
            rec
        } else {
            app.live_dp_recording().unwrap_or(rec)
        };
        draw_dp_pane(frame, dp_area, dp_rec);
    }
}

/// The "you are not looking at live data" banner, `REVERSED` because it is a
/// mode: which retained frame is on screen, 1-based, and its batch when known.
fn draw_scrub_banner(frame: &mut Frame, area: Rect, app: &App, frame_index: usize) {
    let batch_note = app
        .scrub_chunk()
        .map(|c| format!(", batch {c}"))
        .unwrap_or_default();
    let text = format!(
        " SCRUBBED — retained frame {}/{}{batch_note} (not live; `>` returns to now)",
        frame_index + 1,
        app.frames().retained.max(1),
    );
    frame.render_widget(
        Paragraph::new(text)
            .style(Style::default().add_modifier(Modifier::REVERSED | Modifier::BOLD)),
        area,
    );
}

/// What one half of a terminal cell is marked with, for whichever layer
/// (`crate::Layer`) is active.
///
/// **Invariant:** a `Region` mark keeps the density glyph — `░▒▓█▀▄`,
/// half-occupancy intact — and is identified by inversion, so a mark and the
/// distribution under it stay both visible; color (`region_accent`) is redundant
/// reinforcement. `DpNode`/`Tail` are the exception: they *replace* the cell with
/// `O`/`X`, a path node being a point estimate rather than a region.
///
/// Declaration order *is* the resolution order (`Ord`) for the one collision that
/// exists: `DpNode` outranks `Tail` so the DP's chain shows through a tail sharing
/// its cell, and everything outranks `None`.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Mark {
    None,
    /// Curve, ridge or suppressed layer — the density glyph (or, in place of a
    /// space, `\u{b7}` — see `compose_marked`), reversed.
    Region,
    /// Path layer, Pass 2's greedily attached tail — `X`, the standard
    /// excluded/suspect-point convention.
    Tail,
    /// Path layer, the DP's own chosen chain — `O`.
    DpNode,
}

/// The active layer's accent color for a `Region` mark — reinforcement only (see
/// `Mark`), so a named ANSI constant (see `heat_color`) and never one of indices
/// 9-15 (`Light*`/`White`, invisible on a light background).
fn region_accent(layer: Layer) -> Color {
    match layer {
        Layer::Suppressed => Color::DarkGray,
        Layer::Curve => Color::Cyan,
        Layer::Ridge => Color::Magenta,
        // `Region` is never emitted while the active layer is `None`/`Path`, so
        // this arm exists only to keep the match exhaustive.
        Layer::None | Layer::Path => Color::Reset,
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

/// The density field's color. `░▒▓█` is already a luminance ramp, so this only
/// reinforces it with at most two gray steps, never `Rgb`: an `Rgb` color is not
/// remapped by a terminal's theme, so on a light background the ramp inverts and
/// the hottest cells become the *least* visible. The 0.5 cut matches
/// `heat_glyph`'s own mid-scale split.
fn heat_color(v: f32, max: f32) -> Color {
    if max <= 0.0 || v <= 0.0 {
        return Color::Reset;
    }
    let t = (v / max).clamp(0.0, 1.0);
    if t >= 0.5 {
        Color::Gray
    } else {
        Color::DarkGray
    }
}

/// The area's size in cells (`w`, `h`) plus the grid's bin count, in one arg.
#[derive(Clone, Copy)]
struct Dims {
    w: usize,
    h: usize,
    bins: usize,
}

impl Dims {
    /// The vertical display resolution: `h * 2`, the half-block doubling.
    fn disp_rows(&self) -> usize {
        self.h * 2
    }
}

/// Which half of a terminal cell a grid row landed in. Tracked per half so two
/// grid rows sharing a terminal line stay distinct — the doubled resolution
/// `heatmap_cells` promises.
#[derive(Clone, Copy)]
enum Half {
    Upper,
    Lower,
}

impl Half {
    /// Slot offset in a cell's pair of half-rows — the layout `heatmap_cells`'s
    /// `out` and the mark buffer share.
    fn slot(self) -> usize {
        usize::from(matches!(self, Half::Lower))
    }
}

/// Flips a display-row index so grid row 0 — the *lowest* observed RT — lands at
/// the *bottom* of the canvas. Every y flip in this module goes through here.
/// Without it grid row 0 would land at display row 0, the *top* one, and a
/// monotonically increasing calibration would render as a *descending* line.
///
/// The flip also flips which half of the doubled resolution the row falls into
/// (`dr` and `disp_rows - 1 - dr` have opposite parity, `disp_rows` always being
/// even here), so `Upper`/`Lower` must be re-derived from the flipped index.
fn flip_display_row(dr: usize, disp_rows: usize) -> (usize, Half) {
    let flipped = disp_rows.saturating_sub(1).saturating_sub(dr);
    let half = if flipped.is_multiple_of(2) {
        Half::Upper
    } else {
        Half::Lower
    };
    (flipped / 2, half)
}

/// Framed heatmap: a `Block` naming the axes in its title, a left gutter of
/// y-tick labels and an x-tick label row below, with `┤`/`┴` written into the
/// border itself. Falls back to painting straight into `area`, unframed.
fn draw_heatmap(frame: &mut Frame, area: Rect, app: &App) {
    let rec = app.active_recording();
    let geom = rec.geom();
    let Some(layout) = heatmap_layout(area, app) else {
        paint_heatmap(frame, area, rec, app);
        return;
    };

    let title_left = " Fit \u{2014} observed RT (s) \u{2191} vs library RT (s) \u{2192} ";
    // `batch` only — the active layer gets its own subtitle row below.
    let title_right = format!(" b{} ", app.batch());
    frame.render_widget(
        Block::bordered()
            .title_top(Line::from(title_left).left_aligned())
            .title_top(Line::from(title_right).right_aligned()),
        layout.block_area,
    );

    // `BOLD` so the subtitle does not read as one run with the titles above it.
    if let Some(subtitle_area) = layout.subtitle {
        let subtitle = fit_subtitle(app.layer(), subtitle_area.width as usize);
        frame.render_widget(
            Paragraph::new(subtitle).style(Style::default().add_modifier(Modifier::BOLD)),
            subtitle_area,
        );
    }
    paint_heatmap(frame, layout.canvas, rec, app);

    if let Some(gutter) = layout.gutter {
        draw_y_ticks(
            frame,
            gutter,
            layout.block_area.x,
            layout.canvas,
            geom.y_range,
            layout.y_target,
        );
    }
    if let Some(x_labels) = layout.x_labels {
        draw_x_ticks(
            frame,
            x_labels,
            layout.block_area,
            layout.inner_width,
            geom.x_range,
        );
    }
}

/// Where every piece of the framed heatmap goes: the block, the canvas painted
/// inside it, and the three areas — y gutter, x label row, subtitle — that are
/// `None` when the terminal has no room for them. `inner_width` is the canvas width
/// the x ticks are spaced across, `y_target` the tick count the gutter holds.
struct HeatmapLayout {
    block_area: Rect,
    gutter: Option<Rect>,
    x_labels: Option<Rect>,
    subtitle: Option<Rect>,
    canvas: Rect,
    inner_width: u16,
    y_target: usize,
}

/// `area` divided up for `draw_heatmap`, or `None` when it is too small to frame
/// at all and the heatmap is painted straight into it.
fn heatmap_layout(area: Rect, app: &App) -> Option<HeatmapLayout> {
    let (y_lo, y_hi) = app.active_recording().geom().y_range;

    // A framed canvas needs a `Block` with at least one interior row/column.
    // `show_x_axis` reserves one row below the block for x-tick labels.
    let show_x_axis = area.height as usize >= 4;
    let block_h = area.height as usize - usize::from(show_x_axis);
    if area.width < 3 || block_h < 3 {
        return None;
    }
    let canvas_h = block_h - 2;
    let y_target = y_tick_target(canvas_h);
    let gutter_w = y_gutter_width(y_lo, y_hi, y_target);
    // The gutter is only reserved with a few columns of canvas left over.
    let show_y_axis = area.width as usize >= gutter_w + 3;
    let gutter_w = if show_y_axis { gutter_w as u16 } else { 0 };

    let cols = Layout::horizontal([Constraint::Length(gutter_w), Constraint::Min(1)]).split(area);
    let (gutter_area, right_area) = (cols[0], cols[1]);
    let rows = Layout::vertical([
        Constraint::Min(1),
        Constraint::Length(u16::from(show_x_axis)),
    ])
    .split(right_area);
    let (block_area, x_label_area) = (rows[0], rows[1]);

    // Asked of a bare `Block::bordered()`: a `title_top` does not move the
    // interior while the top border is already taking that row.
    let inner = Block::bordered().inner(block_area);
    // The subtitle takes a row only when the canvas has one to spare.
    let (subtitle, canvas) = if inner.height >= 2 {
        let rows = Layout::vertical([Constraint::Length(1), Constraint::Min(0)]).split(inner);
        (Some(rows[0]), rows[1])
    } else {
        (None, inner)
    };

    Some(HeatmapLayout {
        block_area,
        gutter: (show_y_axis && canvas.height > 0).then_some(gutter_area),
        x_labels: (show_x_axis && inner.width > 0).then_some(x_label_area),
        subtitle,
        canvas,
        inner_width: inner.width,
        y_target,
    })
}

/// How many y-axis ticks a `canvas_h`-row canvas gets: one per five rows, never
/// fewer than two (one tick is not a scale) nor more than eight (labels collide).
fn y_tick_target(canvas_h: usize) -> usize {
    const ROWS_PER_TICK: usize = 5;
    (canvas_h / ROWS_PER_TICK).clamp(2, 8)
}

/// The x-axis counterpart. Spacing is far wider because an x label runs *along* its
/// axis, so neighbouring labels collide unless they are a label-width apart.
fn x_tick_target(canvas_w: usize) -> usize {
    const COLS_PER_TICK: usize = 20;
    (canvas_w / COLS_PER_TICK).clamp(2, 10)
}

/// The y-tick labels, right-aligned in the gutter, each with a `┤` written
/// into the block's left border at its row. `target` is passed in because the
/// gutter was already sized from it (`y_gutter_width`), and the two must agree.
fn draw_y_ticks(
    frame: &mut Frame,
    gutter: Rect,
    border_x: u16,
    canvas: Rect,
    (lo, hi): (f64, f64),
    target: usize,
) {
    let (_, decimals) = axis_scale(lo, hi, target);
    for v in axis_ticks(lo, hi, target) {
        let row = value_to_row(v, lo, hi, canvas.height as usize);
        let abs_y = canvas.y + row as u16;
        let label = format!("{v:>width$.decimals$}", width = gutter.width as usize);
        for (i, ch) in label.chars().enumerate() {
            let Some(cell) = frame.buffer_mut().cell_mut((gutter.x + i as u16, abs_y)) else {
                break;
            };
            cell.set_char(ch);
            cell.set_fg(Color::DarkGray);
        }
        if let Some(cell) = frame.buffer_mut().cell_mut((border_x, abs_y)) {
            cell.set_char('\u{2524}'); // ┤
        }
    }
}

/// The x-tick label row below the block, each label centered on its column,
/// with a `┴` written into the block's bottom border above it.
fn draw_x_ticks(
    frame: &mut Frame,
    label_area: Rect,
    block: Rect,
    canvas_w: u16,
    (lo, hi): (f64, f64),
) {
    let target = x_tick_target(canvas_w as usize);
    let (_, decimals) = axis_scale(lo, hi, target);
    let border_y = block.y + block.height - 1;
    let mut label_row = vec![' '; label_area.width as usize];
    for v in axis_ticks(lo, hi, target) {
        // Offset by the block's left border, which the canvas starts after.
        let block_col = 1 + value_to_index(v, lo, hi, canvas_w as usize);
        // Centered on its column and clipped: a label near the edge must not spill.
        let chars: Vec<char> = format!("{v:.decimals$}").chars().collect();
        let start = block_col.saturating_sub(chars.len() / 2);
        for (i, ch) in chars.into_iter().enumerate() {
            if let Some(slot) = label_row.get_mut(start + i) {
                *slot = ch;
            }
        }
        if let Some(cell) = frame
            .buffer_mut()
            .cell_mut((block.x + block_col as u16, border_y))
        {
            cell.set_char('\u{2534}'); // ┴
        }
    }
    let text: String = label_row.into_iter().collect();
    frame.render_widget(
        Paragraph::new(text).style(Style::default().fg(Color::DarkGray)),
        label_area,
    );
}

/// One phrase glossing what a layer's marks mean, appended to its name in
/// `fit_subtitle`. A `Region` mark is an inversion, not a distinct glyph, so there
/// is no glyph legend and this subtitle is the only place that can say what an
/// inverted cell means.
fn layer_gloss(layer: Layer) -> &'static str {
    match layer {
        Layer::None => "density only",
        Layer::Path => "O chosen, X greedy tail",
        Layer::Curve => "fitted calibration",
        Layer::Ridge => "tolerance band",
        Layer::Suppressed => "discarded by non-max suppression",
    }
}

/// The Fit heatmap's subtitle, spelled out as a label ("Showing: ridge —
/// tolerance band") rather than a key hint: a reader who does not know the active
/// layer cannot interpret a single inverted cell. Degrades fits-or-drops — gloss
/// first, then the label truncated — so it can never wrap or panic at `avail == 0`.
fn fit_subtitle(layer: Layer, avail: usize) -> String {
    let label = layer.label();
    let gloss = layer_gloss(layer);
    let full = format!(" Showing: {label} \u{2014} {gloss} ");
    if full.chars().count() <= avail {
        return full;
    }
    let no_gloss = format!(" Showing: {label} ");
    if no_gloss.chars().count() <= avail {
        return no_gloss;
    }
    label.chars().take(avail).collect()
}

/// Paints the half-block heatmap itself — density field plus the active mark
/// layer — into `area`: the framed canvas's inner rect, or the whole Fit-tab body
/// on a terminal too small to frame.
fn paint_heatmap(frame: &mut Frame, area: Rect, rec: &FitRecording, app: &App) {
    if area.is_empty() {
        return;
    }
    // The one zero-bins gate for everything the mark layers below reach: they
    // are only ever called from here, so none of them repeats it.
    let bins = rec.geom().bins;
    if bins == 0 {
        return;
    }

    let w = area.width as usize;
    let h = area.height as usize;
    let dims = Dims { w, h, bins };
    let cells = heatmap_cells(rec, area.width, area.height);
    let max_w = cells.iter().copied().fold(0.0f32, f32::max).max(1e-9);

    // Two mark slots per terminal cell, so a mark on one grid row never clobbers
    // what the *other* grid row sharing that character wanted to show.
    let layer = app.layer();
    let mut marks = vec![Mark::None; w * h * 2];
    match layer {
        Layer::None => {}
        Layer::Path => mark_path(&mut marks, dims, rec),
        Layer::Curve => mark_curve(&mut marks, dims, rec),
        Layer::Ridge => mark_ridge(&mut marks, dims, rec),
        Layer::Suppressed => mark_suppressed(&mut marks, dims, rec),
    }

    for ty in 0..h {
        for tx in 0..w {
            let idx = (ty * w + tx) * 2;
            let (symbol, color, modifier) = compose_cell(
                marks[idx],
                marks[idx + 1],
                cells[idx],
                cells[idx + 1],
                max_w,
                layer,
            );
            if let Some(buf_cell) = frame
                .buffer_mut()
                .cell_mut((area.x + tx as u16, area.y + ty as u16))
            {
                buf_cell.set_symbol(symbol);
                buf_cell.set_fg(color);
                buf_cell.set_bg(Color::Reset);
                buf_cell.modifier = modifier;
            }
        }
    }
}

/// A "nice" axis step (1/2/5 x 10^k) close to `span / target` — 37s at a target of
/// 8 ticks rounds to a step of 5s, not an unreadable 4.625.
fn nice_step(span: f64, target: usize) -> f64 {
    let raw = span / target.max(1) as f64;
    let mag = 10f64.powf(raw.log10().floor());
    let mult = match raw / mag {
        n if n <= 1.0 => 1.0,
        n if n <= 2.0 => 2.0,
        n if n <= 5.0 => 5.0,
        _ => 10.0,
    };
    mult * mag
}

/// Tick values in `[lo, hi]` spaced by `nice_step`, from the first step-aligned
/// value `>= lo`. Empty for a degenerate range or step.
fn axis_ticks(lo: f64, hi: f64, target: usize) -> Vec<f64> {
    let span = hi - lo;
    if !span.is_finite() || span <= 0.0 {
        return Vec::new();
    }
    let step = nice_step(span, target);
    if !step.is_finite() || step <= 0.0 {
        return Vec::new();
    }
    let first = (lo / step).ceil() * step;
    if !first.is_finite() {
        return Vec::new();
    }
    // Bounded so a pathological step can never loop unboundedly.
    let max_ticks = target.saturating_mul(3).clamp(4, 64);
    let mut out = Vec::with_capacity(max_ticks);
    let mut v = first;
    while v <= hi + step * 1e-6 && out.len() < max_ticks {
        out.push(v);
        v += step;
    }
    out
}

/// Decimal places to format an axis tick with, derived from the step size. Capped
/// at 2: RT seconds are never usefully shown finer than hundredths here.
fn axis_decimals(step: f64) -> usize {
    if !step.is_finite() || step <= 0.0 {
        return 0;
    }
    (-step.log10().floor()).clamp(0.0, 2.0) as usize
}

/// An axis's tick step and the precision to print its labels at — paired here
/// because tick placement, labels and the gutter sized for them must agree.
fn axis_scale(lo: f64, hi: f64, target: usize) -> (f64, usize) {
    let step = nice_step((hi - lo).abs().max(EPS), target);
    (step, axis_decimals(step))
}

/// Maps a data value to a 0-indexed slot of `n`, clamped to `0..n`. No flip: the
/// library-RT x axis runs left to right in screen order already.
fn value_to_index(v: f64, lo: f64, hi: f64, n: usize) -> usize {
    if n == 0 {
        return 0;
    }
    let span = hi - lo;
    if !span.is_finite() || span <= 0.0 {
        return 0;
    }
    let t = ((v - lo) / span).clamp(0.0, 1.0);
    ((t * n as f64) as usize).min(n - 1)
}

/// Maps a data value to a 0-indexed row inside a `canvas_h`-row canvas, y flipped
/// (see `flip_display_row`): the largest value lands at row 0 (screen top).
fn value_to_row(v: f64, lo: f64, hi: f64, canvas_h: usize) -> usize {
    if canvas_h == 0 {
        return 0;
    }
    canvas_h - 1 - value_to_index(v, lo, hi, canvas_h)
}

/// Left-gutter width for the y-tick labels: the wider of `lo`/`hi` at the tick
/// step's own precision, plus a column before the border, clamped.
fn y_gutter_width(lo: f64, hi: f64, target: usize) -> usize {
    let (_, decimals) = axis_scale(lo, hi, target);
    let lo_len = format!("{lo:.decimals$}").chars().count();
    let hi_len = format!("{hi:.decimals$}").chars().count();
    (lo_len.max(hi_len) + 1).clamp(4, 8)
}

/// Picks one terminal cell's glyph and color from its two independent half-rows.
/// "Occupied" means a half has nonzero weight or a mark, and the occupancy *shape*
/// (`▀` upper only, `▄` lower only, a density glyph for both, ` ` for neither) is
/// the primary signal — the one that makes the doubled vertical resolution visible
/// in a `.symbol()`-only snapshot. With both halves occupied one character cannot
/// show two marks, so `Mark`'s `Ord` picks the winner.
fn compose_cell(
    upper_mark: Mark,
    lower_mark: Mark,
    upper_heat: f32,
    lower_heat: f32,
    max: f32,
    layer: Layer,
) -> (&'static str, Color, Modifier) {
    let upper_on = upper_mark != Mark::None || upper_heat > 0.0;
    let lower_on = lower_mark != Mark::None || lower_heat > 0.0;
    match (upper_on, lower_on) {
        (false, false) => (" ", Color::Reset, Modifier::empty()),
        (true, false) => compose_marked(upper_mark, upper_heat, max, layer, "\u{2580}"), // ▀
        (false, true) => compose_marked(lower_mark, lower_heat, max, layer, "\u{2584}"), // ▄
        (true, true) => {
            let heat = upper_heat.max(lower_heat);
            let winner = upper_mark.max(lower_mark);
            compose_marked(winner, heat, max, layer, heat_glyph(heat, max))
        }
    }
}

/// Resolves one already-occupied half (or, from `compose_cell`'s both-occupied
/// case, the winning mark for the whole cell) into a glyph/color/modifier.
/// `none_glyph` is what to draw for `Mark::None` — the asymmetric `▀`/`▄` for a
/// single occupied half, or the combined density glyph when both are occupied and
/// neither carries a mark.
///
/// `O`/`X` replace the cell rather than inverting it and lean on `BOLD`,
/// `REVERSED` being reserved for mode and for `Mark::Region`.
///
/// A `Region` mark on a zero-weight cell must not keep `heat_glyph`'s plain space:
/// a *reversed* space is a solid block, the highest-density glyph here, and would
/// lie about the data — routine, since the curve is evaluated at every display
/// column and the ridge band often extends past the last real observation.
/// `\u{b7}` reversed reads instead as "mark here, nothing under it".
fn compose_marked(
    mark: Mark,
    heat: f32,
    max: f32,
    layer: Layer,
    none_glyph: &'static str,
) -> (&'static str, Color, Modifier) {
    match mark {
        Mark::DpNode => ("O", Color::Green, Modifier::BOLD),
        Mark::Tail => ("X", Color::Yellow, Modifier::BOLD),
        Mark::Region => {
            let glyph = match heat_glyph(heat, max) {
                " " => "\u{b7}",
                g => g,
            };
            (glyph, region_accent(layer), Modifier::REVERSED)
        }
        Mark::None => (none_glyph, heat_color(heat, max), Modifier::empty()),
    }
}

/// Raises one half-cell's mark to `level`, keeping whichever ranks higher.
fn raise_mark(marks: &mut [Mark], dims: Dims, tx: usize, ty: usize, half: Half, level: Mark) {
    if tx >= dims.w || ty >= dims.h {
        return;
    }
    let slot_idx = (ty * dims.w + tx) * 2 + half.slot();
    if let Some(slot) = marks.get_mut(slot_idx) {
        *slot = (*slot).max(level);
    }
}

fn mark_grid_indices(marks: &mut [Mark], dims: Dims, indices: &[usize], level: Mark) {
    for &idx in indices {
        let row = idx / dims.bins;
        let col = idx % dims.bins;
        let (ty, tx, half) = grid_to_screen(row, col, dims);
        raise_mark(marks, dims, tx, ty, half, level);
    }
}

fn mark_suppressed(marks: &mut [Mark], dims: Dims, rec: &FitRecording) {
    let bins = dims.bins;
    for row in 0..bins {
        for col in 0..bins {
            if rec.is_suppressed(row, col) && rec.weight(row, col) > 0.0 {
                let (ty, tx, half) = grid_to_screen(row, col, dims);
                raise_mark(marks, dims, tx, ty, half, Mark::Region);
            }
        }
    }
}

/// Marks the two greedily attached tails (`path[..dp_range.start]` and
/// `path[dp_range.end..]`, `Mark::Tail`) apart from the DP's own chain
/// (`path[dp_range]`, `Mark::DpNode`) — DP's choice, or a tail grafted on after?
/// The one layer that draws two mark kinds, so it relies on `Mark`'s `Ord` for the
/// rare cell downsampling rounds one of each into.
fn mark_path(marks: &mut [Mark], dims: Dims, rec: &FitRecording) {
    let path = rec.path_indices();
    if path.is_empty() {
        return;
    }
    let dp_range = rec.dp_range();
    // `FitRecording` only debug-asserts `dp_range.end <= path.len()`, which
    // buys a release build nothing, so clamp rather than trust it.
    let start = dp_range.start.min(path.len());
    let end = dp_range.end.clamp(start, path.len());

    mark_grid_indices(marks, dims, &path[..start], Mark::Tail);
    mark_grid_indices(marks, dims, &path[start..end], Mark::DpNode);
    mark_grid_indices(marks, dims, &path[end..], Mark::Tail);
}

/// Marks the fitted curve at every display column (not just at path nodes),
/// so it renders as a continuous line rather than sparse dots.
///
/// The sweep covers the whole grid x-range, wider than the curve's own span (whose
/// ends are path *cell centers*), so the outermost columns fall out of bounds and
/// stay unmarked: extrapolating or clamping flat would put ink where the fit has no
/// calibrant, and an overlay on the evidence heatmap is a claim about the evidence.
fn mark_curve(marks: &mut [Mark], dims: Dims, rec: &FitRecording) {
    let geom = rec.geom();
    let bins = dims.bins;
    let Some(curve) = rec.curve() else {
        return;
    };
    let (x_lo, x_hi) = geom.x_range;
    let span = (x_hi - x_lo).max(EPS);
    for tx in 0..dims.w {
        // Only the first source column this display column covers is evaluated.
        let col = bin_range(tx, dims.w, bins).start;
        let x = x_lo + (col as f64 + 0.5) / bins as f64 * span;
        let Ok(y) = curve.predict(LibraryRT(x)) else {
            continue;
        };
        let row = bin_of(y.0, geom.y_range, bins);
        let dr = forward_map(row, bins, dims.disp_rows());
        let (ty, half) = flip_display_row(dr, dims.disp_rows());
        raise_mark(marks, dims, tx, ty, half, Mark::Region);
    }
}

/// Brackets each measured column at `curve ± half_width`. Every measurement's
/// `library` is a path cell center, so it is one of the curve's own points and
/// never out of bounds — the `Err` arm is unreachable in practice.
fn mark_ridge(marks: &mut [Mark], dims: Dims, rec: &FitRecording) {
    let geom = rec.geom();
    let bins = dims.bins;
    let Some(curve) = rec.curve() else {
        return;
    };
    for m in rec.ridge() {
        let Ok(center) = curve.predict(m.library) else {
            continue;
        };
        let col = bin_of(m.library.0, geom.x_range, bins);
        for y in [center.0 + m.half_width, center.0 - m.half_width] {
            let row = bin_of(y, geom.y_range, bins);
            let (ty, tx, half) = grid_to_screen(row, col, dims);
            raise_mark(marks, dims, tx, ty, half, Mark::Region);
        }
    }
}

/// The panel a tab falls back to with nothing to show: the explanation inside
/// the same bordered, titled block the populated view uses, so it keeps its frame.
fn empty_panel(frame: &mut Frame, area: Rect, title: &'static str, text: &str) {
    frame.render_widget(
        Paragraph::new(text)
            .wrap(Wrap { trim: true })
            .block(Block::bordered().title(title)),
        area,
    );
}

fn draw_dp_pane(frame: &mut Frame, area: Rect, rec: &FitRecording) {
    let dp = rec.dp();
    if dp.is_empty() {
        empty_panel(
            frame,
            area,
            DP_TITLE,
            "No DP node trace recorded for this batch (re-fit with dp_nodes enabled).",
        );
        return;
    }

    frame.render_widget(
        Paragraph::new(dp_lines(dp))
            .wrap(Wrap { trim: false })
            .block(Block::bordered().title(DP_TITLE)),
        area,
    );
}

/// The DP pane's body: the last node's decision and the edges it weighed, then the
/// whole chain one line per node. `dp` must be non-empty.
fn dp_lines(dp: &[DpDecision]) -> Vec<Line<'static>> {
    let last = &dp[dp.len() - 1];
    let mut lines: Vec<Line> = Vec::with_capacity(dp.len() + 4);
    lines.push(Line::styled(
        format!(
            "last node: chose={} acc_w={:.3}",
            last.chose
                .map(|c| c.to_string())
                .unwrap_or_else(|| "-".to_string()),
            last.acc_weight,
        ),
        Style::default().add_modifier(Modifier::BOLD),
    ));
    if last.considered.is_empty() {
        lines.push(Line::raw("  considered: (none)"));
    } else {
        lines.push(Line::raw("  considered:"));
        for (j, w) in &last.considered {
            lines.push(Line::raw(format!("    j={j} edge_w={w:.3}")));
        }
    }
    lines.push(Line::raw(""));
    lines.push(Line::raw("path:"));
    for d in dp {
        // `root` is the absent predecessor — the chain's start, not a failure.
        let edge = match d.chose {
            Some(j) => format!("<-{j:<3}"),
            None => "root ".to_string(),
        };
        lines.push(Line::raw(format!(
            "i={:>3} lib={:.2} obs={:.2} {edge}",
            d.i, d.library, d.observed
        )));
    }
    lines
}

/// Resamples `rec`'s `bins x bins` weight grid to `area_w * area_h * 2` values: two
/// half-rows per terminal line, so the cell at row `y`, column `x` reads its upper
/// grid row from `cells[(y * area_w + x) * 2]` and its lower one from the next slot.
///
/// Each display cell takes the max weight over the source block mapping to it;
/// those blocks partition `0..bins` in both directions (`bin_range`), so the work
/// is `bins * bins` however the area is shaped, and upsampling replicates.
///
/// Rows are written to the *flipped* display position, through the same
/// `flip_display_row` `grid_to_screen` routes marks through, so the density field
/// and the marks cannot disagree about where a row lives.
fn heatmap_cells(rec: &FitRecording, area_w: u16, area_h: u16) -> Vec<f32> {
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
            let (ty, half) = flip_display_row(dr, disp_rows);
            out[(ty * area_w + dc) * 2 + half.slot()] = m;
        }
    }
    out
}

/// The half-open range of source indices (out of `src_n`) display index `disp_i`
/// (out of `disp_n`) covers — a partition of `0..src_n` whether downsampling or
/// upsampling (where consecutive `disp_i` share one single-element range).
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

/// The forward counterpart of `bin_range`: which display index a source index falls
/// into, so an overlay point lands in the cell `heatmap_cells` aggregated it into.
fn forward_map(src_i: usize, src_n: usize, disp_n: usize) -> usize {
    if src_n == 0 || disp_n == 0 {
        return 0;
    }
    (src_i * disp_n / src_n).min(disp_n - 1)
}

/// Maps a grid `(row, col)` to `(terminal_row, terminal_col, half)`, y flipped via
/// `flip_display_row` (which is also why the `Half` parity is the flipped index's).
fn grid_to_screen(row: usize, col: usize, dims: Dims) -> (usize, usize, Half) {
    let dr = forward_map(row, dims.bins, dims.disp_rows());
    let dc = forward_map(col, dims.bins, dims.w);
    let (ty, half) = flip_display_row(dr, dims.disp_rows());
    (ty, dc, half)
}

/// Grid-bin index of `v` within `range`, for placing an overlay mark on a value
/// that is not itself a grid cell — a curve prediction or a ridge bound.
///
/// Zero bins, a non-finite `v` or a zero-width range map to bin 0 rather than
/// panicking or producing NaN: overlay placement reads geometry straight off a
/// recording, so it cannot assume a caller ruled those out.
fn bin_of(v: f64, range: (f64, f64), bins: usize) -> usize {
    if bins == 0 {
        return 0;
    }
    let (lo, hi) = range;
    let span = hi - lo;
    // A positive check, not `!(span > 0.0)`: a NaN `span` is then unambiguously
    // "not usable" rather than depending on how `!` reads a partial order.
    if !(span.is_finite() && span > 0.0 && v.is_finite()) {
        return 0;
    }
    (((v - lo) / span * bins as f64) as usize).min(bins - 1)
}

// ---------------------------------------------------------------------
// Convergence tab
// ---------------------------------------------------------------------

/// One metric the Convergence tab reports, named once for both places it is
/// labelled — the sparkline title and the batch-table header, both on screen at
/// once. The short form wins: the table column is `TABLE_COL_WIDTH` wide.
/// `in_table`/`sparkline` select each view's subset — the table's bookkeeping
/// columns have nothing worth plotting, and `mean_d` is plotted with no column.
struct MetricColumn {
    label: &'static str,
    value: fn(&BatchMetrics) -> f64,
    in_table: bool,
    sparkline: bool,
    /// Whether a non-finite sample holds the previous value instead of dropping to
    /// zero (see `scaled_u64`), which the sparkline title has to say so a flat run
    /// is not misread as convergence. A count can never be non-finite.
    nan_holds: bool,
}

impl MetricColumn {
    const fn count(label: &'static str, value: fn(&BatchMetrics) -> f64) -> Self {
        Self {
            label,
            value,
            in_table: true,
            sparkline: false,
            nan_holds: false,
        }
    }

    const fn plotted_count(label: &'static str, value: fn(&BatchMetrics) -> f64) -> Self {
        Self {
            label,
            value,
            in_table: true,
            sparkline: true,
            nan_holds: false,
        }
    }

    const fn metric(label: &'static str, value: fn(&BatchMetrics) -> f64) -> Self {
        Self {
            label,
            value,
            in_table: true,
            sparkline: true,
            nan_holds: true,
        }
    }

    const fn plot_only(label: &'static str, value: fn(&BatchMetrics) -> f64) -> Self {
        Self {
            label,
            value,
            in_table: false,
            sparkline: true,
            nan_holds: true,
        }
    }
}

/// Every column, in the order both views render their own subset in.
const METRIC_COLUMNS: &[MetricColumn] = &[
    MetricColumn::count("chunk", |m| m.chunk as f64),
    MetricColumn::count("n", |m| m.n_points as f64),
    MetricColumn::metric("wrmse", |m| m.wrmse),
    MetricColumn::metric("max_d", |m| m.max_delta),
    MetricColumn::plot_only("mean_d", |m| m.mean_delta),
    MetricColumn::plotted_count("path", |m| m.path_nodes as f64),
    MetricColumn::metric("ridge_hw", |m| m.ridge_half_width),
    MetricColumn::count("admit", |m| m.admitted as f64),
    MetricColumn::count("evict", |m| m.evicted as f64),
];

fn draw_convergence_tab(frame: &mut Frame, area: Rect, app: &App) {
    let metrics = app.metrics();
    if metrics.is_empty() {
        empty_panel(
            frame,
            area,
            "Convergence",
            "No batches recorded yet — metrics appear after the first Phase 1 batch.",
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
        app.frames().retained,
        app.frames().keep_every,
        app.frames().dropped,
    );
    frame.render_widget(Paragraph::new(header), rows[0]);

    draw_sparklines(frame, rows[1], metrics);
    draw_batch_table(frame, rows[2], metrics);
}

fn table_columns() -> impl Iterator<Item = &'static MetricColumn> {
    METRIC_COLUMNS.iter().filter(|c| c.in_table)
}

fn draw_sparklines(frame: &mut Frame, area: Rect, metrics: &[BatchMetrics]) {
    // Every series is normalized to its own maximum and they share no unit, so no
    // shape or relative height here carries a magnitude — `spark_title` puts the
    // numbers back. `Fill(1)` spreads rounding leftovers across the panes.
    let columns: Vec<&MetricColumn> = METRIC_COLUMNS.iter().filter(|c| c.sparkline).collect();
    let spark_rows = Layout::vertical(vec![Constraint::Fill(1); columns.len()]).split(area);
    for (col, row) in columns.iter().zip(spark_rows.iter()) {
        let values: Vec<f64> = metrics.iter().map(col.value).collect();
        let data = scaled_u64(&values);
        let spark = Sparkline::default()
            .block(Block::bordered().title(spark_title(col.label, &values, col.nan_holds)))
            .data(data.as_slice());
        frame.render_widget(spark, *row);
    }
}

/// `label  peak <p>  now <n>` — the y-scale and the latest sample, which is what
/// turns a self-normalized shape back into a measurement: a `wrmse` that fell
/// 5.0 → 0.08 otherwise draws the same descent as one that fell 5.0 → 4.9. `now`
/// is the last *finite* sample, matching what `scaled_u64`'s hold actually draws.
fn spark_title(label: &str, values: &[f64], nan_holds: bool) -> String {
    let peak = values
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .fold(f64::NEG_INFINITY, f64::max);
    let now = values.iter().copied().rfind(|v| v.is_finite());
    let holds = if nan_holds { "  (NaN holds)" } else { "" };
    format!(
        "{label}  peak {}  now {}{holds}",
        fmt_metric(peak),
        fmt_metric(now.unwrap_or(f64::NAN)),
    )
}

/// A metric as a short decimal: `—` when non-finite, no fractional part for a
/// whole number (`path` is a count), scientific notation only outside `FIXED`.
fn fmt_metric(v: f64) -> String {
    // Outside this range four decimals print `0.0000` or run past the width.
    const FIXED: Range<f64> = 0.001..100_000.0;
    if !v.is_finite() {
        return "—".to_string();
    }
    let mag = v.abs();
    if v == 0.0 || v.fract() == 0.0 && mag < FIXED.end {
        format!("{v:.0}")
    } else if !FIXED.contains(&mag) {
        format!("{v:.2e}")
    } else {
        format!("{v:.4}")
    }
}

/// The resolution `scaled_u64` quantizes a series to. `Sparkline` scales a column's
/// height to the maximum sample, so this only has to be far finer than any pane.
const SPARK_SCALE: f64 = 1000.0;

/// Scales a metric series into `0..=SPARK_SCALE` for `Sparkline`, which takes
/// `u64` and so cannot mark a sample as "no data" distinctly from `0`.
/// Non-finite samples are routine here — batch 0 has no prior curve to diff
/// against, and a failed fit produces NaN throughout — and mapping them to `0` would
/// draw identically to "the curve stopped moving". So each holds the last finite
/// value, and only a leading run with no finite value yet reports as `0`.
fn scaled_u64(values: &[f64]) -> Vec<u64> {
    let max = values
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .fold(0.0f64, f64::max);
    if max <= 0.0 {
        return vec![0; values.len()];
    }
    let mut last = 0u64;
    values
        .iter()
        .map(|v| {
            if v.is_finite() {
                last = ((v / max) * SPARK_SCALE).round().clamp(0.0, SPARK_SCALE) as u64;
            }
            last
        })
        .collect()
}

/// Wide enough for the longest header (`ridge_hw`) and the widest `fmt_metric`.
const TABLE_COL_WIDTH: u16 = 8;

fn draw_batch_table(frame: &mut Frame, area: Rect, metrics: &[BatchMetrics]) {
    let start = metrics.len().saturating_sub(TABLE_ROWS);
    let recent = &metrics[start..];

    let header = Row::new(table_columns().map(|c| c.label));
    let table_rows: Vec<Row> = recent
        .iter()
        .map(|m| Row::new(table_columns().map(|c| fmt_metric((c.value)(m)))))
        .collect();
    let widths = vec![Constraint::Length(TABLE_COL_WIDTH); table_columns().count()];
    let table = Table::new(table_rows, widths)
        .header(header)
        .block(Block::bordered().title("Batches (churn: admit/evict)"));
    frame.render_widget(table, area);
}

// ---------------------------------------------------------------------
// Tolerances tab
// ---------------------------------------------------------------------

/// A signed tolerance window around zero: `±9.5` when the two sides match,
/// `-8.5 .. +9.5` when they do not. The explicit signs read as a window around each
/// calibrant, not a range of measured values, and collapsing the symmetric case
/// leaves asymmetry visible as a difference in shape.
fn fmt_interval((lo, hi): (f64, f64)) -> String {
    if lo == -hi {
        format!("±{hi:.1}")
    } else {
        format!("{lo:+.1} .. {hi:+.1}")
    }
}

fn draw_tolerances_tab(frame: &mut Frame, area: Rect, app: &App) {
    let Some(rec) = app.real_fit() else {
        empty_panel(
            frame,
            area,
            TOLERANCES_TITLE,
            "Step B (m/z, mobility and RT-residual tolerance estimation) has not run yet. It \
             runs after Phase 2 derives the final calibration from every collected calibrant.",
        );
        return;
    };
    // `trim: false`: the ridge sub-line's indent shows it qualifies the RT-residual
    // line above it, and `trim: true` strips exactly that indent.
    frame.render_widget(
        Paragraph::new(tolerance_lines(rec, app.tolerances()))
            .wrap(Wrap { trim: false })
            .block(Block::bordered().title(TOLERANCES_TITLE)),
        area,
    );
}

/// The Tolerances tab's body once Step B has run: the RT-residual half-width
/// measured off the ridge, then the m/z and mobility windows when wired through.
fn tolerance_lines(
    rec: &FitRecording,
    tolerances: Option<&ToleranceSummary>,
) -> Vec<Line<'static>> {
    let summary = RidgeSummary::of(rec.ridge());
    let mut lines = vec![Line::raw(format!(
        "RT residual: weighted half-width {:.4}s over {} ridge column(s)",
        summary.map_or(f64::NAN, |s| s.weighted_half_width),
        summary.map_or(0, |s| s.n_columns),
    ))];
    lines.push(match summary {
        None => Line::raw("  (no ridge measurements recorded)"),
        Some(s) => Line::raw(format!(
            "  range: {:.4}s .. {:.4}s",
            s.min_half_width, s.max_half_width
        )),
    });
    lines.push(Line::raw(""));
    match tolerances {
        Some(t) => {
            lines.push(Line::raw(format!(
                "m/z tolerance: {} ppm",
                fmt_interval(t.mz_ppm)
            )));
            lines.push(Line::raw(format!(
                "mobility tolerance: {} %",
                fmt_interval(t.mobility_pct)
            )));
            lines.push(Line::raw(format!(
                "RT tol (fallback): ±{:.1}s across {} calibrant(s)",
                t.rt_seconds, t.n_calibrants
            )));
        }
        None => lines.push(Line::raw(
            "m/z and mobility distributions are Step B measurements outside this RT recording; \
             they have not been wired through for this run.",
        )),
    }
    lines
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::frames::FrameSummary;
    use calibrt::{
        CalibrationState,
        ObserveOpts,
        ObservedRTSeconds,
    };
    use ratatui::Terminal;
    use ratatui::backend::TestBackend;
    use ratatui::buffer::Buffer;
    use ratatui::crossterm::event::{
        KeyCode,
        KeyEvent,
        KeyModifiers,
    };
    use std::num::NonZeroUsize;

    /// Five retained frames, nothing dropped, so `<`/`>` have frames to move over.
    const RETAINED_5: FrameSummary = FrameSummary {
        retained: 5,
        keep_every: 1,
        dropped: 0,
    };

    /// The width thresholds the layout branches at, and their neighbours: the
    /// unframed fallback (`area.width < 3`), the y-tick gutter's `gutter_w + 3`
    /// (`y_gutter_width` clamps to 4..=8, so 7 and 11 are the only widths that
    /// boundary sits at), `DP_PANE_MIN.0` and `KEYS_OVERLAY_WIDTH`.
    const SWEEP_WIDTHS: [u16; 12] = [0, 1, 2, 3, 6, 7, 11, 12, 39, 40, 51, 52];

    /// The height thresholds: two chrome rows come off before the body sees anything,
    /// then `area.height >= 4` (the x-tick row), `block_h < 3`, `inner.height >= 2`
    /// (the subtitle row) and `DP_PANE_MIN.1` — one row later again with the banner.
    const SWEEP_HEIGHTS: [u16; 7] = [0, 1, 2, 3, 4, 5, 6];

    /// Both dimensions inside their smallest thresholds, then where two thresholds
    /// cross.
    const SWEEP_TINY: [(u16, u16); 3] = [(0, 0), (1, 1), (2, 2)];
    const SWEEP_CROSSINGS: [(u16, u16); 4] = [(3, 4), (7, 4), (39, 3), (40, 3)];

    /// Larger terminals, including two extreme aspect ratios.
    const SWEEP_OVERSIZED: [(u16, u16); 5] = [(1, 40), (40, 1), (60, 16), (100, 30), (200, 60)];

    /// The short list, for bodies that are a single widget: the degenerate sizes
    /// plus one of each shape.
    const SWEEP_SMOKE: [(u16, u16); 7] =
        [(0, 0), (1, 1), (2, 2), (0, 12), (60, 0), (60, 3), (100, 30)];

    fn press(app: &mut App, c: char) {
        app.handle_key(KeyEvent::new(KeyCode::Char(c), KeyModifiers::NONE));
    }

    /// Cycles to `tab` with the real `l` binding — `App` has no jump-to-tab setter.
    fn goto_tab(app: &mut App, tab: Tab) {
        while app.tab() != tab {
            press(app, 'l');
        }
    }

    /// Cycles to `layer` with the real `m` binding.
    fn goto_layer(app: &mut App, layer: Layer) {
        while app.layer() != layer {
            press(app, 'm');
        }
    }

    /// The drawn glyphs alone, one string per row. Style-blind: use
    /// `render_snapshot` where inversion is part of what is being pinned.
    fn render(app: &mut App, w: u16, h: u16) -> String {
        glyph_grid(&draw_to_buffer(app, w, h))
    }

    /// `render` cropped to the tab body, for snapshots about a tab's content.
    fn render_body(app: &mut App, w: u16, h: u16) -> String {
        glyph_grid(&body_buffer(app, w, h))
    }

    /// The tab body alone — `draw`'s middle rows, without the tab bar or status line
    /// — rebased so row 0 is the body's first row. A tab-content snapshot has no
    /// business regenerating because a key hint was reworded.
    fn body_buffer(app: &mut App, w: u16, h: u16) -> Buffer {
        let buf = draw_to_buffer(app, w, h);
        let body = Rect::new(0, 0, buf.area.width, buf.area.height.saturating_sub(2));
        let mut out = Buffer::empty(body);
        for y in 0..body.height {
            for x in 0..body.width {
                out[(x, y)] = buf[(x, y + 1)].clone();
            }
        }
        out
    }

    fn draw_to_buffer(app: &mut App, w: u16, h: u16) -> Buffer {
        let mut t = Terminal::new(TestBackend::new(w, h)).expect("test terminal");
        t.draw(|f| draw(f, app)).expect("draw");
        t.backend().buffer().clone()
    }

    fn glyph_grid(buf: &Buffer) -> String {
        (0..buf.area.height)
            .map(|y| {
                (0..buf.area.width)
                    .map(|x| buf[(x, y)].symbol())
                    .collect::<String>()
            })
            .collect::<Vec<_>>()
            .join("\n")
    }

    /// The tab body's glyph grid plus a trailing section naming every `REVERSED`
    /// cell — for the snapshots where inversion *is* the payload: a region mark keeps
    /// the density glyph under it (`Mark`), so in a `.symbol()`-only grid a layer
    /// marking the wrong cells is byte-identical to one marking the right ones.
    /// Column spans keep it to one line per affected row.
    fn render_snapshot(app: &mut App, w: u16, h: u16) -> String {
        let buf = body_buffer(app, w, h);
        let mut out = glyph_grid(&buf);
        out.push_str("\n--- REVERSED cells (row: column spans) ---");
        let mut any = false;
        for y in 0..buf.area.height {
            let spans = reversed_spans(&buf, y);
            if spans.is_empty() {
                continue;
            }
            any = true;
            out.push_str(&format!("\n{y:>3}: {spans}"));
        }
        if !any {
            out.push_str("\n    (none)");
        }
        out
    }

    /// The inverted column runs on terminal row `y`, as `12-18` for a run and `7`
    /// for a single column, space separated.
    fn reversed_spans(buf: &Buffer, y: u16) -> String {
        let mut spans: Vec<String> = Vec::new();
        let mut start: Option<u16> = None;
        // `..=width`: the one-past-the-end step closes a run that reaches the
        // right edge.
        for x in 0..=buf.area.width {
            let reversed = x < buf.area.width && buf[(x, y)].modifier.contains(Modifier::REVERSED);
            match (reversed, start) {
                (true, None) => start = Some(x),
                (false, Some(lo)) => {
                    let hi = x - 1;
                    spans.push(if lo == hi {
                        format!("{lo}")
                    } else {
                        format!("{lo}-{hi}")
                    });
                    start = None;
                }
                _ => {}
            }
        }
        spans.join(" ")
    }

    /// An asymmetric ridge: `x` over `(0, 16)`, `y` over `(0, 48)`, bending at the
    /// midpoint (slope 2, then 4), so it is neither symmetric about the diagonal nor
    /// straight. Every other fixture has `library == observed`, so a transposed
    /// heatmap could only show up here.
    ///
    /// The trailing `dip` and `stray` points force the DP-chain versus greedy-tail
    /// split this tab shows: with `lookback == 1` the dip fails the DP's monotonic
    /// edge back to the core chain, so the DP's path never extends past it, while
    /// Pass 2's forward walk re-checks monotonicity against the DP's *chosen*
    /// endpoint, skips the dip and grafts the stray on as a suffix.
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
        // Dip: `observed` well below the core chain's end (36), and in a grid row
        // (7 of 16) no core point occupies, so nothing suppresses it.
        pts.push((LibraryRT(14.0), ObservedRTSeconds(22.0), 3.0));
        // Stray: past the chain's last node, low weight, own row and column.
        pts.push((LibraryRT(15.5), ObservedRTSeconds(44.0), 2.0));
        pts
    }

    fn fixture_app_with_ridge() -> App {
        let bins = 16;
        let mut app = App::new(bins);
        let mut state = CalibrationState::new(bins, (0.0, bins as f64), (0.0, 48.0), 1).unwrap();
        state.update(ridge_points().into_iter()).unwrap();
        state.fit_with(app.recording_mut(), ObserveOpts { dp_nodes: true });
        app.recording_mut().set_fit(&state);
        app
    }

    /// A plain diagonal ridge sized to `bins`, asymmetric in range only (`x` spans
    /// `(0, bins)`, `y` spans `(0, 2*bins)`) — only exercises `heatmap_cells`'s
    /// down/upsampling at arbitrary bin counts.
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

    /// Ten batches with a decaying `max_delta` and some churn, on top of
    /// `fixture_app_with_ridge` so the Fit tab still has something to draw.
    fn fixture_app_with_metrics() -> App {
        let mut app = fixture_app_with_ridge();
        // More batches than `TABLE_ROWS`, so the table's truncation is exercised.
        for i in 0..10u32 {
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
        app.set_frame_summary(FrameSummary {
            retained: 6,
            keep_every: 2,
            dropped: 4,
        });
        app
    }

    /// One snapshot per stop of the `m`/`M` cycle a picture can pin, each in
    /// isolation. `Layer::None` doubles as the plain Fit-tab render.
    ///
    /// `Layer::Suppressed` and `Layer::Ridge` are deliberately not pinned: Suppressed
    /// only marks cells that already carry weight, so its glyph grid is
    /// byte-identical to `Layer::None`'s; Ridge's payload is a *structure* a picture
    /// pins only by also re-pinning calibrt's blur kernel and
    /// `DEFAULT_RIDGE_FRACTION` across a crate boundary. Both are asserted at the
    /// mark-buffer level instead.
    #[test]
    fn fit_tab_renders_each_layer() {
        for layer in [Layer::None, Layer::Path, Layer::Curve] {
            let mut app = fixture_app_with_ridge();
            goto_layer(&mut app, layer);
            // Curve changes no glyph, so inversion is the only evidence it drew.
            let out = if layer == Layer::Curve {
                render_snapshot(&mut app, 100, 30)
            } else {
                render_body(&mut app, 100, 30)
            };
            insta::assert_snapshot!(format!("fit_layer_{:?}", layer), out);
        }
    }

    /// The active layer's name must survive on screen — without it a single inverted
    /// mark is uninterpretable. Checked on the two layers no snapshot pins, so a
    /// change to the degrade order cannot drop it for only some of them.
    #[test]
    fn fit_tab_subtitle_names_the_active_layer() {
        for layer in [Layer::Ridge, Layer::Suppressed] {
            let mut app = fixture_app_with_ridge();
            goto_layer(&mut app, layer);
            let out = render(&mut app, 100, 30);
            assert!(
                out.contains("Showing:") && out.contains(layer.label()),
                "expected the {layer:?} layer's subtitle to name it:\n{out}"
            );
        }
    }

    /// The DP pane's content: the header reports the *last* node's decision, the
    /// considered list one line per weighed edge, the path list one per node.
    #[test]
    fn dp_lines_report_the_last_nodes_decision_and_considered_list() {
        let app = fixture_app_with_ridge();
        let dp = app.recording().dp();
        assert!(
            dp.last().is_some_and(|last| !last.considered.is_empty()),
            "fixture drifted: the last node must have weighed at least one edge"
        );

        let text: Vec<String> = dp_lines(dp).iter().map(line_text).collect();
        insta::assert_snapshot!(text.join("\n"));
    }

    /// `App::set_scrub_recording` is what `CalibDash::sync_scrub` calls once `<`/`>`
    /// have moved `scrub_frame`. Pins that the Fit tab switches to that recording —
    /// a visibly different grid — and draws the "not live" banner.
    #[test]
    fn fit_tab_shows_a_banner_and_a_different_grid_when_scrubbing() {
        let mut app = fixture_app_with_ridge();
        app.set_frame_summary(RETAINED_5);
        let scrubbed = fixture_recording(8); // deliberately a different grid
        app.set_scrub_recording(2, 17, scrubbed);

        // The banner is inverted, so the `REVERSED` section is load-bearing here.
        insta::assert_snapshot!(render_snapshot(&mut app, 100, 30));
    }

    /// A Fit-tab body exactly one row tall must still show the banner rather than
    /// spend that row on a heatmap too short to read: without it the replayed grid is
    /// indistinguishable from live. And clearing the scrub must take the banner back
    /// off — a banner that survives the clear is the same lie inverted.
    #[test]
    fn the_scrub_banner_wins_the_only_body_row_and_goes_away_when_the_scrub_clears() {
        let mut app = fixture_app_with_ridge();
        app.set_frame_summary(RETAINED_5);
        app.set_scrub_recording(2, 17, fixture_recording(8));
        // 1 tab-bar row + 1 body row + 1 status row = 3.
        let out = render(&mut app, 100, 3);
        assert!(
            out.contains("SCRUBBED"),
            "the banner must win the body's only row rather than leave it \
             looking like an unlabeled (and possibly mistaken-for-live) heatmap:\n{out}"
        );

        app.clear_scrub();
        for (w, h) in [(100, 3), (100, 30)] {
            let out = render(&mut app, w, h);
            assert!(
                !out.contains("SCRUBBED"),
                "banner must be gone once scrub is cleared, at {w}x{h}:\n{out}"
            );
        }
    }

    /// Also pins `TABLE_ROWS`' cap: the fixture pushes more batches than the table
    /// shows, so the snapshot catches the window sliding to the wrong end.
    #[test]
    fn convergence_tab_renders_metrics_and_churn() {
        let mut app = fixture_app_with_metrics();
        assert!(
            app.metrics().len() > TABLE_ROWS,
            "the cap must be exercised"
        );
        goto_tab(&mut app, Tab::Convergence);
        insta::assert_snapshot!(render_body(&mut app, 100, 30));
    }

    /// The `?` overlay is the only place a key's spelled-out meaning appears, so a
    /// binding dropped from the table is invisible everywhere else. Also pins the
    /// `Clear`: the Fit tab's heatmap must not bleed through the overlay's rows.
    #[test]
    fn the_keys_overlay_lists_every_binding_over_the_tab_beneath_it() {
        let mut app = fixture_app_with_ridge();
        press(&mut app, '?');
        insta::assert_snapshot!(render(&mut app, 100, 30));
    }

    /// A `contains`, not a snapshot: the panel is literal prose with no computed
    /// content, so a picture of it pins four wrapped lines and rows of border padding
    /// to catch one thing — whether the explanation reaches the screen.
    #[test]
    fn tolerances_tab_explains_itself_during_phase_one() {
        let mut app = fixture_app_with_ridge(); // real_fit is None
        goto_tab(&mut app, Tab::Tolerances);
        let out = render(&mut app, 100, 12);
        assert!(
            out.contains(
                "Step B (m/z, mobility and RT-residual tolerance estimation) has not run yet."
            ),
            "the Tolerances tab must say why it has nothing to show in Phase 1:\n{out}"
        );
    }

    /// Once `App::set_final` has run, the Tolerances tab must render the m/z,
    /// mobility and RT numbers rather than the "not wired through yet" placeholder.
    #[test]
    fn tolerances_tab_renders_the_summary_once_set_final_has_run() {
        let mut app = fixture_app_with_ridge();
        let rec = {
            let mut state = CalibrationState::new(16, (0.0, 16.0), (0.0, 48.0), 1).unwrap();
            state.update(ridge_points().into_iter()).unwrap();
            let mut rec = FitRecording::new(16);
            state.fit_with(&mut rec, ObserveOpts::NONE);
            rec.set_fit(&state);
            rec
        };
        app.set_final(
            rec,
            ToleranceSummary {
                mz_ppm: (-8.5, 9.5),
                mobility_pct: (-3.0, 3.0),
                rt_seconds: 12.5,
                n_calibrants: 42,
            },
        );
        goto_tab(&mut app, Tab::Tolerances);
        insta::assert_snapshot!(render_body(&mut app, 100, 12));
    }

    /// The half-block parity, end to end and independently derived.
    /// `fixture_recording(4)`'s four points land one per grid cell on the diagonal —
    /// `(row, col)` = `(0,0) (1,1) (2,2) (3,3)` with weights 1..4 — so a 4x2 canvas
    /// is exactly two grid rows per terminal line.
    ///
    /// Row 0 (the *lowest* observed RT) must land in the bottom line's *lower* half
    /// and row 3 in the top line's *upper* half; swapping `Half`'s two slots keeps
    /// every glyph a valid block character and moves nothing but this.
    #[test]
    fn heatmap_cells_flips_grid_row_zero_into_the_bottom_lines_lower_half() {
        let rec = fixture_recording(4);
        for (row, col) in [(0, 0), (1, 1), (2, 2), (3, 3)] {
            assert_eq!(
                rec.weight(row, col),
                1.0 + row as f32,
                "fixture drifted: expected weight {} at grid ({row}, {col})",
                1.0 + row as f32
            );
        }

        // One line per terminal row, each `(upper, lower)` pair a cell: line 0 holds
        // grid rows 3 and 2, line 1 holds grid rows 1 and 0.
        #[rustfmt::skip]
        let expected = [
            0.0, 0.0,   0.0, 0.0,   0.0, 3.0,   4.0, 0.0,
            0.0, 1.0,   2.0, 0.0,   0.0, 0.0,   0.0, 0.0,
        ];
        assert_eq!(heatmap_cells(&rec, 4, 2), expected);
    }

    /// The other half of the same parity: a cell occupied only in its upper half must
    /// render as `▀` and one occupied only in its lower half as `▄`. Without this the
    /// only detector of a swapped pair is `fit_layer_None`'s snapshot.
    #[test]
    fn compose_cell_draws_an_upper_only_cell_as_a_top_half_block() {
        let upper = compose_cell(Mark::None, Mark::None, 1.0, 0.0, 1.0, Layer::None);
        assert_eq!(upper.0, "\u{2580}", "upper-half-only must draw ▀");
        let lower = compose_cell(Mark::None, Mark::None, 0.0, 1.0, 1.0, Layer::None);
        assert_eq!(lower.0, "\u{2584}", "lower-half-only must draw ▄");
    }

    /// `Mark`'s ordering through `compose_cell`, exercised directly.
    #[test]
    fn compose_cell_gives_the_dp_chain_mark_priority_over_a_tail_mark_sharing_its_cell() {
        let (glyph, color, modifier) =
            compose_cell(Mark::Tail, Mark::DpNode, 0.0, 0.0, 1.0, Layer::Path);
        assert_eq!(
            glyph, "O",
            "the DP chain mark must win over a tail mark sharing the cell"
        );
        assert_eq!(color, Color::Green);
        assert_eq!(modifier, Modifier::BOLD);

        // Order-independence: not which half happened to carry which mark.
        let (glyph2, ..) = compose_cell(Mark::DpNode, Mark::Tail, 0.0, 0.0, 1.0, Layer::Path);
        assert_eq!(glyph2, "O");
    }

    /// The empty-cell trap `compose_marked` describes: a reversed space renders as a
    /// solid block on a cell with no weight at all, so `\u{b7}` must substitute.
    #[test]
    fn compose_cell_never_renders_a_marked_zero_weight_cell_as_a_reversed_space() {
        let (glyph, _color, modifier) =
            compose_cell(Mark::Region, Mark::None, 0.0, 0.0, 1.0, Layer::Curve);
        assert_ne!(
            glyph, " ",
            "a reversed space reads as maximum density and lies about a \
             zero-weight cell"
        );
        assert_eq!(glyph, "\u{b7}");
        assert!(modifier.contains(Modifier::REVERSED));

        let (glyph_both, _color2, modifier_both) =
            compose_cell(Mark::Region, Mark::Region, 0.0, 0.0, 1.0, Layer::Ridge);
        assert_ne!(glyph_both, " ");
        assert_eq!(glyph_both, "\u{b7}");
        assert!(modifier_both.contains(Modifier::REVERSED));
    }

    /// A `Mark::Region` on a weighted cell changes no glyph, so neither function's
    /// output is visible in a `.symbol()`-only snapshot — this pins what they mark, at
    /// the mark-buffer level. `Layer::Ridge`'s whole payload is the sparse
    /// `center ± half_width` pair per measured column, straddling the curve, rather
    /// than the filled band between them; asserted structurally because a snapshot
    /// would re-pin calibrt's blur kernel across a crate boundary.
    #[test]
    fn mark_suppressed_marks_a_weighted_cell_and_mark_ridge_marks_a_straddling_pair() {
        let app = fixture_app_with_ridge();
        let rec = app.recording();
        let geom = rec.geom();
        let dims = Dims {
            w: 100,
            h: 15,
            bins: geom.bins,
        };

        let mut suppressed = vec![Mark::None; dims.w * dims.h * 2];
        mark_suppressed(&mut suppressed, dims, rec);
        assert!(
            suppressed.contains(&Mark::Region),
            "expected at least one suppressed cell in the fixture"
        );

        let mut ridge = vec![Mark::None; dims.w * dims.h * 2];
        mark_ridge(&mut ridge, dims, rec);

        // The screen half-rows marked in display column `tx`, top to bottom: slot 0
        // is the upper half, so half-row `sr` lives in cell `sr / 2` at slot `sr % 2`.
        let marked_rows = |tx: usize| -> Vec<usize> {
            (0..dims.disp_rows())
                .filter(|sr| ridge[((sr / 2) * dims.w + tx) * 2 + sr % 2] != Mark::None)
                .collect()
        };
        let columns: Vec<(usize, Vec<usize>)> = (0..dims.w)
            .map(|tx| (tx, marked_rows(tx)))
            .filter(|(_, rows)| !rows.is_empty())
            .collect();

        assert!(
            !columns.is_empty(),
            "expected at least one ridge-band cell in the fixture"
        );
        for (tx, rows) in &columns {
            assert!(
                rows.len() <= 2,
                "column {tx} carries {} marks: `center ± half_width` is a pair \
                 of edges, not the filled interval between them: {rows:?}",
                rows.len()
            );
        }

        // And the pair brackets the curve: dropping either the `+` or the `-` term,
        // or collapsing both to `center`, leaves every assertion above green.
        let straddling = rec
            .ridge()
            .iter()
            .filter(|m| {
                let Ok(center) = rec.curve().unwrap().predict(m.library) else {
                    return false;
                };
                let col = bin_of(m.library.0, geom.x_range, dims.bins);
                let row = bin_of(center.0, geom.y_range, dims.bins);
                let (ty, tx, half) = grid_to_screen(row, col, dims);
                let center_sr = ty * 2 + half.slot();
                let rows = marked_rows(tx);
                rows.iter().any(|&sr| sr < center_sr) && rows.iter().any(|&sr| sr > center_sr)
            })
            .count();
        assert!(
            straddling > 0,
            "expected at least one measured column whose marks land both above \
             and below the curve row they bracket"
        );
    }

    /// Below `DP_PANE_MIN` the pane must be *hidden* — the heatmap keeps the whole
    /// body — which is what `>` for `>=` in either comparison would break.
    #[test]
    fn a_body_narrower_than_dp_pane_min_hides_the_pane_instead_of_splitting() {
        let mut app = fixture_app_with_ridge();
        press(&mut app, 'd');
        assert!(
            app.dp_pane(),
            "the pane must be requested for this to mean anything"
        );

        // One `┌` when the heatmap has the body to itself, two once the pane is out.
        let blocks_at = |app: &mut App, w: u16| -> usize {
            render(app, w, 30)
                .lines()
                .nth(1)
                .expect("the Fit tab's block starts on the row under the tab bar")
                .matches('\u{250c}')
                .count()
        };

        assert_eq!(
            blocks_at(&mut app, DP_PANE_MIN.0),
            2,
            "at exactly DP_PANE_MIN width the pane must be carved out"
        );
        assert_eq!(
            blocks_at(&mut app, DP_PANE_MIN.0 - 1),
            1,
            "one column below DP_PANE_MIN the heatmap must keep the whole body, \
             not share it with a sliver of a DP pane"
        );
    }

    // ---- axis helpers: a plausible-but-wrong axis is a snapshot diff a reviewer
    // waves through, so every expected value here is derived by hand. -------------

    /// The tick densities, which nothing else here can see: `nice_step` rounds so hard
    /// that changing `ROWS_PER_TICK` from 5 to 4 leaves every snapshot byte-identical.
    #[test]
    fn tick_targets_scale_with_the_canvas_and_stay_inside_their_clamps() {
        assert_eq!(y_tick_target(20), 4);
        assert_eq!(y_tick_target(25), 5);
        assert_eq!(y_tick_target(0), 2);
        assert_eq!(y_tick_target(9), 2);
        assert_eq!(y_tick_target(1000), 8);
        assert_eq!(x_tick_target(40), 2);
        assert_eq!(x_tick_target(100), 5);
        assert_eq!(x_tick_target(0), 2);
        assert_eq!(x_tick_target(1000), 10);
    }

    #[test]
    fn nice_step_rounds_to_one_two_or_five_times_a_power_of_ten() {
        // The doc comment's own example: 4.625 a tick must read as 5s.
        assert_eq!(nice_step(37.0, 8), 5.0);
        assert_eq!(nice_step(100.0, 10), 10.0);
        // Sub-unit spans keep the 1/2/5 shape one decade down.
        assert_eq!(nice_step(0.8, 4), 0.2);
        assert_eq!(nice_step(0.4, 4), 0.1);
        // `target` of zero is treated as one rather than dividing by it.
        assert_eq!(nice_step(5.0, 0), 5.0);
    }

    #[test]
    fn axis_ticks_are_step_aligned_and_stop_inside_the_range() {
        // Step 5 over [0, 16]: 20 would fall outside, so the axis ends at 15.
        assert_eq!(axis_ticks(0.0, 16.0, 4), vec![0.0, 5.0, 10.0, 15.0]);
        // Aligned to multiples of the step, not to `lo`: the first is 10, not 7.
        assert_eq!(axis_ticks(7.0, 25.0, 4), vec![10.0, 15.0, 20.0, 25.0]);
        // A degenerate range draws nothing rather than looping or dividing by zero.
        assert!(axis_ticks(5.0, 5.0, 4).is_empty());
        assert!(axis_ticks(10.0, 0.0, 4).is_empty());
        assert!(axis_ticks(0.0, f64::NAN, 4).is_empty());
    }

    #[test]
    fn axis_decimals_prints_a_step_at_its_own_precision() {
        assert_eq!(axis_decimals(5.0), 0, "whole-second steps need no decimals");
        assert_eq!(axis_decimals(1.0), 0);
        assert_eq!(axis_decimals(0.5), 1);
        assert_eq!(axis_decimals(0.05), 2);
        // Capped: hundredths are the finest this axis ever prints.
        assert_eq!(axis_decimals(0.0005), 2);
        assert_eq!(axis_decimals(0.0), 0, "a degenerate step must not panic");
    }

    #[test]
    fn axis_scale_pairs_the_step_with_the_precision_its_labels_need() {
        // Sign-insensitive: only the span's magnitude sets the step.
        assert_eq!(axis_scale(16.0, 0.0, 4), (5.0, 0));
        // A zero span falls back to `EPS` rather than producing a NaN step, whose
        // precision would then print as no decimals at all.
        assert_eq!(axis_scale(3.0, 3.0, 4).1, 2);
    }

    #[test]
    fn value_to_row_flips_the_y_axis_and_value_to_index_does_not() {
        // The y flip's two endpoints, the pair an off-by-one in `canvas_h - 1 - dr`
        // moves: `hi` is screen row 0, `lo` the last. The unflipped index inverts.
        assert_eq!(value_to_row(10.0, 0.0, 10.0, 10), 0);
        assert_eq!(value_to_row(0.0, 0.0, 10.0, 10), 9);
        assert_eq!(
            value_to_index(10.0, 0.0, 10.0, 10),
            9,
            "the last column, not 10"
        );
        assert_eq!(value_to_index(0.0, 0.0, 10.0, 10), 0);
        // Halfway up a ten-row canvas is five rows below the top, so row 4.
        assert_eq!(value_to_row(5.0, 0.0, 10.0, 10), 4);
        assert_eq!(value_to_index(5.0, 0.0, 10.0, 10), 5);
        assert_eq!(value_to_row(99.0, 0.0, 10.0, 10), 0);
        assert_eq!(value_to_row(-99.0, 0.0, 10.0, 10), 9);
        assert_eq!(value_to_index(99.0, 0.0, 10.0, 10), 9);
        // A degenerate span or an empty canvas answers index 0 rather than panicking;
        // the flip reads that as the bottom row.
        assert_eq!(value_to_index(5.0, 10.0, 10.0, 10), 0);
        assert_eq!(value_to_row(5.0, 10.0, 10.0, 10), 9);
        assert_eq!(value_to_index(5.0, 0.0, 10.0, 0), 0);
        assert_eq!(value_to_row(5.0, 0.0, 10.0, 0), 0);
    }

    // ---- scaled_u64: NaN must not read as convergence ----

    #[test]
    fn scaled_u64_carries_the_last_finite_value_across_a_nan() {
        // Batch 0's own delta and any failed batch's metrics are NaN — this must not
        // draw as the scale's `0`, which looks identical to "the curve stopped".
        let scaled = scaled_u64(&[10.0, f64::NAN, 10.0]);
        assert_eq!(
            scaled[1], scaled[0],
            "a NaN sample must repeat the prior finite value, not read as 0: {scaled:?}"
        );
        assert_eq!(scaled[2], scaled[0]);

        // A *leading* NaN run has nothing to carry forward — the one case where `0`
        // means "nothing to show yet" rather than "converged".
        let scaled = scaled_u64(&[f64::NAN, f64::NAN, 10.0]);
        assert_eq!(scaled[0], 0, "{scaled:?}");
        assert_eq!(scaled[1], 0, "{scaled:?}");
        assert!(
            scaled[2] > 0,
            "the first real sample must still scale normally: {scaled:?}"
        );
    }

    // ---- status line: per-tab bindings, the count only when pending, `? keys` pinned
    // right, and a degrade that never overflows `area.width`. ---------------------

    fn long_prefix_app() -> App {
        let mut app = App::new(10);
        for c in "123456789".chars() {
            press(&mut app, c);
        }
        app
    }

    fn status_line(app: &mut App, width: u16) -> String {
        render(app, width, 3)
            .lines()
            .last()
            .expect("draw() always renders a status row")
            .to_string()
    }

    fn line_text(line: &Line) -> String {
        line.spans.iter().map(|s| s.content.as_ref()).collect()
    }

    /// Inversion, not a color, is what makes a pending count unmistakable — and the
    /// glyph-only `render` harness cannot see style, so this reads the buffer.
    #[test]
    fn status_line_pending_count_is_reversed() {
        let mut app = App::new(10);
        press(&mut app, '4');
        assert_eq!(app.pending_count(), Some(4));

        let buf = draw_to_buffer(&mut app, 60, 3);
        let status_row = buf.area.height - 1;
        let found = (0..buf.area.width).find_map(|x| {
            let cell = &buf[(x, status_row)];
            (cell.symbol() == "4").then_some(cell.modifier)
        });
        assert_eq!(
            found,
            Some(Modifier::REVERSED),
            "the pending count must render REVERSED so it is unmistakable"
        );
    }

    #[test]
    fn status_line_hides_fit_only_bindings_on_other_tabs() {
        let mut app = App::new(10);
        goto_tab(&mut app, Tab::Convergence);
        let status = status_line(&mut app, 200);
        for word in ["frame", "layer", "dp"] {
            assert!(
                !status.contains(word),
                "Convergence must not advertise Fit-only bindings that do \
                 nothing there: {status:?}"
            );
        }
        assert!(status.contains("next"), "{status:?}");
        assert!(status.contains("run"), "{status:?}");
    }

    #[test]
    fn status_line_shows_the_keys_hint_whenever_there_is_room_for_it() {
        for width in [40, 60, 90, 120, 200] {
            let status = status_line(&mut long_prefix_app(), width);
            assert!(
                status.contains("keys"),
                "`? keys` is pinned to its own fixed-width column and should \
                 survive at width {width}: {status:?}"
            );
        }
    }

    // `fit_status_hints` directly, not through `render`, so the degrade widths do not
    // also depend on `draw_status_line`'s columns.

    /// The four stages, each at the exact width that produces it and one column past
    /// the width above it: full text, keys only, global keys only, nothing.
    #[test]
    fn fit_status_hints_degrades_one_stage_per_column_it_runs_out_of() {
        let tab_local = tab_keys(Tab::Fit);
        let full: Vec<Binding> = tab_local.iter().chain(GLOBAL_KEYS).copied().collect();
        let width_of =
            |bindings: &[Binding], actions| Line::from(binding_spans(bindings, actions)).width();
        let full_w = width_of(&full, true);
        let keys_w = width_of(&full, false);
        let global_w = width_of(GLOBAL_KEYS, false);
        assert!(
            full_w > keys_w && keys_w > global_w && global_w > 0,
            "each stage must be strictly narrower than the one above it or the \
             ladder below asserts nothing: {full_w} / {keys_w} / {global_w}"
        );

        // Stage 1, at exactly the width the full line needs.
        let text = line_text(&fit_status_hints(tab_local, full_w));
        for b in &full {
            let pair = format!("{} {}", b.keys, b.hint);
            assert!(text.contains(&pair), "missing {pair:?}: {text:?}");
        }

        // Stage 2, one column short of that: the action words go, every key stays.
        let text = line_text(&fit_status_hints(tab_local, full_w - 1));
        for b in &full {
            assert!(
                !text.contains(b.hint),
                "{:?} must be gone: {text:?}",
                b.hint
            );
            assert!(text.contains(b.keys), "{:?} must survive: {text:?}", b.keys);
        }

        // Stage 3, one column short again: the tab-local group goes.
        let text = line_text(&fit_status_hints(tab_local, keys_w - 1));
        for b in tab_local {
            assert!(
                !text.contains(b.keys),
                "tab-local {:?} must be gone: {text:?}",
                b.keys
            );
        }
        for b in GLOBAL_KEYS {
            assert!(
                text.contains(b.keys),
                "global {:?} must survive: {text:?}",
                b.keys
            );
        }

        // Stage 4: below the global keys' width there is nothing left, zero included.
        for width in [global_w - 1, 0] {
            let text = line_text(&fit_status_hints(tab_local, width));
            assert!(text.is_empty(), "at width {width}: {text:?}");
        }
    }

    /// Reporting `-8.5 .. +9.5` as `±9.5` understates one side by two ppm and
    /// `±8.5` overstates the other — the reason `mz_ppm` is a pair at all.
    #[test]
    fn fmt_interval_collapses_only_the_symmetric_case() {
        assert_eq!(fmt_interval((-3.0, 3.0)), "±3.0");
        assert_eq!(fmt_interval((-8.5, 9.5)), "-8.5 .. +9.5");
    }

    /// `now` must report what the right edge draws — the last *finite* sample under
    /// `scaled_u64`'s hold. Reporting `—` beside a held bar contradicts the plot.
    #[test]
    fn spark_title_reports_the_value_the_right_edge_actually_draws() {
        let title = spark_title("wrmse", &[2.0, 0.5, f64::NAN], true);
        assert!(title.contains("peak 2"), "{title}");
        assert!(title.contains("now 0.5000"), "{title}");

        // Batch 0's `max_delta` is NaN, so the plot is flat at the floor, where
        // `peak 0  now 0` would read as "converged" rather than "not measured".
        let title = spark_title("max_delta", &[f64::NAN], true);
        assert_eq!(title, "max_delta  peak —  now —  (NaN holds)");
    }

    /// Every tab, every mark layer, both DP-pane states, a scrubbed frame, an
    /// empty app and the keys overlay, drawn at every size the layout branches
    /// differently at.
    ///
    /// A panic here is not cosmetic: `[profile.release]` sets `panic = "abort"`, so a
    /// user who shrinks their terminal during a pause kills the search, and terminal
    /// size is outside the program's control.
    ///
    /// The sizes are the layout's own thresholds and their neighbours, not a product,
    /// plus a tiny-by-tiny block for the `width < 3 || block_h < 3` disjunction with
    /// both sides true.
    #[test]
    fn every_screen_draws_at_any_terminal_size() {
        // A sweep artifact, not a rendering requirement: ratatui memoizes
        // `Layout::split` in a 500-entry thread-local LRU, which thrashes once one
        // test walks thousands of distinct areas.
        Layout::init_cache(NonZeroUsize::new(20_000).expect("nonzero"));

        let thresholds: Vec<(u16, u16)> = SWEEP_WIDTHS
            .map(|w| (w, 12))
            .into_iter()
            .chain(SWEEP_HEIGHTS.map(|h| (60, h)))
            .chain(SWEEP_TINY)
            .chain(SWEEP_CROSSINGS)
            .collect();

        // One app per combination, redrawn at every size: drawing changes none of the
        // state set up here, and rebuilding the fixture per size dominates the cost.
        let draw_all = |app: &mut App, sizes: &[(u16, u16)]| {
            for (w, h) in sizes.iter().copied() {
                draw_to_buffer(app, w, h);
            }
        };
        // Everything that paints the heatmap by hand: every threshold, plus the
        // oversized and lopsided areas its raw indexing could walk off.
        let sweep_painted = |app: &mut App| {
            draw_all(app, &thresholds);
            draw_all(app, &SWEEP_OVERSIZED);
        };

        // The `?` overlay is swept once rather than per tab: it `Clear`s and redraws
        // the same panel over the whole `frame.area()` regardless of what is beneath.
        let mut keys_open = fixture_app_with_metrics();
        press(&mut keys_open, '?');
        sweep_painted(&mut keys_open);

        for dp in [false, true] {
            let mut app = fixture_app_with_metrics();
            if dp {
                press(&mut app, 'd');
            }
            sweep_painted(&mut app);
        }

        // A scrubbed frame is `draw_fit_tab`'s `Length(1), Min(0)` split — the one
        // path that deliberately hands the heatmap zero rows.
        let mut scrubbed = fixture_app_with_ridge();
        scrubbed.set_frame_summary(RETAINED_5);
        scrubbed.set_scrub_recording(2, 17, fixture_recording(8));
        sweep_painted(&mut scrubbed);

        // A mark layer changes the glyphs, never the layout, so it needs the size
        // sweep once rather than once per layout combination above — and only on the
        // Fit tab, the only tab that reads `app.layer()`.
        for layer in Layer::ALL.into_iter().filter(|l| *l != Layer::None) {
            let mut app = fixture_app_with_metrics();
            goto_layer(&mut app, layer);
            sweep_painted(&mut app);
        }

        // The other two tabs, populated and empty, plus the empty Fit tab. Widget
        // layout clips and needs no size guards, so these take the short list.
        for tab in [Tab::Convergence, Tab::Tolerances] {
            let mut app = fixture_app_with_metrics();
            goto_tab(&mut app, tab);
            draw_all(&mut app, &SWEEP_SMOKE);
        }
        for tab in Tab::ALL {
            let mut app = App::new(10); // no frames, no metrics, no fit
            goto_tab(&mut app, tab);
            draw_all(&mut app, &SWEEP_SMOKE);
        }
    }
}
