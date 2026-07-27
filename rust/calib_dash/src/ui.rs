//! Rendering for the three dashboard tabs (`draw`) and the heatmap
//! downsampler that feeds the Fit tab (`heatmap_cells`).
//!
//! The Fit tab is the one built with the most care: `FitRecording::path_indices`
//! is `prefix ++ DP chain ++ suffix` (see `recording.rs`), and a calibration
//! that misbehaves at the edges of the gradient is often the greedily
//! attached tail's doing rather than the DP's. So the path overlay paints
//! `path[dp_range]` with one glyph and the two tails with another — that
//! split is the first thing a user reading this tab needs to see.
//!
//! Exactly one mark layer (`crate::Layer`) is drawn on the heatmap at a
//! time, cycled by `m`/`M`: `none`, `path`, `curve`, `ridge`, `suppressed`.
//! Each layer is internally overlap-free (the DP chain and greedy tails
//! partition the path by construction; the other three layers only ever
//! emit one mark kind), so there is no cross-layer priority to resolve —
//! unlike the composited overlay this replaced, a layer is a single,
//! isolated view.
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
    Layer,
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

    if app.show_keys() {
        draw_keys_overlay(frame, area);
    }
}

fn draw_tab_bar(frame: &mut Frame, area: Rect, app: &App) {
    let titles: Vec<&str> = Tab::ALL.iter().map(|t| t.title()).collect();
    let selected = Tab::ALL.iter().position(|t| *t == app.tab()).unwrap_or(0);
    // `REVERSED`, not a hue: color is reserved for "what a mark is" on the
    // Fit tab (Item 3/4), and `REVERSED` is reserved for "mode" everywhere
    // else — which tab is selected, whether the Fit tab is showing a
    // scrubbed frame, whether a count is pending. Before this, `Yellow` meant
    // "selected tab" here *and* "not live data" on the scrub banner *and*
    // "greedy tail" on the heatmap all at once.
    let tabs = Tabs::new(titles)
        .select(selected)
        .highlight_style(Style::default().add_modifier(Modifier::REVERSED | Modifier::BOLD));
    frame.render_widget(tabs, area);
}

/// One `key`/`action` hint, rendered as the key in `BOLD` and the action in
/// `DarkGray` — weight, not `key:action` punctuation, is what makes the bound
/// letter scannable (Item 2). `action` empty renders the key alone, which is
/// how the overflow degrade below drops action words while keeping keys.
fn hint_spans(pairs: &[(&str, &str)]) -> Vec<Span<'static>> {
    let mut spans = Vec::new();
    for (i, (key, action)) in pairs.iter().enumerate() {
        if i > 0 {
            spans.push(Span::styled(
                " \u{b7} ",
                Style::default().fg(Color::DarkGray),
            ));
        }
        spans.push(Span::styled(
            (*key).to_string(),
            Style::default().add_modifier(Modifier::BOLD),
        ));
        if !action.is_empty() {
            spans.push(Span::raw(" "));
            spans.push(Span::styled(
                (*action).to_string(),
                Style::default().fg(Color::DarkGray),
            ));
        }
    }
    spans
}

fn spans_width(spans: &[Span]) -> usize {
    spans.iter().map(|s| s.content.chars().count()).sum()
}

/// Picks the most detailed hint line (tab-local bindings plus the global
/// ones) that still fits `width`, degrading in the stages Item 2 calls for:
/// full text with action words, then keys only, then only the global
/// bindings' keys, then nothing at all. `? keys` is never part of this line —
/// it lives in its own pinned-right column (`draw_status_line`) and so
/// survives regardless of how far this one degrades.
fn fit_status_hints(
    tab_local: &[(&str, &str)],
    global: &[(&str, &str)],
    width: usize,
) -> Line<'static> {
    let full: Vec<(&str, &str)> = tab_local.iter().chain(global.iter()).copied().collect();
    let with_actions = hint_spans(&full);
    if spans_width(&with_actions) <= width {
        return Line::from(with_actions);
    }
    let keys_only: Vec<(&str, &str)> = full.iter().map(|(k, _)| (*k, "")).collect();
    let keys_line = hint_spans(&keys_only);
    if spans_width(&keys_line) <= width {
        return Line::from(keys_line);
    }
    let global_keys_only: Vec<(&str, &str)> = global.iter().map(|(k, _)| (*k, "")).collect();
    let global_line = hint_spans(&global_keys_only);
    if spans_width(&global_line) <= width {
        return Line::from(global_line);
    }
    Line::default()
}

/// The status line: batch/stage/pending-count on the left, key hints in the
/// middle (degrading as `fit_status_hints` describes when the terminal is
/// narrow), `? keys` pinned to the right. Three `Paragraph`s rather than one
/// `Line`, because a `Line` has a single alignment and `? keys` must stay
/// right-aligned regardless of how much the middle grows or shrinks.
fn draw_status_line(frame: &mut Frame, area: Rect, app: &App) {
    if area.width == 0 || area.height == 0 {
        return;
    }
    let mut state_spans = vec![Span::raw(format!(
        " b{} {} ",
        app.batch(),
        app.stage().label()
    ))];
    // The pending vim-style count only takes a column when there is one to
    // show — unlike the line this replaced, which spent six columns
    // permanently reporting its absence (`cnt:-`). `REVERSED` marks it as a
    // *mode* (a keystroke is being accumulated), the same signal the tab bar
    // and the scrub banner use, not a color — see `draw_tab_bar`.
    if let Some(n) = app.pending_count() {
        state_spans.push(Span::styled(
            n.to_string(),
            Style::default().add_modifier(Modifier::REVERSED),
        ));
        state_spans.push(Span::raw(" "));
    }
    let state_w = (spans_width(&state_spans) as u16).min(area.width);

    // ---- per-tab binding sets (Item 2) ----------------------------------
    // Convergence and Tolerances only ever respond to `n`/`r` — advertising
    // the Fit tab's stage/frame/overlay/dp keys on them (as the line this
    // replaced did) is four bindings that do nothing there.
    let tab_local: &[(&str, &str)] = match app.tab() {
        Tab::Fit => &[
            ("n", "next"),
            ("r", "run"),
            ("[ ]", "stage"),
            ("< >", "frame"),
            ("m M", "layer"),
            ("d", "dp"),
        ],
        Tab::Convergence | Tab::Tolerances => &[("n", "next"), ("r", "run")],
    };
    let global: &[(&str, &str)] = &[("h l", "tab"), ("q", "detach"), ("^C", "abort")];

    let right_w = 7u16.min(area.width);
    let mid_w = area.width.saturating_sub(state_w).saturating_sub(right_w) as usize;
    let mid_line = fit_status_hints(tab_local, global, mid_w);

    let cols = Layout::horizontal([
        Constraint::Length(state_w),
        Constraint::Min(0),
        Constraint::Length(right_w),
    ])
    .split(area);

    frame.render_widget(Paragraph::new(Line::from(state_spans)), cols[0]);
    frame.render_widget(Paragraph::new(mid_line), cols[1]);
    frame.render_widget(
        Paragraph::new(Line::from(hint_spans(&[("?", "keys")]))),
        cols[2],
    );
}

/// The `?` overlay: the full key map, grouped by heading, drawn over
/// whatever the tab underneath was showing. `Clear` first so the overlay is
/// actually opaque — without it, the heatmap/table cells behind would show
/// through anywhere the overlay's own background does not touch. Dismissed
/// on any key (`App::handle_key`'s first check), so nothing here needs to
/// render a "press any key" footer.
fn draw_keys_overlay(frame: &mut Frame, area: Rect) {
    let w = 52.min(area.width);
    let h = 16.min(area.height);
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

    let heading = |s: &'static str| Line::styled(s, Style::default().add_modifier(Modifier::BOLD));
    let lines = vec![
        heading("Every tab"),
        Line::raw("  h l        switch tab"),
        Line::raw("  q          detach"),
        Line::raw("  ^C         abort"),
        Line::raw("  ?          this screen"),
        Line::raw(""),
        heading("Fit tab"),
        Line::raw("  n          next batch"),
        Line::raw("  r          run to end"),
        Line::raw("  [ ]        step pipeline stage"),
        Line::raw("  < >        scrub retained frames"),
        Line::raw("  m M        cycle mark layer"),
        Line::raw("  d          toggle DP pane"),
        Line::raw(""),
        heading("Convergence / Tolerances"),
        Line::raw("  n          next batch"),
        Line::raw("  r          run to end"),
    ];
    frame.render_widget(
        Paragraph::new(lines).block(Block::bordered().title(" Keys ")),
        popup,
    );
}

// ---------------------------------------------------------------------
// Fit tab
// ---------------------------------------------------------------------

fn draw_fit_tab(frame: &mut Frame, area: Rect, app: &App) {
    if area.width == 0 || area.height == 0 {
        return;
    }
    // The recording actually on screen: the scrubbed frame while `<`/`>`
    // have selected one, the live batch otherwise — see
    // `App::active_recording`.
    let rec = app.active_recording();
    if rec.geom().bins == 0 {
        frame.render_widget(Paragraph::new("No grid recorded yet."), area);
        return;
    }

    // When scrubbing, carve out one banner row so it is never ambiguous that
    // this is a replayed batch rather than the live one — the whole reason
    // `<`/`>` exist is to compare an earlier batch against the current one,
    // which only works if it's always obvious which is which. Skipped
    // entirely for the live view (the common case); otherwise the banner
    // always wins the one row it needs, even when that leaves the heatmap
    // zero rows (`draw_heatmap` already handles a zero-height area without
    // panicking) — a heatmap too short to show anything useful is a worse
    // outcome than an unlabeled one that a user might mistake for live.
    let (banner, area) = match app.scrub_frame() {
        Some(i) => {
            let rows = Layout::vertical([Constraint::Length(1), Constraint::Min(0)]).split(area);
            (Some((i, rows[0])), rows[1])
        }
        None => (None, area),
    };
    if let Some((i, banner_area)) = banner {
        let total = app.retained_frames().max(1);
        let batch_note = app
            .scrub_chunk()
            .map(|c| format!(", batch {c}"))
            .unwrap_or_default();
        let text = format!(
            " SCRUBBED — retained frame {}/{}{batch_note} (not live; `>` returns to now)",
            i + 1,
            total
        );
        // `REVERSED`, not a hue — this is a *mode* indicator ("you are not
        // looking at live data"), the same signal `draw_tab_bar`'s selected
        // tab and the status line's pending count use. Before Item 4, this
        // was `Yellow`, the same color the heatmap's greedy-tail mark now
        // also uses — one hue meaning three different things.
        frame.render_widget(
            Paragraph::new(text)
                .style(Style::default().add_modifier(Modifier::REVERSED | Modifier::BOLD)),
            banner_area,
        );
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
        // The DP pane needs a recording actually fit with `dp_nodes: true`,
        // which `rec` (`active_recording()`) is not guaranteed to be: the
        // live batch's own recording is deliberately fit with
        // `ObserveOpts::NONE` (see `CalibDash::refit_live`'s doc comment) to
        // keep that cost off the per-batch hot path. While scrubbed, `rec`
        // already qualifies — `refit_frame` always observes DP nodes — so
        // only the live case needs the on-demand recording
        // `CalibDash::sync_dp` maintains. Falling back to `rec` when that is
        // `None` (pane just switched on and not yet synced, or the batch was
        // too degenerate to refit) still renders a real recording rather
        // than nothing — `draw_dp_pane` already says so when its `dp()` is
        // empty, instead of leaving the pane blank.
        let dp_rec = if app.scrub_frame().is_some() {
            rec
        } else {
            app.live_dp_recording().unwrap_or(rec)
        };
        draw_dp_pane(frame, dp_area, dp_rec);
    }
}

/// What one half of a terminal cell is marked with, for whichever layer
/// (`crate::Layer`) is currently active. At most one layer ever paints into
/// `marks` per frame — see `paint_heatmap` — so unlike the priority ladder
/// this replaced, there is no cross-layer case to resolve: `Region` is the
/// only mark kind the curve, ridge and suppressed layers ever emit (each is
/// a single, uniform region of the grid), and `DpNode`/`Tail` are the only
/// two the path layer emits, mutually exclusive by construction
/// (`FitRecording::dp_range` partitions `path_indices` into exactly those two
/// groups — see `mark_path`).
///
/// **Invariant (superseding the old "every mark is identified by its glyph"
/// rule this module used to state):** a `Region` mark keeps the density
/// glyph — `░▒▓█▀▄`, half-occupancy intact — and is identified by
/// `Modifier::REVERSED`, not by replacing the glyph. Inverting rather than
/// overwriting is what lets a mark and the distribution it sits on both stay
/// visible at once: the shape underneath answers "how much weight is here",
/// and the inversion answers "is this cell marked", and neither question
/// erases the other's answer. `Modifier::REVERSED` also survives every
/// condition `mark_style`'s old comment cared about — a colorblind reader, a
/// monochrome screenshot, any terminal theme — for the same reason
/// `draw_tab_bar`/the scrub banner/the pending count already lean on it: it
/// is a swap of whatever foreground/background are already there, not a
/// specific hue that could fail to contrast or fail to register. Color
/// remains attached (`region_accent`) as pure redundant reinforcement, same
/// as before — a colorblind reader loses nothing by losing it.
///
/// A path-layer node mark is the one exception: `DpNode`/`Tail` still
/// *replace* the cell with `O`/`X` (Item 3's original marker set), because a
/// path node is a point estimate — one grid cell chosen by the DP or the
/// greedy walk — not a region being highlighted against its own density.
#[derive(Clone, Copy, PartialEq, Eq, Default)]
enum Mark {
    #[default]
    None,
    /// Curve, ridge, or suppressed layer — a region mark, rendered as the
    /// underlying density glyph (or `\u{b7}` in place of a space — see
    /// `marked_glyph`) reversed.
    Region,
    /// Path layer, DP's own chosen chain — replaces the cell with `O`.
    DpNode,
    /// Path layer, Pass 2's greedily attached tail — replaces the cell with
    /// `X`, the standard excluded/suspect-point convention: these are the
    /// nodes a greedy pass grafted on after the DP's own search ended.
    Tail,
}

/// Resolution order when both halves of one terminal cell carry different
/// marks (only possible within the path layer, where a DP node and a tail
/// node can round to the same display cell): a node mark always wins over a
/// region mark or no mark at all — `mark_glyph`'s doc comment above explains
/// why a path node replaces the cell rather than composing with it — and
/// `DpNode` wins over `Tail` between the two node kinds, so evidence (the
/// DP's own chain) still shows through a tail node sharing its cell. `Region`
/// only ever contends with `None` in practice (a layer never emits two
/// different mark kinds at once other than the path layer's two node kinds),
/// but is ranked above it for the same reason: a mark must never lose to "no
/// mark" just because it shares a cell with an unmarked half.
fn mark_rank(m: Mark) -> u8 {
    match m {
        Mark::None => 0,
        Mark::Region => 1,
        Mark::Tail => 2,
        Mark::DpNode => 3,
    }
}

/// `O`/`X` and their color/modifier — see `Mark::DpNode`/`Mark::Tail`'s doc
/// comments for why these replace the cell instead of inverting it.
/// `Modifier::REVERSED` is deliberately not used here — Item 4 reserves it
/// for *mode* elsewhere and `Mark::Region` for "this cell is marked" here, so
/// node glyphs lean on `BOLD` instead for emphasis, exactly as before this
/// pass.
fn node_glyph_and_style(m: Mark) -> (&'static str, Color, Modifier) {
    match m {
        Mark::DpNode => ("O", Color::Green, Modifier::BOLD),
        Mark::Tail => ("X", Color::Yellow, Modifier::BOLD),
        _ => (" ", Color::Reset, Modifier::empty()),
    }
}

/// The active layer's accent color for a `Region` mark — reinforcement only
/// (see `Mark`'s doc comment), so, as before, a named ANSI constant rather
/// than `Rgb` (see `heat_color`'s own doc comment for why) and never one of
/// indices 9-15 (`Light*`/`White`, invisible on a light background).
fn region_accent(layer: Layer) -> Color {
    match layer {
        Layer::Suppressed => Color::DarkGray,
        Layer::Curve => Color::Cyan,
        Layer::Ridge => Color::Magenta,
        // `Region` is never emitted while the active layer is `None`/`Path`
        // (see `paint_heatmap`), so this arm is unreachable in practice; it
        // exists only so the match is exhaustive rather than partial.
        Layer::None | Layer::Path => Color::Reset,
    }
}

/// `heat_glyph`, but a marked cell must never render as a reversed space.
/// `heat_glyph` returns a plain space for zero weight, and a *reversed*
/// space is a solid block — the single highest-density glyph this module
/// has — which would lie about the data wherever a region mark lands on a
/// zero-weight cell. That happens routinely: the curve is evaluated at every
/// display column regardless of whether any calibrant actually falls there,
/// and the ridge's `half_width` band frequently extends past the last real
/// observation. `\u{b7}` (a single dot — the least ink of any glyph this
/// module draws) reversed reads instead as "mark here, nothing under it",
/// which is the truth.
fn marked_glyph(heat: f32, max: f32) -> &'static str {
    match heat_glyph(heat, max) {
        " " => "\u{b7}",
        g => g,
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

/// The density field's color (Item 4). `░▒▓█` is already a luminance ramp —
/// the glyph alone carries the density — so this only reinforces it with at
/// most two gray steps, never `Rgb`: an `Rgb` terminal color is not remapped
/// by a terminal's theme, so on a light background the ramp inverts and the
/// hottest cells become the *least* visible, exactly backwards. The 0.5 cut
/// matches `heat_glyph`'s own mid-scale split (`▓`/`█` are `Gray`, `░`/`▒`
/// are the dimmer `DarkGray`).
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

/// Flips a display-row index so grid row 0 — the *lowest* observed RT —
/// lands at the *bottom* of the canvas rather than the top (Item 0: the
/// y-axis inversion fix).
///
/// `forward_map` alone would place grid row 0 at display row 0, which reads
/// naturally as an array index but is backwards for a plot: terminal rows
/// grow downward, so display row 0 is the terminal's *top* row, and a
/// monotonically increasing calibration (higher library RT -> higher
/// observed RT) would render as a *descending* line. This negates the index
/// within `0..disp_rows` instead.
///
/// Flipping the index also flips which half of the doubled half-block
/// resolution it falls into — `dr` and `disp_rows - 1 - dr` have opposite
/// parity whenever `disp_rows` is even (always true here: it is `h * 2`) —
/// so the `Upper`/`Lower` parity must be re-derived from the flipped index,
/// not carried over from the pre-flip one, or the `▀`/`▄` glyphs would paint
/// the wrong half of the terminal cell even though the row landed in the
/// right one.
fn flip_display_row(dr: usize, disp_rows: usize) -> (usize, Half) {
    let flipped = disp_rows.saturating_sub(1).saturating_sub(dr);
    let half = if flipped.is_multiple_of(2) {
        Half::Upper
    } else {
        Half::Lower
    };
    (flipped / 2, half)
}

/// Framed heatmap: a `Block` around the painted canvas (Item 1) stating the
/// axes in its title, a left gutter of y-tick labels, and an x-tick label
/// row below it, with `┤`/`┴` tick marks written into the border itself.
/// Falls back to painting straight into `area` with no frame at all when
/// there is not enough room for the frame to have any interior — see the
/// size guard below — rather than spending a tiny terminal's one or two
/// rows entirely on borders.
fn draw_heatmap(frame: &mut Frame, area: Rect, app: &App) {
    if area.width == 0 || area.height == 0 {
        return;
    }
    let rec = app.active_recording();
    let bins = rec.geom().bins;
    if bins == 0 {
        frame.render_widget(Paragraph::new("Grid has zero bins."), area);
        return;
    }
    let geom = rec.geom();
    let (x_lo, x_hi) = geom.x_range;
    let (y_lo, y_hi) = geom.y_range;

    // A framed canvas needs a `Block` with at least one interior row/column
    // to mean anything. `show_x_axis` (one row below the block for x-tick
    // labels) is only reserved when there is a row to spare for it.
    let show_x_axis = area.height as usize >= 4;
    let block_h = area.height as usize - usize::from(show_x_axis);
    if area.width < 3 || block_h < 3 {
        paint_heatmap(frame, area, rec, app);
        return;
    }
    let canvas_h = block_h - 2;
    let y_target = (canvas_h / 5).clamp(2, 8);
    let gutter_w = y_gutter_width(y_lo, y_hi, y_target);
    // The gutter (y-tick labels) is only reserved when there is width to
    // spare for it *and* at least a few columns of canvas left over —
    // below that, ticks would crowd the one thing this panel actually
    // exists to show, so this drops them instead.
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

    let title_left = " Fit \u{2014} observed RT (s) \u{2191} vs library RT (s) \u{2192} ";
    // `batch`/`stage` only — the active layer used to share this corner too,
    // but now has its own dedicated subtitle row (below) and no longer
    // competes with these for width. See `heatmap_title_right`'s doc comment
    // for why `stage` is still the one dropped here under a narrow DP pane.
    let title_right = heatmap_title_right(
        app.batch(),
        app.stage().label(),
        (block_area.width as usize).saturating_sub(title_left.chars().count()),
    );
    let block = Block::bordered()
        .title_top(Line::from(title_left).left_aligned())
        .title_top(Line::from(title_right).right_aligned());
    let inner = block.inner(block_area);
    frame.render_widget(block, block_area);

    // The active mark layer's dedicated subtitle row — see `fit_subtitle`'s
    // doc comment for why this can no longer share a line with anything
    // else. Styled `BOLD` (the file's existing convention for emphasis —
    // `heading` above, the node glyphs) specifically so it never blends into
    // one run of text with the unstyled `title_left`/`title_right` above it.
    // Only skipped when there is not even one row to spare for it beyond the
    // canvas's own minimum of one row — a heatmap with zero rows would defeat
    // the tab's entire purpose, so this is the one case the subtitle still
    // yields to something else, and it is a height, not a width, shortage.
    let show_subtitle = inner.height >= 2;
    let canvas = if show_subtitle {
        let rows = Layout::vertical([Constraint::Length(1), Constraint::Min(0)]).split(inner);
        let subtitle_area = rows[0];
        let subtitle = fit_subtitle(app.layer(), subtitle_area.width as usize);
        frame.render_widget(
            Paragraph::new(subtitle).style(Style::default().add_modifier(Modifier::BOLD)),
            subtitle_area,
        );
        rows[1]
    } else {
        inner
    };
    paint_heatmap(frame, canvas, rec, app);

    if show_y_axis && canvas.height > 0 {
        let y_step = nice_step((y_hi - y_lo).abs().max(1e-9), y_target);
        let y_decimals = axis_decimals(y_step);
        for v in axis_ticks(y_lo, y_hi, y_target) {
            let row = value_to_row(v, y_lo, y_hi, canvas.height as usize);
            let abs_y = canvas.y + row as u16;
            let label = format!(
                "{:>width$}",
                fmt_axis_value(v, y_decimals),
                width = gutter_area.width as usize
            );
            write_row_text(frame, gutter_area.x, abs_y, &label, Color::DarkGray);
            if let Some(cell) = frame.buffer_mut().cell_mut((block_area.x, abs_y)) {
                cell.set_char('\u{2524}'); // ┤
            }
        }
    }

    if show_x_axis && inner.width > 0 {
        let x_target = (inner.width as usize / 20).clamp(2, 10);
        let x_step = nice_step((x_hi - x_lo).abs().max(1e-9), x_target);
        let x_decimals = axis_decimals(x_step);
        let mut label_row = vec![' '; x_label_area.width as usize];
        for v in axis_ticks(x_lo, x_hi, x_target) {
            let col = value_to_col(v, x_lo, x_hi, inner.width as usize);
            let block_col = 1 + col; // offset for the block's left border
            place_centered(&mut label_row, block_col, &fmt_axis_value(v, x_decimals));
            if let Some(cell) = frame.buffer_mut().cell_mut((
                block_area.x + block_col as u16,
                block_area.y + block_area.height - 1,
            )) {
                cell.set_char('\u{2534}'); // ┴
            }
        }
        let text: String = label_row.into_iter().collect();
        frame.render_widget(
            Paragraph::new(text).style(Style::default().fg(Color::DarkGray)),
            x_label_area,
        );
    }
}

/// Picks the most detailed right-side block title that still fits `avail`
/// columns, degrading by dropping the pipeline stage label: `Stage` is
/// already shown on the global status line (`draw_status_line`) for every
/// tab, so losing it here under a narrow DP pane (that pane's grid area is
/// only 65% of the terminal, further reduced by the y-axis gutter — see
/// `draw_fit_tab`/`y_gutter_width`) costs nothing not already visible
/// elsewhere.
///
/// The active layer used to be a third field crammed into this same corner,
/// dropped first of the three when space ran out — exactly backwards, since
/// it has no other home on screen at all. It now has its own dedicated
/// subtitle row (`fit_subtitle`, rendered just below this title) that never
/// competes with `batch`/`stage` for width, so there is nothing left here to
/// rank it against.
///
/// This exists because `ratatui`'s `Block` draws each of its titles
/// independently into the border row with no collision detection — two
/// titles whose combined width exceeds the block's simply overlap and
/// garble each other rather than one yielding space to the other.
fn heatmap_title_right(batch: u32, stage: &str, avail: usize) -> String {
    let full = format!(" b{batch} \u{b7} {stage} ");
    if full.chars().count() <= avail {
        return full;
    }
    format!(" b{batch} ")
}

/// One phrase glossing what an active layer's marks mean, appended to its
/// name in `fit_subtitle`. Needed because a `Region` mark (`Mark`'s doc
/// comment) is identified purely by `Modifier::REVERSED` now, not by a
/// distinct glyph — there is no glyph-coded legend left to read the marks
/// against, so the subtitle is the only place left to say what an inverted
/// cell on the active layer actually means. Lives beside `heatmap_title_right`/
/// `fit_subtitle` rather than on `Layer` itself (`app.rs`): this is about how
/// `ui.rs` renders the layer, not a fact about the layer.
fn layer_gloss(layer: Layer) -> &'static str {
    match layer {
        Layer::None => "density only",
        Layer::Path => "O chosen, X greedy tail",
        Layer::Curve => "fitted calibration",
        Layer::Ridge => "tolerance band",
        Layer::Suppressed => "discarded by non-max suppression",
    }
}

/// The Fit heatmap's dedicated subtitle: which mark layer is active, spelled
/// out as a label ("Showing: ridge — tolerance band") rather than a key hint
/// ("L: path") — a reader who does not know the active layer cannot
/// interpret a single inverted cell on this tab at all, so this states it in
/// full rather than abbreviating.
///
/// Degrades in the same fits-or-drops style as `heatmap_title_right`, but
/// with the gloss as the one thing shed first and the layer name itself
/// truncated (never wrapped) only as the last resort, character by
/// character, so this can never wrap or panic — including at `avail == 0`,
/// where it renders as an empty line rather than nothing at all crashing.
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

/// Paints the half-block heatmap itself — density field plus the
/// the active mark layer (`app.layer()`) — into `area`: the framed canvas's
/// inner rect in the common case, or the whole Fit-tab body on a terminal
/// too small to frame at all (`draw_heatmap`'s size guard above).
fn paint_heatmap(frame: &mut Frame, area: Rect, rec: &FitRecording, app: &App) {
    if area.width == 0 || area.height == 0 {
        return;
    }
    let bins = rec.geom().bins;
    if bins == 0 {
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
    // `half` 0 for upper, 1 for lower — so a mark on one grid row never
    // clobbers what the *other* grid row sharing that character wanted to
    // show. Exactly one of these runs per draw — the whole point of Layer
    // cycling instead of independently toggled overlays — so `marks` only
    // ever carries one layer's marks at a time.
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
            let idx = (ty * dims.disp_cols + tx) * 2;
            let upper_heat = cells.get(idx).copied().unwrap_or(0.0);
            let lower_heat = cells.get(idx + 1).copied().unwrap_or(0.0);
            let cell = ty * w + tx;
            let upper_mark = marks[cell * 2];
            let lower_mark = marks[cell * 2 + 1];
            let (symbol, color, modifier) =
                compose_cell(upper_mark, lower_mark, upper_heat, lower_heat, max_w, layer);
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

/// A "nice" axis step (1/2/5 x 10^k) close to `span / target` — e.g. seconds
/// on a 37s canvas at a target of 8 ticks rounds to a step of 5s, not an
/// ugly, unreadable 4.625.
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

/// Tick values in `[lo, hi]` spaced by `nice_step`, starting at the first
/// step-aligned value `>= lo`. Empty for a degenerate (zero-width or
/// non-finite) range or a non-finite/non-positive step — the caller then
/// simply draws no ticks rather than dividing by zero or looping on a `NaN`
/// step.
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
    // `target` already bounds how many ticks make sense on this canvas; a
    // small multiple covers float rounding at the high end without ever
    // looping unboundedly on a pathological step.
    let max_ticks = target.saturating_mul(3).clamp(4, 64);
    let mut out = Vec::with_capacity(max_ticks);
    let mut v = first;
    while v <= hi + step * 1e-6 && out.len() < max_ticks {
        out.push(v);
        v += step;
    }
    out
}

/// Decimal places to format an axis tick with, derived from the step size —
/// a step of `5.0` needs none, a step of `0.25` needs two. Capped at 2: RT
/// seconds are never usefully shown finer than hundredths here.
fn axis_decimals(step: f64) -> usize {
    if !step.is_finite() || step <= 0.0 {
        return 0;
    }
    (-step.log10().floor()).clamp(0.0, 2.0) as usize
}

fn fmt_axis_value(v: f64, decimals: usize) -> String {
    format!("{v:.decimals$}")
}

/// Maps a data value to a 0-indexed row inside a `canvas_h`-row canvas, with
/// the y axis flipped the same way `flip_display_row` flips the grid — the
/// largest value lands at row 0 (screen top), matching a normal plot rather
/// than terminal rows growing downward from grid row 0.
fn value_to_row(v: f64, lo: f64, hi: f64, canvas_h: usize) -> usize {
    if canvas_h == 0 {
        return 0;
    }
    let span = hi - lo;
    if !span.is_finite() || span <= 0.0 {
        return 0;
    }
    let t = ((v - lo) / span).clamp(0.0, 1.0);
    let dr = ((t * canvas_h as f64) as usize).min(canvas_h - 1);
    canvas_h - 1 - dr
}

/// The x-axis counterpart of `value_to_row` — no flip; the library-RT axis
/// runs left to right in screen order already.
fn value_to_col(v: f64, lo: f64, hi: f64, canvas_w: usize) -> usize {
    if canvas_w == 0 {
        return 0;
    }
    let span = hi - lo;
    if !span.is_finite() || span <= 0.0 {
        return 0;
    }
    let t = ((v - lo) / span).clamp(0.0, 1.0);
    ((t * canvas_w as f64) as usize).min(canvas_w - 1)
}

/// Left-gutter width for the y-axis tick labels: enough for the wider of
/// `lo`/`hi` formatted at the tick step's own precision, plus one column of
/// breathing room before the axis border, clamped to a sane range so an
/// extreme value never eats the whole canvas.
fn y_gutter_width(lo: f64, hi: f64, target: usize) -> usize {
    let step = nice_step((hi - lo).abs().max(1e-9), target);
    let decimals = axis_decimals(step);
    let lo_len = fmt_axis_value(lo, decimals).chars().count();
    let hi_len = fmt_axis_value(hi, decimals).chars().count();
    (lo_len.max(hi_len) + 1).clamp(4, 8)
}

/// Writes `text` starting at terminal column `x`, row `y`, one buffer cell
/// per `char` — used for the y-axis gutter labels, which are a handful of
/// characters on an otherwise-untouched row rather than a whole `Paragraph`.
fn write_row_text(frame: &mut Frame, x: u16, y: u16, text: &str, color: Color) {
    for (i, ch) in text.chars().enumerate() {
        let Some(cell) = frame.buffer_mut().cell_mut((x + i as u16, y)) else {
            break;
        };
        cell.set_char(ch);
        cell.set_fg(color);
    }
}

/// Centers `label` on column `col` within `row` (a fixed-width buffer of
/// chars), clipping at either edge rather than wrapping — an x-tick label
/// near the canvas edge must not spill onto the next line, panic, or wrap.
fn place_centered(row: &mut [char], col: usize, label: &str) {
    let width = row.len();
    if width == 0 {
        return;
    }
    let chars: Vec<char> = label.chars().collect();
    let start = col.saturating_sub(chars.len() / 2);
    for (i, ch) in chars.into_iter().enumerate() {
        let pos = start + i;
        if pos < width {
            row[pos] = ch;
        }
    }
}

/// Picks one terminal cell's glyph and color from its two independent
/// half-rows.
///
/// "Occupied" means either half has a nonzero weight or a mark.  Occupancy
/// shape is the primary, `TestBackend`-visible signal (`▀` upper only, `▄`
/// lower only, `█`/a density glyph when both, ` ` when neither) — this is
/// what makes `heatmap_cells`'s two-half-rows-per-line resolution actually
/// show up in a snapshot instead of being silently collapsed to one glyph
/// per cell.
///
/// When only one half is occupied, a mark on that half decides the glyph via
/// `compose_marked` (below); with no mark, that half's own heat intensity
/// picks the color and the shape is the plain `▀`/`▄` half-block (there is
/// no standard partial-block character for "the upper half, at quarter
/// density").
///
/// When *both* halves are occupied, `heatmap_cells` collapsing them to one
/// glyph is unavoidable — one character cannot show two independent marks at
/// once — so `winning_mark` (below) picks whichever of the two halves' marks
/// outranks the other, and the cell falls back to a density glyph over the
/// combined (max) heat only when *neither* half carries a mark at all.
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
            let winner = winning_mark(upper_mark, lower_mark);
            let heat = upper_heat.max(lower_heat);
            compose_marked(winner, heat, max, layer, heat_glyph(heat, max))
        }
    }
}

/// Resolves one already-occupied half (or, from `compose_cell`'s
/// both-occupied case, an already-picked winning mark for the whole cell)
/// into a glyph/color/modifier. `none_glyph` is what to draw when `mark` is
/// `Mark::None` — the asymmetric `▀`/`▄` half-block for a single occupied
/// half, or the combined density glyph `compose_cell` passes when both
/// halves are occupied but neither carries a mark.
fn compose_marked(
    mark: Mark,
    heat: f32,
    max: f32,
    layer: Layer,
    none_glyph: &'static str,
) -> (&'static str, Color, Modifier) {
    match mark {
        Mark::DpNode | Mark::Tail => node_glyph_and_style(mark),
        Mark::Region => (
            marked_glyph(heat, max),
            region_accent(layer),
            Modifier::REVERSED,
        ),
        Mark::None => (none_glyph, heat_color(heat, max), Modifier::empty()),
    }
}

/// Which of two marks sharing one terminal cell wins the whole cell — see
/// `mark_rank`'s doc comment for the ordering rationale.
fn winning_mark(a: Mark, b: Mark) -> Mark {
    if mark_rank(b) > mark_rank(a) { b } else { a }
}

fn raise_mark(marks: &mut [Mark], dims: Dims, tx: usize, ty: usize, half: Half, level: Mark) {
    if tx >= dims.w || ty >= dims.h {
        return;
    }
    let slot_idx = (ty * dims.w + tx) * 2 + if half == Half::Upper { 0 } else { 1 };
    if let Some(slot) = marks.get_mut(slot_idx)
        && mark_rank(level) > mark_rank(*slot)
    {
        *slot = level;
    }
}

fn mark_grid_indices(marks: &mut [Mark], dims: Dims, indices: &[usize], level: Mark) {
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

/// Marks `path[..dp_range.start]` and `path[dp_range.end..]` (Pass 2's
/// greedily attached tails, `Mark::Tail`) and `path[dp_range]` (the DP's own
/// chosen chain, `Mark::DpNode`) — the split the brief calls out as the
/// answer to the first question a user asks about a misbehaving edge of the
/// fit. The two are disjoint by construction (`dp_range` partitions
/// `path_indices`), so unlike the other three layers this one draws two
/// distinct marks and relies on `mark_rank` to arbitrate the rare terminal
/// cell where downsampling rounds one of each into the same slot.
fn mark_path(marks: &mut [Mark], dims: Dims, rec: &FitRecording) {
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

    mark_grid_indices(marks, dims, &path[..start], Mark::Tail);
    mark_grid_indices(marks, dims, &path[start..end], Mark::DpNode);
    mark_grid_indices(marks, dims, &path[end..], Mark::Tail);
}

/// Marks the fitted curve at every display column (not just at path nodes),
/// so it renders as a continuous line rather than sparse dots.
fn mark_curve(marks: &mut [Mark], dims: Dims, rec: &FitRecording) {
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
        // Item 0: re-derive the flipped row/half from `dr` via
        // `flip_display_row` rather than computing `Half` from `dr`'s own
        // parity — `grid_to_screen`'s doc comment explains why the parity
        // must come from the flipped index, not the pre-flip one.
        let (ty, half) = flip_display_row(dr, dims.disp_rows);
        raise_mark(marks, dims, tx, ty, half, Mark::Region);
    }
}

fn mark_ridge(marks: &mut [Mark], dims: Dims, rec: &FitRecording) {
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
            raise_mark(marks, dims, tx, ty, half, Mark::Region);
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
///
/// `dr` walks source row ranges in grid order (low to high), but is written
/// to the *flipped* display position — Item 0's y-axis fix: grid row 0 (the
/// lowest observed RT) must land at the bottom of the canvas, not the top,
/// the same flip `grid_to_screen` applies to the overlay marks so the two
/// stay in agreement.
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
            let flipped = disp_rows - 1 - dr;
            let cell_index = (flipped / 2) * area_w + dc;
            let half = flipped % 2;
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

/// Maps a grid `(row, col)` to `(terminal_row, terminal_col, half)`, flipping
/// the row so grid row 0 (the lowest observed RT) lands at the bottom of the
/// canvas — see `flip_display_row`'s doc comment for why this must be a
/// flip, not a swap, and why the `Half` parity has to be re-derived from the
/// flipped index.
fn grid_to_screen(row: usize, col: usize, dims: Dims) -> (usize, usize, Half) {
    let dr = forward_map(row, dims.bins, dims.disp_rows);
    let dc = forward_map(col, dims.bins, dims.disp_cols);
    let (ty, half) = flip_display_row(dr, dims.disp_rows);
    (ty, dc, half)
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

    // Titles for the four series that pass through `scaled_u64` say what a
    // flat run actually means now that a NaN sample holds the last finite
    // value instead of dropping to 0 — see that function's doc comment for
    // why: a flat run at batch 0, or across a failed batch, must not read
    // identically to "the curve stopped moving."
    let series: [(&str, Vec<u64>); 5] = [
        ("wrmse (NaN holds)", scaled_u64(&wrmse)),
        ("max_delta (NaN holds)", scaled_u64(&max_delta)),
        ("mean_delta (NaN holds)", scaled_u64(&mean_delta)),
        ("path_nodes", path_nodes),
        ("ridge_half_width (NaN holds)", scaled_u64(&ridge_hw)),
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
/// `u64` and so cannot mark a sample as "no data" distinctly from `0`.
/// NaN/infinite samples are a real possibility here — batch 0 always has a
/// NaN `max_delta`/`mean_delta` (no prior curve to compare against yet), and
/// so does any batch whose fit fails outright (`wrmse`, `ridge_half_width`
/// too) — and mapping every one of them to the scale's `0` would draw
/// identically to "the curve stopped moving," on the one tab whose job is
/// showing when it actually did.
///
/// Instead, each non-finite sample repeats the last finite value seen so far
/// (holds rather than drops), and only a run of non-finite samples with no
/// finite value before them yet — only possible at the very start of the
/// series — reports as `0`, which is genuinely "nothing to show yet," not
/// "converged." `draw_sparklines`'s titles say a flat run means this.
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
                last = ((v / max) * 1000.0).round().clamp(0.0, 1000.0) as u64;
            }
            last
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
            match app.tolerances() {
                Some(t) => {
                    lines.push(Line::raw(format!(
                        "m/z tolerance: ({:.1}, {:.1}) ppm",
                        t.mz_ppm.0, t.mz_ppm.1
                    )));
                    lines.push(Line::raw(format!(
                        "mobility tolerance: ({:.1}, {:.1}) %",
                        t.mobility_pct.0, t.mobility_pct.1
                    )));
                    // `rt_seconds` is symmetric by construction (both tuple
                    // elements are always equal — see `ToleranceSummary`'s
                    // doc comment), so this reports it as a single ± value
                    // rather than a `(lo, hi)` pair that would read as if the
                    // two sides could differ.
                    lines.push(Line::raw(format!(
                        "RT tolerance (symmetric): ±{:.1}s across {} calibrant(s)",
                        t.rt_seconds.0, t.n_calibrants
                    )));
                }
                None => {
                    lines.push(Line::raw(
                        "m/z and mobility distributions are Step B measurements outside this \
                         RT recording; they have not been wired through for this run.",
                    ));
                }
            }
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

    /// One snapshot per stop of the `m`/`M` cycle — the direct replacement
    /// for the old per-`Stage` snapshots. Unlike those (`Stage::Ridge` used
    /// to imply every overlay was on at once, via `sync_overlays_to_stage`),
    /// each of these renders in isolation: `Layer::Path` shows only the DP
    /// chain/tail, not also the suppressed mask, and so on.
    #[test]
    fn fit_tab_renders_each_layer() {
        for layer in Layer::ALL {
            let mut app = fixture_app_with_ridge();
            app.set_layer(layer);
            insta::assert_snapshot!(format!("fit_layer_{:?}", layer), render(&mut app, 100, 30));
        }
    }

    /// Pins the Fit tab's title-line home for the active layer's name (Layer's
    /// doc comment: `stage` and `layer` are independent axes now). Stepping
    /// the pipeline stage must not move the mark layer, and cycling the mark
    /// layer must not move the pipeline stage — a regression back toward the
    /// old `sync_overlays_to_stage` coupling would fail this.
    #[test]
    fn fit_tab_stage_and_layer_are_independent_axes() {
        let mut app = fixture_app_with_ridge();
        app.set_layer(Layer::Curve);
        app.set_stage(Stage::Ridge);
        assert_eq!(
            app.layer(),
            Layer::Curve,
            "changing the pipeline stage must not change the mark layer"
        );

        let out = render(&mut app, 100, 30);
        assert!(
            out.contains("Ridge") && out.contains("curve"),
            "the title states both axes independently:\n{out}"
        );
    }

    /// Guards the degrade order directly: the active layer used to share
    /// `heatmap_title_right`'s corner with `batch`/`stage` and was the first
    /// of the three dropped under a narrow DP pane — exactly backwards, since
    /// a reader who doesn't know the active layer cannot interpret a single
    /// inverted mark on this tab at all. It now has its own subtitle row
    /// (`fit_subtitle`) instead, spelled out in full with a gloss of what its
    /// marks mean, not abbreviated to a key hint. Checked across more than
    /// one layer so a future change to the degrade order (or to
    /// `fit_subtitle`/`layer_gloss` themselves) can't silently drop this for
    /// only some layers.
    #[test]
    fn fit_tab_subtitle_names_and_glosses_the_active_layer() {
        for (layer, gloss) in [
            (Layer::Ridge, "tolerance band"),
            (Layer::Suppressed, "discarded by non-max suppression"),
        ] {
            let mut app = fixture_app_with_ridge();
            app.set_layer(layer);
            let out = render(&mut app, 100, 30);
            assert!(
                out.contains("Showing:") && out.contains(layer.label()) && out.contains(gloss),
                "expected the {layer:?} layer's subtitle to name and gloss it:\n{out}"
            );
        }
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

    /// `App::set_scrub_recording` is what `CalibDash::sync_scrub` calls once
    /// `<`/`>` have moved `scrub_frame` — this pins that the Fit tab actually
    /// switches to that recording (a different, distinguishable grid from
    /// the live one) and draws the "not live" banner, rather than silently
    /// continuing to show the live batch.
    #[test]
    fn fit_tab_shows_a_banner_and_a_different_grid_when_scrubbing() {
        let mut app = fixture_app_with_ridge();
        app.set_frame_summary(5, 1, 0);
        let scrubbed = fixture_recording(8); // deliberately a different grid
        app.set_scrub_recording(2, Some(17), scrubbed);

        let out = render(&mut app, 100, 30);
        assert!(
            out.contains("SCRUBBED"),
            "must announce this is a replayed batch:\n{out}"
        );
        assert!(
            out.contains("3/5"),
            "must show a 1-based frame position out of the retained total:\n{out}"
        );
        assert!(
            out.contains("batch 17"),
            "must show the scrubbed frame's original batch number:\n{out}"
        );
        insta::assert_snapshot!(out);
    }

    /// A Fit-tab body exactly one row tall must still show the banner rather
    /// than spend that one row on a heatmap too short to read anyway — a
    /// review of the original fix found the old `area.height > 1` guard
    /// silently dropped the banner here, leaving the replayed grid
    /// indistinguishable from live in that one edge case.
    #[test]
    fn the_scrub_banner_still_shows_when_the_fit_tab_body_is_one_row_tall() {
        let mut app = fixture_app_with_ridge();
        app.set_frame_summary(5, 1, 0);
        app.set_scrub_recording(2, Some(17), fixture_recording(8));
        // 1 tab-bar row + 1 body row + 1 status row = 3.
        let out = render(&mut app, 100, 3);
        assert!(
            out.contains("SCRUBBED"),
            "the banner must win the body's only row rather than leave it \
             looking like an unlabeled (and possibly mistaken-for-live) heatmap:\n{out}"
        );
    }

    /// Clearing the scrub (as `>` past the last retained frame does) must
    /// drop the banner and go back to rendering the live recording.
    #[test]
    fn clearing_the_scrub_returns_to_the_live_view() {
        let mut app = fixture_app_with_ridge();
        app.set_frame_summary(5, 1, 0);
        app.set_scrub_recording(2, Some(17), fixture_recording(8));
        app.clear_scrub();

        let out = render(&mut app, 100, 30);
        assert!(
            !out.contains("SCRUBBED"),
            "banner must be gone once scrub is cleared:\n{out}"
        );
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

    /// Once `App::set_final` (the same setter `CalibDash::show_final` uses)
    /// has run, the Tolerances tab must actually render the m/z, mobility
    /// and RT numbers rather than the "not wired through yet" placeholder —
    /// that placeholder text is exactly what regressed before this fix (the
    /// summary was computed and thrown away).
    #[test]
    fn tolerances_tab_renders_the_summary_once_set_final_has_run() {
        let mut app = fixture_app_with_ridge();
        let rec = {
            let mut state = CalibrationState::new(16, (0.0, 16.0), (0.0, 48.0), 1).unwrap();
            state.update(ridge_points().into_iter()).unwrap();
            let mut rec = FitRecording::new(16);
            state.fit_with(&mut rec, ObserveOpts::NONE);
            state.measure_ridge_width_with(0.3, &mut rec);
            rec
        };
        app.set_final(
            rec,
            crate::ToleranceSummary {
                mz_ppm: (-8.5, 9.5),
                mobility_pct: (-3.0, 3.0),
                rt_seconds: (12.5, 12.5),
                n_calibrants: 42,
            },
        );
        app.set_tab(Tab::Tolerances);
        let out = render(&mut app, 100, 30);
        assert!(
            out.contains("-8.5") && out.contains("9.5"),
            "m/z tolerance missing:\n{out}"
        );
        assert!(
            out.contains("-3.0") && out.contains("3.0"),
            "mobility tolerance missing:\n{out}"
        );
        assert!(out.contains("12.5"), "RT tolerance missing:\n{out}");
        assert!(
            out.contains('±'),
            "RT tolerance must read as symmetric:\n{out}"
        );
        assert!(out.contains("42"), "n_calibrants missing:\n{out}");
        insta::assert_snapshot!(out);
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

    /// Item 1 adds several new size guards (the framed heatmap's gutter, its
    /// x-tick-label row, the `?` overlay's `Clear`) around specific
    /// thresholds (`area.width < 3`, `block_h < 3`, `gutter_w + 3`,
    /// `area.height >= 4`). An exhaustive sweep of every width/height up to
    /// 12 — comfortably past all of those — is cheap and catches an
    /// off-by-one at any of them that a handful of hand-picked sizes could
    /// miss. Every tab and both DP-pane states are covered too, since the
    /// same small area feeds `draw_dp_pane` when it's on.
    #[test]
    fn every_small_terminal_size_survives_every_tab() {
        use ratatui::crossterm::event::{
            KeyCode,
            KeyEvent,
            KeyModifiers,
        };
        let mut app = fixture_app_with_ridge();
        app.handle_key(KeyEvent::new(KeyCode::Char('d'), KeyModifiers::NONE)); // dp_pane on
        for tab in Tab::ALL {
            app.set_tab(tab);
            for w in 0..=12u16 {
                for h in 0..=12u16 {
                    render(&mut app, w, h);
                }
            }
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

    /// Item 0 regression: `grid_to_screen` used to map grid row 0 (the
    /// *lowest* observed RT) to screen row 0, but terminal rows grow
    /// downward — so a monotonically increasing calibration
    /// (`fixture_recording`'s ridge is `observed = 2 * library`, unambiguously
    /// increasing and with lookback 5, so the DP chains every point into one
    /// `MARK_DP_CHAIN` ('O') path) rendered as a *descending* line instead of
    /// an ascending one. The fix must make the marked cells' screen-row index
    /// descend (move toward the top of the canvas, i.e. get smaller) as the
    /// column index grows, not the reverse.
    #[test]
    fn increasing_ridge_renders_ascending_not_descending() {
        let bins = 20;
        let mut app = App::new(bins);
        *app.recording_mut() = fixture_recording(bins);
        app.set_layer(Layer::Path);

        let out = render(&mut app, 100, 30);
        // `O` is the DP-chain glyph (Item 3); collect the screen (row, col)
        // of every one, keyed by column, then walk columns left to right.
        // Skips the Fit subtitle row: the Path layer's gloss (`fit_subtitle`/
        // `layer_gloss`) spells out the glyphs it explains ("O chosen, X
        // greedy tail"), so it legitimately contains a literal `O` that a
        // whole-output scan can't tell apart from an actual marked cell.
        let mut by_col: std::collections::BTreeMap<usize, usize> =
            std::collections::BTreeMap::new();
        for (row, line) in out.lines().enumerate() {
            if line.contains("Showing:") {
                continue;
            }
            for (col, ch) in line.chars().enumerate() {
                if ch == 'O' {
                    by_col.entry(col).or_insert(row);
                }
            }
        }
        assert!(
            by_col.len() >= 2,
            "expected at least two distinct columns of DP-chain marks:\n{out}"
        );
        let mut rows: Vec<usize> = by_col.into_values().collect();
        let first = *rows.first().unwrap();
        let last = *rows.last().unwrap();
        assert!(
            last < first,
            "an increasing ridge must ascend on screen (later columns at a \
             smaller screen-row index) — first column's row {first}, last \
             column's row {last}:\n{out}"
        );
        // Non-increasing throughout, not just at the two ends — a real
        // transpose/flip bug would show up as a reversal somewhere in the
        // middle even if the endpoints happened to differ correctly.
        rows.dedup();
        assert!(
            rows.windows(2).all(|w| w[0] >= w[1]),
            "screen-row index must never increase as the column index grows:\n{out}"
        );
    }

    /// `compose_cell`'s per-cell mark-priority rule, exercised directly
    /// rather than through a fixture engineered to make two grid rows
    /// collapse into one terminal cell. Before this pass, that fixture
    /// forced a *suppressed* mark and a *DP-chain* mark into the same cell —
    /// the one collision reachable at all, since `Stage::Path` used to show
    /// both overlays at once. Suppressed and DP-chain can no longer collide:
    /// they are two different, mutually exclusive `Layer`s now, so the only
    /// same-layer collision left is a DP node against a tail node within the
    /// `Path` layer (`mark_path` marks both), and this pins that the DP
    /// chain wins regardless of which half carried which mark.
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

        // Order-independence: the winner must not depend on which half
        // happened to carry which mark.
        let (glyph2, ..) = compose_cell(Mark::DpNode, Mark::Tail, 0.0, 0.0, 1.0, Layer::Path);
        assert_eq!(glyph2, "O");
    }

    /// The empty-cell trap the brief calls out by name: `heat_glyph` returns
    /// a plain space for zero weight, and a *reversed* space renders as a
    /// solid block — the single highest-density glyph this module draws —
    /// which would lie about a zero-weight cell a region mark (curve/ridge/
    /// suppressed) lands on. `compose_cell` must substitute the minimum-ink
    /// `\u{b7}` glyph instead, still reversed.
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

        // Both halves marked and both zero weight must not let the combined
        // glyph regress to a bare space either.
        let (glyph_both, _color2, modifier_both) =
            compose_cell(Mark::Region, Mark::Region, 0.0, 0.0, 1.0, Layer::Ridge);
        assert_ne!(glyph_both, " ");
        assert_eq!(glyph_both, "\u{b7}");
        assert!(modifier_both.contains(Modifier::REVERSED));
    }

    /// End-to-end counterpart of the two `compose_cell` unit tests above:
    /// the real Fit tab, with a real fixture and the Curve layer active,
    /// must not paint any reversed space anywhere on screen, and must
    /// exercise the `\u{b7}` fallback at least once — the fitted curve is
    /// evaluated at every display column regardless of whether any
    /// calibrant actually falls there, so with only 16 bins spread across a
    /// ~90-column canvas, several columns land on zero-weight grid cells.
    #[test]
    fn curve_layer_never_renders_a_reversed_space_on_a_zero_weight_cell() {
        let mut app = fixture_app_with_ridge();
        app.set_layer(Layer::Curve);
        let mut t = Terminal::new(TestBackend::new(100, 30)).expect("test terminal");
        t.draw(|f| draw(f, &mut app)).expect("draw");
        let buf = t.backend().buffer().clone();

        let mut saw_reversed_dot = false;
        for y in 0..buf.area.height {
            for x in 0..buf.area.width {
                let cell = &buf[(x, y)];
                if cell.modifier.contains(Modifier::REVERSED) {
                    assert_ne!(
                        cell.symbol(),
                        " ",
                        "a reversed space at ({x}, {y}) reads as a solid block \
                         (maximum density), lying about a zero-weight cell the \
                         curve mark landed on"
                    );
                    saw_reversed_dot |= cell.symbol() == "\u{b7}";
                }
            }
        }
        assert!(
            saw_reversed_dot,
            "expected the curve to cross at least one zero-weight column and \
             render the reversed minimum-ink glyph there"
        );
    }

    /// The path layer's two node kinds (Item 3's `O`/`X` marker set) must
    /// both still be visible at once — that is the entire reason they are
    /// glyph-distinguished within one layer rather than folded into a single
    /// "path" mark. `fixture_app_with_ridge`'s dip/stray points (see that
    /// fixture's own doc comment) are built specifically so the DP chain
    /// stops short and a tail gets greedily attached, so both glyphs are
    /// guaranteed to appear.
    #[test]
    fn path_layer_shows_the_dp_chain_and_tail_glyphs_distinctly() {
        let mut app = fixture_app_with_ridge();
        app.set_layer(Layer::Path);
        let out = render(&mut app, 100, 30);
        assert!(out.contains('O'), "missing the DP-chain glyph:\n{out}");
        assert!(out.contains('X'), "missing the tail glyph:\n{out}");
    }

    /// `mark_suppressed`/`mark_ridge` write `Mark::Region`, which — unlike
    /// the old glyph-replacing overlay — is invisible in a plain-text
    /// (`.symbol()`-only) snapshot whenever it lands on a nonzero-weight
    /// cell. That is actually the *common* case for `mark_suppressed`, since
    /// it only ever marks cells with `weight > 0.0` to begin with (see its
    /// own guard), so `fit_layer_Suppressed`'s snapshot can look
    /// byte-identical to `fit_layer_None`'s in plain text and still be
    /// correct — the density glyph is preserved on purpose, and only the
    /// (untested-by-`.symbol()`) `REVERSED` style differs. This pins that
    /// the two functions actually mark something in the standard fixture,
    /// directly, rather than relying on a snapshot that cannot see it.
    #[test]
    fn mark_suppressed_and_mark_ridge_mark_at_least_one_cell() {
        let app = fixture_app_with_ridge();
        let rec = app.recording();
        let dims = Dims {
            w: 100,
            h: 15,
            bins: rec.geom().bins,
            disp_rows: 30,
            disp_cols: 100,
        };

        let mut suppressed = vec![Mark::None; dims.w * dims.h * 2];
        mark_suppressed(&mut suppressed, dims, rec);
        assert!(
            suppressed.contains(&Mark::Region),
            "expected at least one suppressed cell in the fixture"
        );

        let mut ridge = vec![Mark::None; dims.w * dims.h * 2];
        mark_ridge(&mut ridge, dims, rec);
        assert!(
            ridge.contains(&Mark::Region),
            "expected at least one ridge-band cell in the fixture"
        );
    }

    // ---- scaled_u64: NaN must not read as convergence ----

    #[test]
    fn scaled_u64_carries_the_last_finite_value_across_a_nan() {
        // Batch 0's own delta and any failed batch's metrics are NaN — this
        // must not draw as the scale's `0`, which would look identical to
        // "the curve stopped moving."
        let values = [10.0, f64::NAN, 10.0];
        let scaled = scaled_u64(&values);
        assert_eq!(
            scaled[1], scaled[0],
            "a NaN sample must repeat the prior finite value, not read as 0: {scaled:?}"
        );
        assert_eq!(scaled[2], scaled[0]);
    }

    #[test]
    fn scaled_u64_reports_zero_only_before_any_finite_value_is_seen() {
        // With no prior finite value to hold, a leading run of NaN (the
        // series' first samples) has nothing to carry forward — this is the
        // one case where `0` is the honest answer ("nothing to show yet"),
        // not a misread of "converged".
        let values = [f64::NAN, f64::NAN, 10.0];
        let scaled = scaled_u64(&values);
        assert_eq!(scaled[0], 0);
        assert_eq!(scaled[1], 0);
        assert!(
            scaled[2] > 0,
            "the first real sample must still scale normally"
        );
    }

    #[test]
    fn scaled_u64_of_an_all_nan_series_is_flat_zero_not_a_panic() {
        let scaled = scaled_u64(&[f64::NAN, f64::NAN, f64::NAN]);
        assert_eq!(scaled, vec![0, 0, 0]);
    }

    // ---- status line (Item 2): weight not punctuation, per-tab bindings,
    // `cnt` only when pending, `? keys` pinned right, and a degrade sequence
    // that drops action words, then the tab-local group, before anything
    // ever overflows `area.width`. -----------------------------------------

    fn long_prefix_app() -> App {
        use ratatui::crossterm::event::{
            KeyCode,
            KeyEvent,
            KeyModifiers,
        };
        let mut app = App::new(10);
        app.set_stage(Stage::Suppressed); // the longest stage label
        for c in "123456789".chars() {
            app.handle_key(KeyEvent::new(KeyCode::Char(c), KeyModifiers::NONE));
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

    #[test]
    fn status_line_shows_batch_stage_and_hides_the_permanent_cnt_column() {
        let mut app = App::new(10);
        app.set_stage(Stage::Ridge);
        let status = status_line(&mut app, 100);
        assert!(status.contains("b0"), "{status:?}");
        assert!(status.contains("Ridge"), "{status:?}");
        // The line this replaced spent six columns permanently reporting
        // `cnt:-`; there is no pending count here, so nothing should.
        assert!(!status.contains("cnt"), "{status:?}");
    }

    /// `REVERSED` — a mode signal, like the selected tab and the scrub
    /// banner (Item 4) — not a color, is what makes a pending count
    /// unmistakable. The plain-text `render` harness can't see style, so
    /// this reads the `TestBackend` buffer directly.
    #[test]
    fn status_line_pending_count_is_reversed() {
        use ratatui::crossterm::event::{
            KeyCode,
            KeyEvent,
            KeyModifiers,
        };
        let mut app = App::new(10);
        app.handle_key(KeyEvent::new(KeyCode::Char('4'), KeyModifiers::NONE));
        assert_eq!(app.pending_count(), Some(4));

        let mut t = Terminal::new(TestBackend::new(60, 3)).expect("test terminal");
        t.draw(|f| draw(f, &mut app)).expect("draw");
        let buf = t.backend().buffer().clone();
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
        app.set_tab(Tab::Convergence);
        let status = status_line(&mut app, 200);
        for word in ["stage", "frame", "layer", "dp"] {
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

    #[test]
    fn status_line_never_panics_across_a_range_of_widths() {
        for width in [10, 20, 30, 44, 60, 78, 90, 100, 120, 150, 200] {
            let _ = status_line(&mut long_prefix_app(), width);
        }
    }

    // `fit_status_hints` itself, exercised directly rather than through the
    // full `render` pipeline: the exact widths at which each stage of the
    // degrade kicks in are arithmetic on the two fixtures below, easiest to
    // pin without also depending on how wide `draw_status_line`'s left/right
    // columns happen to be.

    #[test]
    fn fit_status_hints_shows_full_text_with_room() {
        let tab_local: &[(&str, &str)] = &[("n", "next"), ("r", "run")];
        let global: &[(&str, &str)] = &[("h l", "tab"), ("q", "detach"), ("^C", "abort")];
        let text = line_text(&fit_status_hints(tab_local, global, 200));
        assert!(text.contains("n next"), "{text:?}");
        assert!(text.contains("h l tab"), "{text:?}");
        assert!(text.contains("^C abort"), "{text:?}");
    }

    #[test]
    fn fit_status_hints_drops_action_words_before_the_tab_local_group() {
        let tab_local: &[(&str, &str)] = &[("n", "next"), ("r", "run")];
        let global: &[(&str, &str)] = &[("h l", "tab"), ("q", "detach"), ("^C", "abort")];
        // Full text ("n next · r run · h l tab · q detach · ^C abort") is 46
        // chars; the keys-only line ("n · r · h l · q · ^C") is 20. A width
        // of 30 fits the second but not the first.
        let text = line_text(&fit_status_hints(tab_local, global, 30));
        assert!(!text.contains("next"), "{text:?}");
        assert!(!text.contains("detach"), "{text:?}");
        assert!(text.contains('n'), "{text:?}");
        assert!(text.contains("h l"), "{text:?}");
    }

    #[test]
    fn fit_status_hints_drops_the_tab_local_group_before_the_global_one() {
        let tab_local: &[(&str, &str)] = &[("n", "next"), ("r", "run")];
        let global: &[(&str, &str)] = &[("h l", "tab"), ("q", "detach"), ("^C", "abort")];
        // The keys-only line (20 chars) does not fit 15; the global-only
        // keys line ("h l · q · ^C", 12 chars) does.
        let text = line_text(&fit_status_hints(tab_local, global, 15));
        assert!(!text.contains('n'), "tab-local `n` must be gone: {text:?}");
        assert!(text.contains("h l"), "global keys must survive: {text:?}");
    }

    #[test]
    fn fit_status_hints_is_empty_when_nothing_fits() {
        let tab_local: &[(&str, &str)] = &[("n", "next")];
        let global: &[(&str, &str)] = &[("h l", "tab")];
        let text = line_text(&fit_status_hints(tab_local, global, 0));
        assert!(text.is_empty(), "{text:?}");
    }
}
