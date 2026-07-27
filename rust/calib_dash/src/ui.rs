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
            ("s p c w", "overlay"),
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
        Line::raw("  s p c w    toggle suppressed/path/curve/ridge"),
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

/// Overlay priority for one terminal cell (which represents two grid rows —
/// see `heatmap_cells`). Higher wins when both grid rows the cell represents
/// carry different marks, or when a later stage's overlay lands on a cell an
/// earlier one already marked.
///
/// Ordered so evidence beats derived fit (Item 3): the curve and ridge are
/// *computed from* the suppressed grid and the DP/greedy path, so a chain or
/// tail node must still show through even where the curve or ridge happens
/// to land on the same cell — the point of showing the DP chain and greedy
/// tails at all is to let a user see where the fit disagrees with what it
/// was derived from, which a curve painted on top of them would hide.
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
/// that is a uniform block of one glyph regardless of data cannot show that,
/// or show one stage's overlay differing from another's. So the glyph
/// itself carries the signal here — heat density, suppression, path/DP-chain
/// membership, curve and ridge.
///
/// **Invariant:** every mark is identified by its glyph; color is redundant
/// reinforcement only. That is what keeps the panel readable for a
/// colorblind reader, in a monochrome screenshot, and in this module's own
/// symbol-only snapshot harness — never make two marks share a glyph and
/// rely on color alone to tell them apart.
const MARK_NONE: u8 = 0;
const MARK_SUPPRESSED: u8 = 1;
const MARK_CURVE: u8 = 2;
const MARK_RIDGE: u8 = 3;
const MARK_TAIL: u8 = 4;
const MARK_DP_CHAIN: u8 = 5;

/// The glyph half of the invariant above — matplotlib's own vocabulary for
/// these, since that is what this tool's user reads daily. All narrow-width
/// (single-cell) safe.
fn mark_glyph(m: u8) -> Option<&'static str> {
    match m {
        MARK_SUPPRESSED => Some("\u{b7}"), // ·
        MARK_CURVE => Some("\u{2500}"),    // ─ — a quarter the ink of a solid
        // block, and its staircase gaps correctly signal steepness where the
        // curve is steep, information a filled glyph would destroy.
        MARK_RIDGE => Some("~"),
        MARK_TAIL => Some("+"),
        MARK_DP_CHAIN => Some("O"),
        _ => None,
    }
}

/// The color/modifier half — reinforcement only, per the invariant above.
/// Every value is a named ANSI constant (never `Rgb`, which most terminals
/// do not remap for a light background — see `heat_color`) and none is one
/// of indices 9-15 (`Light*`/`White`), which are invisible on a light
/// background. `Modifier::REVERSED` is deliberately not used here — Item 4
/// reserves it for *mode* (selected tab, scrub banner, pending count), not
/// "what a mark is", so tail/DP-chain lean on `BOLD` instead for emphasis.
fn mark_style(m: u8) -> Option<(Color, Modifier)> {
    match m {
        MARK_SUPPRESSED => Some((Color::DarkGray, Modifier::empty())),
        MARK_CURVE => Some((Color::Cyan, Modifier::empty())),
        MARK_RIDGE => Some((Color::Magenta, Modifier::empty())),
        MARK_TAIL => Some((Color::Yellow, Modifier::BOLD)),
        MARK_DP_CHAIN => Some((Color::Green, Modifier::BOLD)),
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
    let title_right = format!(" b{} \u{b7} {} ", app.batch(), app.stage().label());
    let block = Block::bordered()
        .title_top(Line::from(title_left).left_aligned())
        .title_top(Line::from(title_right).right_aligned());
    let inner = block.inner(block_area);
    frame.render_widget(block, block_area);
    paint_heatmap(frame, inner, rec, app);

    if show_y_axis && inner.height > 0 {
        let y_step = nice_step((y_hi - y_lo).abs().max(1e-9), y_target);
        let y_decimals = axis_decimals(y_step);
        for v in axis_ticks(y_lo, y_hi, y_target) {
            let row = value_to_row(v, y_lo, y_hi, inner.height as usize);
            let abs_y = inner.y + row as u16;
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

/// Paints the half-block heatmap itself — density field plus the
/// suppressed/path/curve/ridge overlays — into `area`: the framed canvas's
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
            let (symbol, color, modifier) =
                compose_cell(upper_mark, lower_mark, upper_heat, lower_heat, max_w);
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
) -> (&'static str, Color, Modifier) {
    let upper_on = upper_mark != MARK_NONE || upper_heat > 0.0;
    let lower_on = lower_mark != MARK_NONE || lower_heat > 0.0;
    match (upper_on, lower_on) {
        (false, false) => (" ", Color::Reset, Modifier::empty()),
        (true, false) => match (mark_glyph(upper_mark), mark_style(upper_mark)) {
            (Some(sym), Some((col, m))) => (sym, col, m),
            _ => (
                "\u{2580}", // ▀
                heat_color(upper_heat, max),
                Modifier::empty(),
            ),
        },
        (false, true) => match (mark_glyph(lower_mark), mark_style(lower_mark)) {
            (Some(sym), Some((col, m))) => (sym, col, m),
            _ => (
                "\u{2584}", // ▄
                heat_color(lower_heat, max),
                Modifier::empty(),
            ),
        },
        (true, true) => {
            let winner = upper_mark.max(lower_mark);
            match (mark_glyph(winner), mark_style(winner)) {
                (Some(sym), Some((col, m))) => (sym, col, m),
                _ => {
                    let heat = upper_heat.max(lower_heat);
                    (
                        heat_glyph(heat, max),
                        heat_color(heat, max),
                        Modifier::empty(),
                    )
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
        // Item 0: re-derive the flipped row/half from `dr` via
        // `flip_display_row` rather than computing `Half` from `dr`'s own
        // parity — `grid_to_screen`'s doc comment explains why the parity
        // must come from the flipped index, not the pre-flip one.
        let (ty, half) = flip_display_row(dr, dims.disp_rows);
        raise_mark(marks, dims, tx, ty, half, MARK_CURVE);
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
    /// marks (`O`/`+` — DP chain/tail, Item 3's marker set) while
    /// `show_curve` (untouched) keeps the curve (`─`) rendering, proving the
    /// two are separate state rather than one stage-derived flag read two
    /// ways.
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
            !out.contains('O') && !out.contains('+'),
            "path overlay was toggled off, so no DP-chain or tail marks should \
             remain:\n{out}"
        );
        assert!(
            out.contains('\u{2500}'),
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
        app.set_stage(Stage::Path);

        let out = render(&mut app, 100, 30);
        // `O` is the DP-chain glyph (Item 3); collect the screen (row, col)
        // of every one, keyed by column, then walk columns left to right.
        let mut by_col: std::collections::BTreeMap<usize, usize> =
            std::collections::BTreeMap::new();
        for (row, line) in out.lines().enumerate() {
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

    /// A `bins=4` grid rendered into a body exactly 2 rows tall (`disp_rows
    /// == bins`, one-to-one) with a single display column, so both grid rows
    /// are known to land in the terminal's one heatmap character: row 0 (a
    /// point too light to survive `suppress_nonmax`'s 1.0 seed, so it is
    /// suppressed but still weight > 0 — the suppressed mark), row 1 (the
    /// sole survivor, so it is trivially the whole DP chain). Item 0's
    /// y-axis flip puts both of them in the terminal's *second* (bottom)
    /// heatmap row rather than its first — grid row 0 is the lowest observed
    /// RT, which now renders at the bottom of the canvas — so this is also a
    /// regression check for that flip, not just for `compose_cell`'s
    /// placement rule: the higher-priority mark wins the whole cell when
    /// both halves are occupied (`MARK_DP_CHAIN > MARK_SUPPRESSED`), so the
    /// character must be `O` (Item 3's DP-chain glyph), not `\u{b7}`
    /// (suppressed), and it must be on the *second* heatmap line.
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

        // 1 tab-bar row + 2 heatmap rows + 1 status row; 1 column. Body
        // width 1 is below `draw_heatmap`'s framing threshold, so this
        // renders the plain (unframed) heatmap — the same layout the old
        // assertion assumed.
        let out = render(&mut app, 1, 4);
        let heatmap_char = out.lines().nth(2).and_then(|l| l.chars().next());
        assert_eq!(
            heatmap_char,
            Some('O'),
            "the DP chain mark must win the whole cell over the suppressed mark \
             sharing its other half, got:\n{out}"
        );
        assert_eq!(
            out.lines().nth(1).and_then(|l| l.chars().next()),
            Some(' '),
            "grid row 0 is the lowest observed RT, so after the y-axis flip \
             nothing should land in the *first* heatmap row here:\n{out}"
        );
        insta::assert_snapshot!(out);
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
        for word in ["stage", "frame", "overlay", "dp"] {
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
