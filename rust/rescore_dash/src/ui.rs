//! Rendering.
//!
//! Every panel is a lookup into the precomputed [`crate::precompute::Dashboard`]:
//! histogram bins, axis ranges, curves, threshold tables, panel titles and
//! table cells were all computed before the TUI opened. A frame indexes arrays,
//! fills two reused point buffers, and draws — a keystroke costs a redraw and
//! not a re-scan.

use crate::app::{
    App,
    Scratch,
    Tab,
};
use crate::precompute::{
    Dashboard,
    FeatureColumn,
};
use crate::transform::{
    Axis,
    YTransform,
};
use ratatui::layout::{
    Constraint,
    Direction,
    Layout,
    Rect,
};
use ratatui::style::{
    Color,
    Modifier,
    Style,
};
use ratatui::text::Line;
// Aliased: `Axis` in this crate is a histogram's x transform, which every panel
// below also passes around.
use ratatui::widgets::{
    Axis as PlotAxis,
    Block,
    Borders,
    Chart,
    Dataset,
    GraphType,
    Paragraph,
    Row,
    Table,
    Tabs,
};
use ratatui::{
    Frame,
    symbols,
};

const TARGET_COLOR: Color = Color::Cyan;
const DECOY_COLOR: Color = Color::Magenta;

/// One line series. Every plot on every tab is drawn with these settings.
fn series<'a>(name: &'a str, color: Color, data: &'a [(f64, f64)]) -> Dataset<'a> {
    Dataset::default()
        .name(name)
        .marker(symbols::Marker::Braille)
        .graph_type(GraphType::Line)
        .style(Style::default().fg(color))
        .data(data)
}

/// The axis settings a histogram panel reads off `App`.
///
/// Gathered before the `App` is destructured for disjoint field borrows, which
/// is why they travel as one value rather than as `&App`.
#[derive(Clone, Copy)]
struct HistAxes {
    x: Axis,
    y: YTransform,
    clip: bool,
}

impl HistAxes {
    fn of(app: &App) -> Self {
        Self {
            x: app.x(),
            y: app.y(),
            clip: app.clip(),
        }
    }
}

/// A bordered panel carrying a message instead of a plot.
fn empty_panel(frame: &mut Frame, area: Rect, title: &str, message: &str) {
    frame.render_widget(
        Paragraph::new(message).block(Block::default().borders(Borders::ALL).title(title)),
        area,
    );
}

/// Top-level: tab bar, body, help line.
pub fn draw(frame: &mut Frame, app: &mut App) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Min(0),
            Constraint::Length(1),
        ])
        .split(frame.area());

    let titles: Vec<Line> = Tab::ALL.iter().map(|t| Line::from(t.title())).collect();
    let selected = Tab::ALL.iter().position(|t| *t == app.tab()).unwrap_or(0);
    // Full width, which is why the run-level basis goes here and not into a
    // histogram subtitle that the panel border would truncate.
    let banner = format!("rescore — {}", app.dash.basis());
    frame.render_widget(
        Tabs::new(titles)
            .select(selected)
            .block(Block::default().borders(Borders::ALL).title(banner))
            .highlight_style(Style::default().add_modifier(Modifier::BOLD)),
        chunks[0],
    );

    match app.tab() {
        Tab::Overview => draw_overview(frame, app, chunks[1]),
        Tab::Fdr => draw_fdr(frame, app, chunks[1]),
        Tab::Calibration => draw_calibration(frame, app, chunks[1]),
        Tab::Features => draw_features(frame, app, chunks[1]),
    }

    frame.render_widget(Paragraph::new(help_line(app)), chunks[2]);
}

/// The bottom key-hint line; a different message while the filter box is open.
///
/// Per-tab rather than one exhaustive list. Listing every binding everywhere
/// both overflowed a narrow terminal and advertised keys that do nothing on the
/// current tab — `s sort` on the FDR curve, `z zoom` on the feature table.
fn help_line(app: &App) -> String {
    if app.filter_editing() {
        return format!("/{}  Enter apply | Esc clear", app.filter());
    }
    let hist = || {
        format!(
            "x/X {}  y/Y {}  c clip:{}",
            app.x().label(),
            app.y().label(),
            if app.clip() { "on" } else { "off" },
        )
    };
    // `None` rather than an empty string: an empty tail leaves the separators on
    // either side of it doubled up.
    let tail = match app.tab() {
        Tab::Overview => Some(hist()),
        Tab::Fdr => Some(format!("z/Z q<={}", app.dash.q_curve(app.q_zoom()).zoom)),
        Tab::Calibration => None,
        Tab::Features => Some(format!(
            "j/k rows  {}  s/S sort:{}  / filter",
            hist(),
            app.sort_summary()
        )),
    };
    match tail {
        Some(tail) => format!("h/l tabs  {tail}  q quit"),
        None => "h/l tabs  q quit".to_string(),
    }
}

/// Render column `column`'s stored histogram.
///
/// Takes the column rather than the looked-up pieces: the title, the subtitle
/// and the bins are three lookups keyed the same way, and having each caller do
/// them was three chances to draw one panel's title over another's counts.
///
/// Tolerates an empty or all-zero histogram: the y bound falls back to `1.0`
/// rather than degenerating to `[0.0, 0.0]`.
///
/// Title on the top border, subtitle on the bottom, because the counts they
/// carry are disjoint: `dropped` is what the transform excluded, the
/// out-of-range count is what survived it and then fell outside the axis.
fn draw_hist(
    frame: &mut Frame,
    area: Rect,
    dash: &Dashboard,
    column: usize,
    axes: HistAxes,
    scratch: &mut Scratch,
) {
    let HistAxes { x, y, clip } = axes;
    let hist = dash.hist(column, x, clip);
    scratch.target.clear();
    scratch.decoy.clear();
    let mut ymax = 0.0f64;
    for (i, (&t, &d)) in hist.target.iter().zip(hist.decoy).enumerate() {
        let px = hist.bin_center(i);
        let (ty, dy) = (y.apply(t, hist.n_target), y.apply(d, hist.n_decoy));
        ymax = ymax.max(ty).max(dy);
        scratch.target.push((px, ty));
        scratch.decoy.push((px, dy));
    }
    let ymax = if ymax > 0.0 { ymax } else { 1.0 };

    let datasets = vec![
        series("target", TARGET_COLOR, &scratch.target),
        series("decoy", DECOY_COLOR, &scratch.decoy),
    ];

    // Axis meanings go on the block's borders, not into `PlotAxis::title`, which
    // ratatui draws inline over the plot's last row and column.
    let mid = (hist.lo + hist.hi) / 2.0;
    let chart = Chart::new(datasets)
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(dash.title(column, x))
                .title_bottom(dash.subtitle(column, x, clip)),
        )
        .x_axis(PlotAxis::default().bounds([hist.lo, hist.hi]).labels([
            format!("{:.3}", hist.lo),
            format!("{mid:.3}"),
            format!("{:.3}", hist.hi),
        ]))
        .y_axis(PlotAxis::default().bounds([0.0, ymax]));
    frame.render_widget(chart, area);
}

/// Score summary line plus the discriminant-score histogram.
fn draw_overview(frame: &mut Frame, app: &mut App, area: Rect) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Length(3), Constraint::Min(0)])
        .split(area);

    let axes = HistAxes::of(app);
    // Disjoint field borrows: the dashboard is read while the scratch buffers
    // are written. Going through `&mut self` accessors would borrow the whole
    // `App` twice over.
    let App { dash, scratch, .. } = app;

    frame.render_widget(Paragraph::new(dash.overview_header()), chunks[0]);
    let score = dash.score_column();
    draw_hist(frame, chunks[1], dash, score, axes, scratch);
}

/// Targets-passing-vs-q-value curve, plus a threshold table at fixed FDR
/// cutoffs.
fn draw_fdr(frame: &mut Frame, app: &mut App, area: Rect) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Percentage(60), Constraint::Percentage(40)])
        .split(area);

    let dash = &app.dash;
    let curve = dash.q_curve(app.q_zoom());

    // Both axes are labelled here where the histogram panel labels only x: the
    // zoom changes both bounds at once, and an unlabelled curve gives no way to
    // tell which zoom is on screen.
    let chart = Chart::new(vec![series("targets", TARGET_COLOR, &curve.points)])
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(curve.title.as_str())
                .title_bottom("x: q-value   y: targets passing   z/Z zoom"),
        )
        .x_axis(
            PlotAxis::default()
                .bounds([0.0, curve.zoom])
                .labels(curve.x_labels.iter().map(String::as_str)),
        )
        .y_axis(
            PlotAxis::default()
                .bounds([0.0, curve.ymax])
                .labels(curve.y_labels.iter().map(String::as_str)),
        );
    frame.render_widget(chart, chunks[0]);

    let rows: Vec<Row> = dash
        .fdr_rows()
        .iter()
        .map(|r| Row::new(r.iter().map(String::as_str)))
        .collect();
    let table = Table::new(
        rows,
        [
            Constraint::Length(8),
            Constraint::Length(10),
            Constraint::Length(10),
            Constraint::Length(10),
        ],
    )
    .header(Row::new(vec!["q", "total", "targets", "decoys"]))
    .block(Block::default().borders(Borders::ALL).title("thresholds"));
    frame.render_widget(table, chunks[1]);
}

/// Target-vs-decoy PP plot (empirical CDF of one against the other) with the
/// y = x reference line.
fn draw_calibration(frame: &mut Frame, app: &mut App, area: Rect) {
    // The plot's axes are (decoy CDF, target CDF): good target/decoy separation
    // pulls the curve BELOW the y = x diagonal (targets need a higher score
    // threshold before as many of them have "arrived" as decoys have), so the
    // title carries that reading hint.
    let title = "Calibration (below y=x = separation)";

    let curve = &app.dash.pp_curve;
    if curve.is_empty() {
        empty_panel(frame, area, title, "PP plot needs both targets and decoys");
        return;
    }

    let reference = [(0.0, 0.0), (1.0, 1.0)];
    let datasets = vec![
        series("PP", TARGET_COLOR, curve),
        series("y = x", Color::DarkGray, &reference),
    ];
    let chart = Chart::new(datasets)
        .block(Block::default().borders(Borders::ALL).title(title))
        .x_axis(PlotAxis::default().title("decoy CDF").bounds([0.0, 1.0]))
        .y_axis(PlotAxis::default().title("target CDF").bounds([0.0, 1.0]));
    frame.render_widget(chart, area);
}

/// Sortable feature table plus the selected feature's stored histogram.
fn draw_features(frame: &mut Frame, app: &mut App, area: Rect) {
    let chunks = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(55), Constraint::Percentage(45)])
        .split(area);

    let axes = HistAxes::of(app);
    let sort = app.sort_summary();
    // Destructured for the same reason as `draw_overview`, but three ways: the
    // table widget holds `dash`'s cells, the table state needs `&mut`, and the
    // histogram panel writes `scratch`.
    let App {
        dash,
        scratch,
        table_state,
        visible,
        cursor,
        ..
    } = app;
    let selected = visible.get(*cursor).copied();

    // Feature names share long prefixes (`ms1_*`, `ms2_*`, `lazyscore_*`), so
    // a fixed-width name column shows an identical prefix on every row. `Min`
    // grows into whatever the numeric columns leave. Not `Fill(1)` or a larger
    // `Min`: once the name width plus six `Length(9)` exceeds the panel, the
    // solver shrinks the `Length` columns instead, truncating the numbers.
    // `9` fits the widest realistic cell (`-307.0000`).
    let widths = [
        Constraint::Min(16),
        Constraint::Length(9),
        Constraint::Length(9),
        Constraint::Length(9),
        Constraint::Length(9),
        Constraint::Length(9),
        Constraint::Length(9),
    ];
    let title = format!(
        "Features ({}/{}) sort:{sort}",
        visible.len(),
        dash.n_features()
    );

    // Rows borrow the dashboard's pre-formatted cells, so building the table
    // copies no strings.
    let rows: Vec<Row> = visible
        .iter()
        .map(|&j| Row::new(dash.cells()[j].iter().map(String::as_str)))
        .collect();
    let table = Table::new(rows, widths)
        .header(Row::new(FeatureColumn::ALL.map(FeatureColumn::label)))
        .row_highlight_style(Style::default().add_modifier(Modifier::REVERSED))
        .block(Block::default().borders(Borders::ALL).title(title));

    table_state.select(selected.is_some().then_some(*cursor));
    frame.render_stateful_widget(table, chunks[0], table_state);

    match selected {
        Some(j) => draw_hist(frame, chunks[1], dash, j, axes, scratch),
        None => empty_panel(frame, chunks[1], "Features", "no feature selected"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::app::tests::{
        code,
        dashboard,
        key,
    };
    use crate::app::{
        App,
        Tab,
    };
    use crate::precompute::tests::fixture_from_columns;
    use ratatui::Terminal;
    use ratatui::backend::TestBackend;
    use ratatui::crossterm::event::KeyCode;

    fn render(app: &mut App, width: u16, height: u16) -> String {
        let backend = TestBackend::new(width, height);
        let mut terminal = Terminal::new(backend).expect("test terminal");
        terminal.draw(|f| draw(f, app)).expect("draw");
        terminal
            .backend()
            .buffer()
            .content()
            .iter()
            .map(|cell| cell.symbol())
            .collect()
    }

    fn go_to(app: &mut App, tab: Tab) {
        while app.tab() != tab {
            app.handle_key(key('l'));
        }
    }

    /// 40 rows x 2 features: targets score high, decoys low, and one feature is
    /// wholly missing.
    fn fixture_app() -> App {
        let alpha: Vec<f64> = (0..40)
            .map(|i| if i % 2 == 0 { i as f64 } else { -(i as f64) })
            .collect();
        App::new(dashboard(
            &["alpha_score", "all_nan_feature"],
            &[&alpha, &[f64::NAN; 40]],
        ))
    }

    #[test]
    fn features_tab_shows_feature_names() {
        let mut app = fixture_app();
        go_to(&mut app, Tab::Features);
        assert!(render(&mut app, 120, 40).contains("alpha_score"));
    }

    /// A smoke test over the whole lookup surface — every axis and both clip
    /// settings on every tab. It asserts only that a frame draws, and that the
    /// tab that drew it says which one it is: a slot the precompute left
    /// unplottable or an axis range it could not derive panics on index here
    /// instead of in front of a user. What is *on* the frame is each panel's own
    /// test.
    #[test]
    fn no_axis_or_clip_setting_panics_on_any_tab() {
        for tab in Tab::ALL {
            let mut app = fixture_app();
            go_to(&mut app, tab);
            for _ in 0..crate::transform::Axis::ALL.len() {
                for _ in 0..2 {
                    let text = render(&mut app, 120, 40);
                    assert!(text.contains(tab.title()), "{tab:?} did not name itself");
                    app.handle_key(key('c'));
                }
                app.handle_key(key('x'));
            }
        }
    }

    /// Each FDR zoom must draw, and must say on screen which range it is: the
    /// curve shape alone does not distinguish `q <= 0.05` from `q <= 1`, so an
    /// unlabelled zoom is worse than none.
    #[test]
    fn every_fdr_zoom_renders_and_names_its_range() {
        let mut app = fixture_app();
        go_to(&mut app, Tab::Fdr);
        for _ in 0..app.dash.n_q_zooms() {
            let want = format!("q <= {}", app.dash.q_curve(app.q_zoom()).zoom);
            let text = render(&mut app, 120, 40);
            assert!(text.contains(&want), "expected {want:?} in: {text}");
            app.handle_key(key('z'));
        }
    }

    /// The features tab is the only one with a sort, and `S` only flips the
    /// direction — which has to appear on screen, or the key silently changes
    /// what the table means.
    #[test]
    fn the_features_tab_names_its_sort_key_and_direction() {
        let mut app = fixture_app();
        go_to(&mut app, Tab::Features);
        let text = render(&mut app, 160, 30);
        assert!(text.contains("name a-z"), "sort direction missing: {text}");

        app.handle_key(key('S'));
        let flipped = render(&mut app, 160, 30);
        assert!(
            flipped.contains("name z-a"),
            "S must change what is shown: {flipped}"
        );
    }

    /// The calibration tab has no per-tab keys, so its help line has no tail.
    /// It must not leave a hole where one would go.
    #[test]
    fn a_tab_without_its_own_keys_has_no_doubled_separator() {
        let mut app = fixture_app();
        go_to(&mut app, Tab::Calibration);
        assert_eq!(help_line(&app), "h/l tabs  q quit");
    }

    /// Selecting the wholly-NaN feature must draw the panel and say why it is
    /// empty, rather than plotting a fabricated range.
    #[test]
    fn selecting_an_all_nan_feature_explains_the_empty_panel() {
        let mut app = fixture_app();
        go_to(&mut app, Tab::Features);
        // The table is name-sorted, so walk to the feature rather than
        // assuming where the default sort put it.
        while app.selected_feature() != Some(1) {
            app.handle_key(key('j'));
        }
        let text = render(&mut app, 160, 40);
        assert!(
            text.contains("nothing this transform"),
            "expected the unplottable diagnostic, got: {text}"
        );
    }

    /// A filter matching nothing leaves no selection; the histogram panel has
    /// to say so rather than indexing past the end of `visible`.
    #[test]
    fn rendering_survives_a_filter_that_matches_nothing() {
        let mut app = fixture_app();
        go_to(&mut app, Tab::Features);
        app.handle_key(key('/'));
        for c in "zzz".chars() {
            app.handle_key(key(c));
        }
        app.handle_key(code(KeyCode::Enter));
        assert!(render(&mut app, 120, 40).contains("no feature selected"));
    }

    /// 100x20 backend, 40 features, ~13-row viewport: reproduces the
    /// scenario where `TableState::default()` every frame pinned the
    /// selected row to the bottom edge instead of letting the viewport
    /// scroll. The offset must persist (and advance) across renders as the
    /// cursor moves.
    #[test]
    fn features_table_offset_persists_and_advances_across_renders() {
        let n_features = 40;
        let n_rows = 4;
        let names: Vec<String> = (0..n_features).map(|i| format!("feature_{i:03}")).collect();
        let columns: Vec<Vec<f64>> = (0..n_features)
            .map(|j| (0..n_rows).map(|r| (r * n_features + j) as f64).collect())
            .collect();
        let name_refs: Vec<&str> = names.iter().map(String::as_str).collect();
        let column_refs: Vec<&[f64]> = columns.iter().map(Vec::as_slice).collect();
        let mut app = App::new(dashboard(&name_refs, &column_refs));
        go_to(&mut app, Tab::Features);

        let backend = TestBackend::new(100, 20);
        let mut terminal = Terminal::new(backend).expect("test terminal");
        terminal.draw(|f| draw(f, &mut app)).expect("first draw");
        let offset_before = app.table_state.offset();

        for _ in 0..30 {
            app.handle_key(key('j'));
        }
        terminal.draw(|f| draw(f, &mut app)).expect("second draw");
        let offset_after = app.table_state.offset();

        assert!(
            offset_after > offset_before,
            "offset must advance as the cursor scrolls past the viewport: {offset_before} -> {offset_after}"
        );
    }

    /// At 160x30 with real, ALL-lane-length feature names (21 and 19 chars —
    /// both would have collapsed to the same 16-char prefix under the old
    /// `Constraint::Length(16)` name column, `"lazyscore_vs_bas"` for the
    /// first), and a numeric value whose formatted output is wider than the old
    /// fixed 6-char column, neither may be silently truncated.
    ///
    /// `name_b` ("ms1_precursor_trace") is not the selected feature (that's
    /// `name_a`, at cursor `0`), so it appears nowhere else in the frame —
    /// the histogram panel's title only ever names the selection. Its full
    /// presence in the render is therefore evidence about the *table
    /// column's* width specifically, not a coincidental match against some
    /// other panel.
    #[test]
    fn feature_table_does_not_truncate_long_names_or_wide_numbers() {
        let name_a = "lazyscore_vs_baseline";
        let name_b = "ms1_precursor_trace";
        // Feature 0's targets are rows 0 and 2, both at -307.0.
        let mut app = App::new(dashboard(
            &[name_a, name_b],
            &[&[-307.0, 1.0, -307.0, 2.0], &[1.0, 2.0, 3.0, 4.0]],
        ));
        go_to(&mut app, Tab::Features);
        let text = render(&mut app, 160, 30);

        assert!(
            text.contains(name_b),
            "name column truncated a real feature name that appears only in the table: {text}"
        );
        assert!(
            text.contains("-307.0000"),
            "numeric column truncated the formatted cell: {text}"
        );
    }

    /// A frame has to say, somewhere the reader will see it, that the
    /// histograms came from a sample rather than every row.
    #[test]
    fn a_sampled_dashboard_says_so_on_every_tab() {
        let column: Vec<f64> = (0..2_000).map(|i| i as f64).collect();
        let dash = fixture_from_columns(&["a"], &[&column]).build_with(500);
        let mut app = App::new(dash);
        for tab in Tab::ALL {
            go_to(&mut app, tab);
            let text = render(&mut app, 160, 30);
            assert!(
                text.contains("2000 rows, histograms from a 500 sample"),
                "{tab:?} did not show the sampling basis: {text}"
            );
        }
    }
}
