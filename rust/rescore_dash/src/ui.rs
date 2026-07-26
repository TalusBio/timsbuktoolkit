//! Rendering. Every panel recomputes from the `App` each frame; the expensive
//! parts (per-feature summaries) are cached on the `App` itself.

use crate::app::{
    App,
    SortKey,
    Tab,
};
use crate::stats::{
    self,
    Hist,
    N_BINS,
};
use crate::{
    curves,
    transform,
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
use ratatui::widgets::{
    Axis,
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
    frame.render_widget(
        Tabs::new(titles)
            .select(selected)
            .block(Block::default().borders(Borders::ALL).title("rescore"))
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
fn help_line(app: &App) -> String {
    if app.filter_editing() {
        format!("/{}  Enter apply | Esc clear", app.filter())
    } else {
        format!(
            "h/l tabs  j/k rows  x/X {}  y/Y {}  c clip:{}  s sort:{}  / filter  q quit",
            app.x().label(),
            app.y().label(),
            if app.clip() { "on" } else { "off" },
            app.sort_key().label(),
        )
    }
}

/// Per-bin `(x, y)` points for the target and decoy series of a histogram, `y`
/// mapped through the current [`transform::YTransform`].
#[allow(clippy::type_complexity)]
fn hist_datasets(hist: &Hist, y: transform::YTransform) -> (Vec<(f64, f64)>, Vec<(f64, f64)>) {
    let n_target = hist.n_target();
    let n_decoy = hist.n_decoy();
    let target = (0..hist.target.len())
        .map(|i| (hist.bin_center(i), y.apply(hist.target[i], n_target)))
        .collect();
    let decoy = (0..hist.decoy.len())
        .map(|i| (hist.bin_center(i), y.apply(hist.decoy[i], n_decoy)))
        .collect();
    (target, decoy)
}

/// Render a target/decoy histogram as a `Chart`. Tolerates an all-empty or
/// all-zero `hist`: the y bound falls back to `1.0` rather than degenerating
/// to `[0.0, 0.0]`.
///
/// `title` (what the panel shows: the score or feature name plus the
/// transform label) goes on the top border; `subtitle` (the drop/out-of-range
/// diagnostics) goes on the bottom border via `Block::title_bottom`. These
/// used to be crammed into one parenthetical on `title`, which read as if
/// "out of range" were a subset of "dropped" — it is not: `dropped` counts
/// values the transform excluded (including non-finite ones), and the
/// out-of-range count applies only to values that survived the transform and
/// were then excluded at binning time. Keeping them on separate lines, worded
/// as the separate things they are, is what `column_hist`'s two call sites
/// both do now.
fn draw_hist(
    frame: &mut Frame,
    area: Rect,
    title: &str,
    subtitle: &str,
    hist: &Hist,
    y: transform::YTransform,
) {
    let (target, decoy) = hist_datasets(hist, y);
    let ymax = target
        .iter()
        .chain(decoy.iter())
        .map(|&(_, v)| v)
        .fold(0.0f64, f64::max);
    let ymax = if ymax > 0.0 { ymax } else { 1.0 };

    let datasets = vec![
        Dataset::default()
            .name("target")
            .marker(symbols::Marker::Braille)
            .graph_type(GraphType::Line)
            .style(Style::default().fg(TARGET_COLOR))
            .data(&target),
        Dataset::default()
            .name("decoy")
            .marker(symbols::Marker::Braille)
            .graph_type(GraphType::Line)
            .style(Style::default().fg(DECOY_COLOR))
            .data(&decoy),
    ];

    // No axis title here: it would duplicate the block title passed in above
    // and get drawn inline over the last row of plot data.
    let mid = (hist.lo + hist.hi) / 2.0;
    let chart = Chart::new(datasets)
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(title)
                .title_bottom(subtitle),
        )
        .x_axis(Axis::default().bounds([hist.lo, hist.hi]).labels([
            format!("{:.3}", hist.lo),
            format!("{mid:.3}"),
            format!("{:.3}", hist.hi),
        ]))
        .y_axis(Axis::default().bounds([0.0, ymax]));
    frame.render_widget(chart, area);
}

/// Transform `values` with the current x-transform, bin the survivors into a
/// [`Hist`], and report how many the transform dropped.
///
/// `transform::transform_column` drops values (non-finite, or out of the
/// transform's domain), which breaks row alignment between the survivors and
/// `view.is_target`. The label vector is rebuilt in lockstep by re-applying
/// the same pointwise predicate (`XTransform::apply`) rather than reusing the
/// original `is_target`, so a dropped row's label is dropped with it instead
/// of silently sliding onto the next surviving value.
fn column_hist(app: &mut App, values: &[f64]) -> (Hist, usize) {
    let view = app.view();
    let x = app.x();
    let clip = app.clip();

    let (transformed, dropped) = transform::transform_column(x, values);
    let labels: Vec<bool> = values
        .iter()
        .zip(view.is_target.iter())
        .filter_map(|(&v, &t)| x.apply(v).map(|_| t))
        .collect();

    // Clipped: the 0.5/99.5 percentile range needs a sort, which is the whole
    // point (trim outlier tails). Unclipped: the range is just the finite
    // min/max, and `finite_range` gets there in one O(n) pass with no
    // allocation, rather than sorting the whole column just to read off its
    // endpoints (`percentile_range(.., 0.0, 100.0)` is a full sort for that).
    let range = if clip {
        stats::percentile_range(transformed.iter().copied(), 0.5, 99.5)
    } else {
        stats::finite_range(transformed.iter().copied()).map(|(lo, hi)| {
            // Mirrors `percentile_range`'s degenerate-range fallback: a
            // single distinct value (or none) would otherwise collapse the
            // histogram to `[lo, lo]`.
            if hi > lo { (lo, hi) } else { (lo, lo + 1.0) }
        })
    };

    let hist = match range {
        Some((lo, hi)) => stats::histogram(transformed.into_iter(), &labels, lo, hi, N_BINS),
        None => stats::histogram(std::iter::empty(), &[], 0.0, 1.0, N_BINS),
    };
    (hist, dropped)
}

/// The two histogram panels' shared diagnostics line: values the transform
/// excluded (`dropped`, with the non-finite subset broken out as
/// `nan_count`), and, separately, finite values that survived the transform
/// but landed outside the plotted range at binning time (`n_out`). These are
/// two different exclusions at two different pipeline steps — `n_out` is not
/// a subset of `dropped` — so they get worded as separate clauses rather than
/// nested in one parenthetical.
fn hist_subtitle(dropped: usize, nan_count: usize, n_out: usize) -> String {
    format!(
        "dropped {dropped} by transform ({nan_count} non-finite) | {n_out} finite values outside plotted range"
    )
}

/// Score summary line plus the discriminant-score histogram.
fn draw_overview(frame: &mut Frame, app: &mut App, area: Rect) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Length(3), Constraint::Min(0)])
        .split(area);

    let view = app.view();
    let n_targets = view.is_target.iter().filter(|&&t| t).count();
    let n_decoys = view.is_target.len() - n_targets;
    let thresholds = curves::threshold_table(&view.qvalue, &view.is_target, &[0.01, 0.05, 0.10]);
    let counts = thresholds
        .iter()
        .map(|&(q, total, targets, decoys)| {
            format!("q<={:.0}%: {total} (t{targets}/d{decoys})", q * 100.0)
        })
        .collect::<Vec<_>>()
        .join("   ");
    let score: Vec<f64> = view.score.iter().map(|&s| s as f64).collect();
    let auc = app.auc();
    let header = format!("{n_targets} targets   {n_decoys} decoys   AUC {auc:.4}\n{counts}");
    frame.render_widget(Paragraph::new(header), chunks[0]);

    let nan_count = score.iter().filter(|v| !v.is_finite()).count();
    let (hist, dropped) = column_hist(app, &score);
    let title = format!("discriminant_score [{}]", app.x().label());
    let subtitle = hist_subtitle(dropped, nan_count, hist.n_out);
    draw_hist(frame, chunks[1], &title, &subtitle, &hist, app.y());
}

/// Targets-passing-vs-q-value curve, plus a threshold table at fixed FDR
/// cutoffs.
fn draw_fdr(frame: &mut Frame, app: &mut App, area: Rect) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Percentage(60), Constraint::Percentage(40)])
        .split(area);

    let curve = app.qvalue_curve();
    let view = app.view();
    let n_targets = view.is_target.iter().filter(|&&t| t).count() as f64;
    let ymax = if n_targets > 0.0 { n_targets } else { 1.0 };

    let dataset = Dataset::default()
        .name("targets")
        .marker(symbols::Marker::Braille)
        .graph_type(GraphType::Line)
        .style(Style::default().fg(TARGET_COLOR))
        .data(&curve);
    let chart = Chart::new(vec![dataset])
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title("q-value curve"),
        )
        .x_axis(Axis::default().title("q-value").bounds([0.0, 1.0]))
        .y_axis(Axis::default().title("targets").bounds([0.0, ymax]));
    frame.render_widget(chart, chunks[0]);

    let thresholds =
        curves::threshold_table(&view.qvalue, &view.is_target, &[0.01, 0.05, 0.1, 0.5, 1.0]);
    let rows: Vec<Row> = thresholds
        .iter()
        .map(|&(q, total, targets, decoys)| {
            Row::new(vec![
                format!("{q:.2}"),
                total.to_string(),
                targets.to_string(),
                decoys.to_string(),
            ])
        })
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
    // The plot's axes are (decoy CDF, target CDF): good target/decoy
    // separation pulls the curve BELOW the y = x diagonal (targets need a
    // higher score threshold before as many of them have "arrived" as
    // decoys have), so the title carries that reading hint rather than
    // leaving it to the viewer to work out from two bare axis labels.
    let title = "Calibration (below y=x = separation)";

    let curve = app.pp_curve();
    if curve.is_empty() {
        frame.render_widget(
            Paragraph::new("PP plot needs both targets and decoys")
                .block(Block::default().borders(Borders::ALL).title(title)),
            area,
        );
        return;
    }

    let reference = [(0.0, 0.0), (1.0, 1.0)];
    let datasets = vec![
        Dataset::default()
            .name("PP")
            .marker(symbols::Marker::Braille)
            .graph_type(GraphType::Line)
            .style(Style::default().fg(TARGET_COLOR))
            .data(&curve),
        Dataset::default()
            .name("y = x")
            .marker(symbols::Marker::Braille)
            .graph_type(GraphType::Line)
            .style(Style::default().fg(Color::DarkGray))
            .data(&reference),
    ];
    let chart = Chart::new(datasets)
        .block(Block::default().borders(Borders::ALL).title(title))
        .x_axis(Axis::default().title("decoy CDF").bounds([0.0, 1.0]))
        .y_axis(Axis::default().title("target CDF").bounds([0.0, 1.0]));
    frame.render_widget(chart, area);
}

/// `{:.4}`, with `NaN` printed as `-` rather than the literal text `NaN`.
fn fmt4(v: f64) -> String {
    if v.is_nan() {
        "-".to_string()
    } else {
        format!("{v:.4}")
    }
}

/// Sortable feature table plus a histogram of the selected feature's column.
fn draw_features(frame: &mut Frame, app: &mut App, area: Rect) {
    let chunks = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(55), Constraint::Percentage(45)])
        .split(area);

    let visible = app.visible().to_vec();
    let total = app.view().n_features();

    let header = Row::new(vec![
        SortKey::Name.label(),
        SortKey::TargetMean.label(),
        SortKey::DecoyMean.label(),
        SortKey::Auc.label(),
        SortKey::CohensD.label(),
        SortKey::NanFrac.label(),
        SortKey::Gain.label(),
    ]);

    let mut rows = Vec::with_capacity(visible.len());
    for &j in &visible {
        let s = app.summary(j);
        let view = app.view();
        let name = view.feature_names[j].to_string();
        let gain = view.mean_gain(&view.feature_names[j]) as f64;
        rows.push(Row::new(vec![
            name,
            fmt4(s.target_mean),
            fmt4(s.decoy_mean),
            fmt4(s.auc),
            fmt4(s.cohens_d),
            fmt4(s.nan_frac),
            fmt4(gain),
        ]));
    }

    // ALL-lane feature names share long prefixes by construction (`ms1_*`,
    // `ms2_*`, `lazyscore_*`), so a name column fixed at the old 16 chars
    // renders an indistinguishable prefix for every row. `Min(16)` keeps that
    // as a floor (so a narrow terminal never does worse than before) but lets
    // the column grow to fill whatever space the numeric columns do not need
    // — e.g. 26 chars at a 160-column terminal, comfortably showing names
    // like `lazyscore_vs_baseline` in full. A plain `Fill(1)` or a larger
    // fixed `Min` would instead pull width away from the numeric columns
    // below (once the fixed name width plus six `Length(9)` columns exceeds
    // the panel's content width, the layout solver shrinks the `Length`
    // columns to make room, which is exactly the truncation this is fixing).
    // The numeric columns themselves are widened to fit `fmt4`'s widest
    // realistic output (e.g. `-307.0000`) without truncating.
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
        "Features ({}/{}) sort:{}",
        visible.len(),
        total,
        app.sort_key().label()
    );
    let table = Table::new(rows, widths)
        .header(header)
        .row_highlight_style(Style::default().add_modifier(Modifier::REVERSED))
        .block(Block::default().borders(Borders::ALL).title(title));

    frame.render_stateful_widget(table, chunks[0], app.table_state());

    match app.selected_feature() {
        Some(j) => {
            let values: Vec<f64> = app.view().feature_column(j).collect();
            let name = app.view().feature_names[j].to_string();
            let nan_count = values.iter().filter(|v| !v.is_finite()).count();
            let (hist, dropped) = column_hist(app, &values);
            let title = format!("{name} [{}]", app.x().label());
            let subtitle = hist_subtitle(dropped, nan_count, hist.n_out);
            draw_hist(frame, chunks[1], &title, &subtitle, &hist, app.y());
        }
        None => {
            frame.render_widget(
                Paragraph::new("no feature selected")
                    .block(Block::default().borders(Borders::ALL).title("Features")),
                chunks[1],
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::app::{
        App,
        Tab,
    };
    use crate::view::RescoreView;
    use ratatui::Terminal;
    use ratatui::backend::TestBackend;
    use ratatui::crossterm::event::{
        KeyCode,
        KeyEvent,
        KeyModifiers,
    };
    use std::sync::Arc;

    /// 40 rows x 2 features: targets score high, decoys low.
    fn fixture() -> (Vec<f64>, Vec<bool>, Vec<f32>, Vec<f32>) {
        let mut matrix = Vec::new();
        let mut is_target = Vec::new();
        let mut score = Vec::new();
        let mut qvalue = Vec::new();
        for i in 0..40 {
            let t = i % 2 == 0;
            matrix.push(if t { i as f64 } else { -(i as f64) });
            matrix.push(f64::NAN); // a wholly missing feature
            is_target.push(t);
            score.push(if t { i as f32 } else { -(i as f32) });
            qvalue.push(if t { 0.001 } else { 0.9 });
        }
        (matrix, is_target, score, qvalue)
    }

    fn render_tab(tab: Tab) -> String {
        let (matrix, is_target, score, qvalue) = fixture();
        let view = RescoreView {
            feature_names: vec![Arc::from("alpha_score"), Arc::from("all_nan_feature")],
            features: &matrix,
            is_target,
            score,
            qvalue,
            importance: Vec::new(),
        };
        let mut app = App::new(&view);
        while app.tab() != tab {
            app.handle_key(KeyEvent::new(KeyCode::Char('l'), KeyModifiers::NONE));
        }
        let backend = TestBackend::new(120, 40);
        let mut terminal = Terminal::new(backend).expect("test terminal");
        terminal.draw(|f| draw(f, &mut app)).expect("draw");
        let buffer = terminal.backend().buffer().clone();
        buffer
            .content()
            .iter()
            .map(|cell| cell.symbol())
            .collect::<String>()
    }

    #[test]
    fn every_tab_renders_with_its_title() {
        for tab in Tab::ALL {
            let text = render_tab(tab);
            assert!(
                text.contains(tab.title()),
                "{} missing from render",
                tab.title()
            );
        }
    }

    #[test]
    fn features_tab_shows_feature_names() {
        let text = render_tab(Tab::Features);
        assert!(text.contains("alpha_score"));
    }

    #[test]
    fn rendering_survives_an_empty_view() {
        let matrix: [f64; 0] = [];
        let view = RescoreView {
            feature_names: Vec::new(),
            features: &matrix,
            is_target: Vec::new(),
            score: Vec::new(),
            qvalue: Vec::new(),
            importance: Vec::new(),
        };
        // Every tab, not just the one `App::new` starts on: each tab's draw
        // function reads the view/cache independently, and an empty matrix
        // exercises the "no feature selected" / "needs both classes" paths
        // that a non-empty fixture never reaches.
        for tab in Tab::ALL {
            let mut app = App::new(&view);
            while app.tab() != tab {
                app.handle_key(KeyEvent::new(KeyCode::Char('l'), KeyModifiers::NONE));
            }
            let backend = TestBackend::new(80, 24);
            let mut terminal = Terminal::new(backend).expect("test terminal");
            terminal
                .draw(|f| draw(f, &mut app))
                .unwrap_or_else(|e| panic!("{:?} must still draw on an empty view: {e}", tab));
        }
    }

    /// 100x20 backend, 40 features, ~13-row viewport: reproduces the
    /// scenario where `TableState::default()` every frame pinned the
    /// selected row to the bottom edge instead of letting the viewport
    /// scroll. The offset must persist (and advance) across renders as the
    /// cursor moves, since `TableState::offset` is what makes a `Table`
    /// widget scroll at all.
    #[test]
    fn features_table_offset_persists_and_advances_across_renders() {
        let n_features = 40;
        let n_rows = 4;
        let feature_names: Vec<Arc<str>> = (0..n_features)
            .map(|i| Arc::from(format!("feature_{i:03}").as_str()))
            .collect();
        let mut matrix = Vec::new();
        for r in 0..n_rows {
            for j in 0..n_features {
                matrix.push((r * n_features + j) as f64);
            }
        }
        let is_target: Vec<bool> = (0..n_rows).map(|i| i % 2 == 0).collect();
        let score: Vec<f32> = (0..n_rows)
            .map(|i| if i % 2 == 0 { 1.0 } else { -1.0 })
            .collect();
        let qvalue: Vec<f32> = (0..n_rows)
            .map(|i| if i % 2 == 0 { 0.001 } else { 0.9 })
            .collect();
        let view = RescoreView {
            feature_names,
            features: &matrix,
            is_target,
            score,
            qvalue,
            importance: Vec::new(),
        };
        let mut app = App::new(&view);
        while app.tab() != Tab::Features {
            app.handle_key(KeyEvent::new(KeyCode::Char('l'), KeyModifiers::NONE));
        }

        let backend = TestBackend::new(100, 20);
        let mut terminal = Terminal::new(backend).expect("test terminal");
        terminal.draw(|f| draw(f, &mut app)).expect("first draw");
        let offset_before = app.table_state().offset();

        for _ in 0..30 {
            app.handle_key(KeyEvent::new(KeyCode::Char('j'), KeyModifiers::NONE));
        }
        terminal.draw(|f| draw(f, &mut app)).expect("second draw");
        let offset_after = app.table_state().offset();

        assert!(
            offset_after > offset_before,
            "offset must advance as the cursor scrolls past the viewport: {offset_before} -> {offset_after}"
        );
    }

    /// At 160x30 with real, ALL-lane-length feature names (21 and 19 chars —
    /// both would have collapsed to the same 16-char prefix under the old
    /// `Constraint::Length(16)` name column, `"lazyscore_vs_bas"` for the
    /// first), and a numeric value whose `fmt4` output is wider than the old
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
        // Feature 0: target mean -307.0 (one target row), decoy mean 307.0.
        let matrix = [-307.0, 1.0, 307.0, 2.0];
        let view = RescoreView {
            feature_names: vec![Arc::from(name_a), Arc::from(name_b)],
            features: &matrix,
            is_target: vec![true, false],
            score: vec![1.0, -1.0],
            qvalue: vec![0.001, 0.9],
            importance: Vec::new(),
        };
        let mut app = App::new(&view);
        while app.tab() != Tab::Features {
            app.handle_key(KeyEvent::new(KeyCode::Char('l'), KeyModifiers::NONE));
        }
        let backend = TestBackend::new(160, 30);
        let mut terminal = Terminal::new(backend).expect("test terminal");
        terminal.draw(|f| draw(f, &mut app)).expect("draw");
        let buffer = terminal.backend().buffer().clone();
        let text: String = buffer.content().iter().map(|c| c.symbol()).collect();

        assert!(
            text.contains(name_b),
            "name column truncated a real feature name that appears only in the table: {text}"
        );
        assert!(
            text.contains("-307.0000"),
            "numeric column truncated fmt4's output: {text}"
        );
    }
}
