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
    TableState,
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
fn draw_hist(
    frame: &mut Frame,
    area: Rect,
    title: &str,
    hist: &Hist,
    y: transform::YTransform,
    x_label: &str,
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

    let mid = (hist.lo + hist.hi) / 2.0;
    let chart = Chart::new(datasets)
        .block(Block::default().borders(Borders::ALL).title(title))
        .x_axis(
            Axis::default()
                .title(x_label)
                .bounds([hist.lo, hist.hi])
                .labels([
                    format!("{:.3}", hist.lo),
                    format!("{mid:.3}"),
                    format!("{:.3}", hist.hi),
                ]),
        )
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

    let range = if clip {
        stats::percentile_range(transformed.iter().copied(), 0.5, 99.5)
    } else {
        stats::percentile_range(transformed.iter().copied(), 0.0, 100.0)
    };

    let hist = match range {
        Some((lo, hi)) => stats::histogram(transformed.into_iter(), &labels, lo, hi, N_BINS),
        None => stats::histogram(std::iter::empty(), &[], 0.0, 1.0, N_BINS),
    };
    (hist, dropped)
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
    let auc = stats::auc_exact(view.score.iter().map(|&s| s as f64), &view.is_target);
    let thresholds = curves::threshold_table(&view.qvalue, &view.is_target, &[0.01, 0.05, 0.10]);
    let counts = thresholds
        .iter()
        .map(|&(q, total, targets, decoys)| {
            format!("q<={:.0}%: {total} (t{targets}/d{decoys})", q * 100.0)
        })
        .collect::<Vec<_>>()
        .join("   ");
    let header = format!("{n_targets} targets   {n_decoys} decoys   AUC {auc:.4}\n{counts}");
    frame.render_widget(Paragraph::new(header), chunks[0]);

    let score: Vec<f64> = view.score.iter().map(|&s| s as f64).collect();
    let nan_count = score.iter().filter(|v| !v.is_finite()).count();
    let (hist, dropped) = column_hist(app, &score);
    let title = format!(
        "discriminant_score [{}] (dropped {dropped}, of which {nan_count} non-finite)",
        app.x().label()
    );
    draw_hist(frame, chunks[1], &title, &hist, app.y(), "score");
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
    let curve = app.pp_curve();
    if curve.is_empty() {
        frame.render_widget(
            Paragraph::new("PP plot needs both targets and decoys")
                .block(Block::default().borders(Borders::ALL).title("Calibration")),
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
        .block(Block::default().borders(Borders::ALL).title("Calibration"))
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
    let cursor = app.cursor();

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

    let widths = [
        Constraint::Length(16),
        Constraint::Length(6),
        Constraint::Length(6),
        Constraint::Length(4),
        Constraint::Length(4),
        Constraint::Length(4),
        Constraint::Length(5),
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

    let mut state = TableState::default();
    state.select(if visible.is_empty() {
        None
    } else {
        Some(cursor)
    });
    frame.render_stateful_widget(table, chunks[0], &mut state);

    match app.selected_feature() {
        Some(j) => {
            let values: Vec<f64> = app.view().feature_column(j).collect();
            let name = app.view().feature_names[j].to_string();
            let nan_count = values.iter().filter(|v| !v.is_finite()).count();
            let (hist, dropped) = column_hist(app, &values);
            let title = format!(
                "{name} [{}] (dropped {dropped}, of which {nan_count} non-finite)",
                app.x().label()
            );
            draw_hist(frame, chunks[1], &title, &hist, app.y(), &name);
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
        let mut app = App::new(&view);
        let backend = TestBackend::new(80, 24);
        let mut terminal = Terminal::new(backend).expect("test terminal");
        terminal
            .draw(|f| draw(f, &mut app))
            .expect("empty view must still draw");
    }
}
