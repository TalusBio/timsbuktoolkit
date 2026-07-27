//! Every string on screen, formatted once before the TUI opens; a frame borrows
//! these and formats nothing. Kept apart from [`crate::precompute`] so that
//! "what fits on a panel border" and "what the numbers are" stay separate
//! questions.

use crate::precompute::{
    FeatureColumn,
    N_FEATURE_COLUMNS,
    Slot,
};
use crate::transform::Axis;
use crate::view::{
    RescoreView,
    ThresholdRow,
};

/// How many of the caller's thresholds, from the front, the Overview header
/// summarizes. The FDR tab tabulates all of them.
const OVERVIEW_SHOWN: usize = 3;

/// One zoom level of the FDR curve, with its axes resolved at init.
///
/// The labels depend only on the zoom and the target count, both fixed once the
/// run is loaded, so a redraw borrows these strings.
pub(crate) struct QCurve {
    pub(crate) points: Vec<(f64, f64)>,
    /// Upper q-value bound of this view.
    pub(crate) zoom: f64,
    /// Targets passing at the loosest threshold shown, floored at 1 so the
    /// axis never degenerates to `[0, 0]`.
    pub(crate) ymax: f64,
    pub(crate) title: String,
    pub(crate) x_labels: [String; 3],
    pub(crate) y_labels: [String; 3],
}

impl QCurve {
    pub(crate) fn new(points: Vec<(f64, f64)>, zoom: f64) -> Self {
        // Non-decreasing in y, so the last point is the maximum. Scaling to it
        // rather than to the run's total target count is the other half of
        // zooming.
        let ymax = points.last().map_or(0.0, |&(_, y)| y).max(1.0);
        Self {
            points,
            zoom,
            ymax,
            // `{zoom}` and not a fixed precision: `Display` is the shortest
            // round-tripping form, so 1.0 prints "1" and 0.05 prints "0.05".
            title: format!("q-value curve (q <= {zoom})"),
            x_labels: [
                "0".to_string(),
                format!("{}", zoom / 2.0),
                format!("{zoom}"),
            ],
            y_labels: [
                "0".to_string(),
                fmt_count((ymax / 2.0) as usize),
                fmt_count(ymax as usize),
            ],
        }
    }
}

/// Every label the dashboard hands to the renderer.
///
/// `titles` and `subtitles` are indexed by the histogram store's layout, which
/// [`crate::precompute::Dashboard`] owns; its accessors do that indexing so the
/// layout keeps one definition.
pub(crate) struct Labels {
    /// `[column][transform]`.
    pub(crate) titles: Vec<String>,
    /// `[column][transform][clip]`.
    pub(crate) subtitles: Vec<String>,
    /// What everything on screen was computed over, drawn once across the full
    /// terminal width rather than in every panel's subtitle.
    pub(crate) basis: String,
    pub(crate) overview_header: String,
    /// `[q, total, targets, decoys]` per threshold row.
    pub(crate) fdr_rows: Vec<[String; 4]>,
    /// One row of seven per feature.
    pub(crate) cells: Vec<[String; N_FEATURE_COLUMNS]>,
}

impl Labels {
    pub(crate) fn build(
        view: &RescoreView<'_>,
        slots: &[Slot],
        values: &[[f64; N_FEATURE_COLUMNS]],
        score_auc: f64,
        n_sampled: usize,
    ) -> Self {
        let n_rows = view.n_rows();
        let n_targets = view.is_target.iter().filter(|&&t| t).count();
        Self {
            titles: build_titles(view),
            subtitles: build_subtitles(slots, n_sampled),
            basis: build_basis(n_sampled, n_rows),
            overview_header: build_overview_header(
                n_targets,
                n_rows - n_targets,
                score_auc,
                view.thresholds
                    .get(..OVERVIEW_SHOWN)
                    .unwrap_or(view.thresholds),
            ),
            fdr_rows: build_fdr_rows(view.thresholds),
            cells: build_cells(view, values),
        }
    }
}

/// One title per `(column, transform)`, the discriminant score last — it is
/// titled like a feature because it is stored and queried like one.
fn build_titles(view: &RescoreView<'_>) -> Vec<String> {
    view.feature_names
        .iter()
        .map(|n| &**n)
        .chain(std::iter::once("discriminant_score"))
        .flat_map(|name| Axis::ALL.map(|t| format!("{name} [{}]", t.label())))
        .collect()
}

/// One diagnostics line per stored histogram, as a fraction of the sample.
///
/// Two exclusions at two steps, neither a subset of the other: `dropped` is what
/// the transform refused, `outside range` what survived it and then fell outside
/// the axis. Kept short — the panel is under half the terminal width and its
/// border truncates the rest; run-level facts go to [`Labels::basis`].
fn build_subtitles(slots: &[Slot], n_sampled: usize) -> Vec<String> {
    let pct = |n: u32| {
        if n_sampled == 0 {
            0.0
        } else {
            100.0 * n as f64 / n_sampled as f64
        }
    };
    slots
        .iter()
        .map(|s| {
            if !s.plottable {
                return "nothing this transform can plot".to_string();
            }
            format!(
                "dropped {:.1}% by transform | {:.1}% outside range",
                pct(s.dropped),
                pct(s.n_out),
            )
        })
        .collect()
}

fn build_overview_header(
    n_targets: usize,
    n_decoys: usize,
    score_auc: f64,
    shown: &[ThresholdRow],
) -> String {
    let counts = shown
        .iter()
        .map(|r| {
            format!(
                "q<={:.0}%: {} (t{}/d{})",
                r.q * 100.0,
                r.total(),
                r.targets,
                r.decoys
            )
        })
        .collect::<Vec<_>>()
        .join("   ");
    format!("{n_targets} targets   {n_decoys} decoys   AUC {score_auc:.4}\n{counts}")
}

fn build_fdr_rows(thresholds: &[ThresholdRow]) -> Vec<[String; 4]> {
    thresholds
        .iter()
        .map(|r| {
            [
                format!("{:.2}", r.q),
                r.total().to_string(),
                r.targets.to_string(),
                r.decoys.to_string(),
            ]
        })
        .collect()
}

/// What every panel on screen was computed over, for the frame's top border.
fn build_basis(n_sampled: usize, n_rows: usize) -> String {
    if n_sampled < n_rows {
        format!(
            "{} rows, histograms from a {} sample",
            fmt_count(n_rows),
            fmt_count(n_sampled)
        )
    } else {
        format!("all {} rows", fmt_count(n_rows))
    }
}

fn build_cells(
    view: &RescoreView<'_>,
    values: &[[f64; N_FEATURE_COLUMNS]],
) -> Vec<[String; N_FEATURE_COLUMNS]> {
    values
        .iter()
        .enumerate()
        .map(|(j, row)| {
            std::array::from_fn(|k| match FeatureColumn::ALL[k] {
                FeatureColumn::Name => view.feature_names[j].to_string(),
                _ => fmt_stat(row[k]),
            })
        })
        .collect()
}

/// `{:.4}`, with `NaN` printed as `-` rather than the literal text `NaN`.
fn fmt_stat(v: f64) -> String {
    if v.is_nan() {
        "-".to_string()
    } else {
        format!("{v:.4}")
    }
}

fn fmt_count(n: usize) -> String {
    match n {
        0..=9_999 => n.to_string(),
        10_000..=999_999 => format!("{:.0}k", n as f64 / 1e3),
        _ => format!("{:.2}M", n as f64 / 1e6),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precompute::tests::fixture;
    use crate::transform::XTransform;

    /// A reader has to be able to tell an all-rows histogram from a sampled
    /// one. It is a fact about the run, so it lives on the dashboard rather
    /// than being repeated in every panel's subtitle.
    #[test]
    fn the_basis_says_whether_the_histograms_were_sampled() {
        let f = fixture(80, &[("a", &|i, _| i as f64)]);
        assert_eq!(f.build().basis(), "all 80 rows");

        let big = fixture(2_000, &[("a", &|i, _| i as f64)]);
        assert_eq!(
            big.build_with(500).basis(),
            "2000 rows, histograms from a 500 sample"
        );
        // `fmt_count`'s unit thresholds, on the path that uses them.
        let huge = fixture(10_500, &[("a", &|i, _| i as f64)]);
        assert_eq!(
            huge.build_with(500).basis(),
            "10k rows, histograms from a 500 sample"
        );
    }

    /// The subtitle has to fit the histogram panel, which is under half of a
    /// terminal's width; anything longer is silently cut off by the border it is
    /// drawn on. `log10` rejects the non-positive half, and the subtitle is the
    /// only place that says so — a histogram that silently shrank would read as
    /// "those rows do not exist".
    #[test]
    fn subtitles_fit_a_narrow_panel_and_report_what_a_transform_refused() {
        let dash = fixture(100, &[("half_negative", &|i, _| i as f64 - 50.0)]).build();
        for t in Axis::ALL {
            for clip in [false, true] {
                let s = dash.subtitle(0, t, clip);
                assert!(s.len() <= 60, "{t:?} clip={clip}: {} chars: {s}", s.len());
            }
        }
        let sub = dash.subtitle(0, Axis::Value(XTransform::Log10), true);
        assert!(sub.contains("dropped 51.0% by transform"), "{sub}");
    }

    #[test]
    fn table_cells_hold_the_formatted_column_values() {
        let dash = fixture(
            80,
            &[("alpha", &|i, _| i as f64), ("beta", &|_, _| f64::NAN)],
        )
        .build();
        assert_eq!(dash.cells().len(), 2);
        assert_eq!(dash.cells()[0][0], "alpha");
        assert_eq!(dash.cells()[1][3], "-", "a NaN AUC prints as a dash");
        assert_eq!(dash.cells()[1][5], "1.0000", "the NaN fraction is 1");
    }
}
