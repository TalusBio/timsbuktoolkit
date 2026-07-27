//! Everything the dashboard shows, computed once before the TUI opens.
//!
//! The matrix is touched exactly twice and is then dead:
//!
//! **Pass A** is one row-major sweep with every column's [`ColumnStats`]
//! accumulator live. It is exact over all rows, O(1) per value, and
//! bandwidth-bound. It is the only thing here that is linear in rows.
//!
//! **Pass B** gathers `min(n_rows, `[`DEFAULT_SAMPLE`]`)` rows *by random
//! index* and, per column, sorts once. That one sort answers three questions
//! that would otherwise each need their own: the Mann-Whitney AUC (a
//! mid-rank walk), every transform's clip range (a `partition_point` for the
//! domain suffix, then two order statistics), and `RankPercentile` (sorted
//! position *is* the percentile). Binning is order-independent, so it reads the
//! sorted array directly. Pass B is flat in rows.
//!
//! What comes out is a lookup table: `(column, transform, clip)` indexes a
//! stored 512-bin histogram with a stored axis range and a pre-formatted title
//! and subtitle. A frame computes nothing over the data. It still allocates:
//! ratatui's own widget model, plus three axis-label strings per plotted panel.
//!
//! Random indices are load-bearing. The rescoring pipeline hands the dashboard
//! rows sorted descending by main score, so a contiguous or strided sample
//! would be all high-scoring targets.
//!
//! # What is exact and what is sampled
//!
//! Exact over all rows: every [`ColumnStats`] field (means, variances,
//! min/max, the positive and absolute floors, zero and non-finite counts), the
//! unclipped axis ranges derived from them, the discriminant score's AUC, and
//! all four curve/threshold tables.
//!
//! From the sample: every histogram's bin counts, every clipped axis range, and
//! the per-feature AUC. Sample size is the only accuracy knob; the
//! `sample_size_holds_its_accuracy_claim` test measures the error it buys.

use crate::stats::{
    self,
    ColumnStats,
    HistView,
    N_BINS,
};
use crate::transform::{
    Axis,
    XTransform,
};
use crate::view::{
    RescoreView,
    ViewError,
};
use crate::{
    curves,
    cycle,
};
use rayon::prelude::*;
use std::sync::Arc;

/// Rows drawn for pass B.
///
/// Accuracy improves as 1/sqrt(m) while pass B's cost grows linearly in it, so
/// past 250k the accuracy stops paying for the time. Measured at 1M x 131:
/// 100k / 250k / 500k give worst-case KS distance .015 / .010 / .009 and AUC
/// error .0052 / .0036 / .0026.
pub const DEFAULT_SAMPLE: usize = 250_000;

/// Percentiles bounding a clipped axis. The whole point is to keep a lone
/// outlier from stretching the range until the bulk collapses into one bin.
const CLIP_LO: f64 = 0.5;
const CLIP_HI: f64 = 99.5;

const N_AXES: usize = Axis::ALL.len();
/// Unclipped and clipped, in that order.
const N_CLIPS: usize = 2;
/// `target` bins then `decoy` bins.
const BINS_PER_SLOT: usize = 2 * N_BINS;
const SLOTS_PER_COLUMN: usize = N_AXES * N_CLIPS;

/// How many of the caller's thresholds, from the front, the Overview header
/// summarizes. The FDR tab tabulates all of them.
const OVERVIEW_SHOWN: usize = 3;
const Q_CURVE_POINTS: usize = 100;
const PP_CURVE_POINTS: usize = 200;

/// Upper q-value bounds the FDR curve can be drawn over, widest first.
///
/// Every one of these gets its own `Q_CURVE_POINTS` grid — see
/// [`curves::qvalue_curves`] for why a zoom cannot just be a slice of the wide
/// curve — which is affordable because they share one sort.
pub const Q_ZOOMS: [f64; 4] = [1.0, 0.1, 0.05, 0.01];

/// Which of [`Q_ZOOMS`] the FDR tab opens on.
///
/// Not the widest. `q <= 1` is the view that says the least: past the FDR
/// cutoffs anyone reports, the curve is flat, so the default is the tightest
/// range that still contains all of the usual reporting thresholds.
pub const DEFAULT_Q_ZOOM: usize = 1;

/// One stored histogram's axis and totals. The counts themselves live in a
/// single flat `Vec<u32>` so the whole store is one allocation.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Slot {
    lo: f64,
    hi: f64,
    n_target: u32,
    n_decoy: u32,
    /// Sampled values this transform refused (non-finite, or out of domain).
    dropped: u32,
    /// Sampled survivors that landed outside `[lo, hi]`.
    n_out: u32,
    /// False when the column has no value this transform can plot, in which
    /// case `lo`/`hi` are placeholders and every count is zero.
    plottable: bool,
}

impl Slot {
    const EMPTY: Self = Self {
        lo: 0.0,
        hi: 1.0,
        n_target: 0,
        n_decoy: 0,
        dropped: 0,
        n_out: 0,
        plottable: false,
    };
}

/// A column of the features table.
///
/// [`Self::ALL`] is the table's column order, the cycle `s` steps through, and
/// the order cells are built in — one definition, so the header, the numbers
/// under it and the sort key cannot disagree about what column 4 is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FeatureColumn {
    Name,
    TargetMean,
    DecoyMean,
    Auc,
    CohensD,
    NanFrac,
    Gain,
}

pub const N_FEATURE_COLUMNS: usize = FeatureColumn::ALL.len();

impl FeatureColumn {
    pub const ALL: [FeatureColumn; 7] = [
        Self::Name,
        Self::TargetMean,
        Self::DecoyMean,
        Self::Auc,
        Self::CohensD,
        Self::NanFrac,
        Self::Gain,
    ];

    pub fn next(self) -> Self {
        cycle::step(&Self::ALL, self, 1)
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Name => "name",
            Self::TargetMean => "target mean",
            Self::DecoyMean => "decoy mean",
            Self::Auc => "AUC",
            Self::CohensD => "|d|",
            Self::NanFrac => "NaN%",
            Self::Gain => "gain",
        }
    }
}

/// One zoom level of the FDR curve, with its axes resolved at init.
///
/// The labels are formatted here rather than in the panel because they depend
/// only on the zoom and the target count, both of which are fixed once the run
/// is loaded; a redraw borrows these strings.
pub struct QCurve {
    pub points: Vec<(f64, f64)>,
    /// Upper q-value bound of this view.
    pub zoom: f64,
    /// Targets passing at the loosest threshold shown, floored at 1 so the
    /// axis never degenerates to `[0, 0]`.
    pub ymax: f64,
    pub title: String,
    pub x_labels: [String; 3],
    pub y_labels: [String; 3],
}

impl QCurve {
    fn new(points: Vec<(f64, f64)>, zoom: f64) -> Self {
        // The curve is non-decreasing in y, so its last point is its maximum.
        // Scaling to that rather than to the run's total target count is the
        // other half of zooming: a `q <= 0.01` panel scaled to the q = 1 count
        // would still draw the interesting part flat against the axis.
        let ymax = points.last().map_or(0.0, |&(_, y)| y).max(1.0);
        Self {
            points,
            zoom,
            ymax,
            // `{zoom}` and not a fixed precision: f64's `Display` is the
            // shortest round-tripping form, so 1.0 prints "1" and 0.05 prints
            // "0.05" without a trailing-zero rule per level.
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

/// The materialized dashboard. Owns everything on screen; borrows nothing.
pub struct Dashboard {
    pub(crate) feature_names: Vec<Arc<str>>,

    /// Pass A output, one per column (features, then the score).
    pub(crate) stats: Vec<ColumnStats>,
    /// Per-feature Mann-Whitney AUC over the sample; the score's entry is
    /// replaced by the exact all-rows value.
    pub(crate) auc: Vec<f64>,
    pub(crate) gain: Vec<f32>,

    /// `[column][transform][clip]` -> `2 * N_BINS` counts, target then decoy.
    counts: Vec<u32>,
    /// Aligned with `counts`' slot index.
    slots: Vec<Slot>,
    /// `[column][transform]` panel titles.
    titles: Vec<String>,
    /// `[column][transform][clip]` diagnostics lines.
    subtitles: Vec<String>,
    /// What everything on screen was computed over, drawn once across the full
    /// terminal width rather than repeated in every panel's subtitle.
    basis: String,

    /// Exact, all-rows AUC of the discriminant score. The headline separation
    /// number, so it is not sampled: it is one column and one sort, run
    /// alongside pass A rather than after it.

    /// The FDR curve at every zoom in [`Q_ZOOMS`], same order.
    q_curves: Vec<QCurve>,
    pub(crate) pp_curve: Vec<(f64, f64)>,
    /// Overview header, formatted once.
    pub(crate) overview_header: String,
    /// FDR threshold table cells, `[q, total, targets, decoys]` per row.
    pub(crate) fdr_rows: Vec<[String; 4]>,

    /// Pre-formatted feature-table cells, one row of seven per feature. Table
    /// rows borrow these, so a frame formats nothing.
    pub(crate) cells: Vec<[String; N_FEATURE_COLUMNS]>,
}

impl Dashboard {
    /// Materialize everything the TUI shows.
    ///
    /// `sample` is the pass-B row count — [`DEFAULT_SAMPLE`] unless the caller
    /// has a reason. At or above `view.n_rows()` every row is taken in order,
    /// which makes the histograms, clip ranges and per-feature AUCs exact
    /// rather than sampled.
    pub fn build(view: &RescoreView<'_>, sample: usize) -> Result<Self, ViewError> {
        view.validate()?;
        let n_features = view.n_features();
        let n_cols = n_features + 1;
        let n_rows = view.n_rows();

        // Pass A and the all-rows curve work are independent, and both are
        // internally parallel; rayon nests them without oversubscribing.
        let t0 = std::time::Instant::now();
        let (stats, all_rows) = rayon::join(|| pass_a(view), || all_rows_tables(view));
        let t_pass_a = t0.elapsed();

        let t0 = std::time::Instant::now();
        let (sample_matrix, sample_labels) = gather_sample(view, sample);
        let n_sampled = sample_labels.len();

        let mut counts = vec![0u32; n_cols * SLOTS_PER_COLUMN * BINS_PER_SLOT];
        let mut slots = vec![Slot::EMPTY; n_cols * SLOTS_PER_COLUMN];
        let mut auc = vec![f64::NAN; n_cols];

        counts
            .par_chunks_mut(SLOTS_PER_COLUMN * BINS_PER_SLOT)
            .zip(slots.par_chunks_mut(SLOTS_PER_COLUMN))
            .zip(auc.par_iter_mut())
            .enumerate()
            .for_each(|(j, ((count_slab, slot_slab), auc_j))| {
                let mut buf = ColumnBuf::with_capacity(n_sampled);
                // Strided straight into the sort buffer; the column is never
                // materialized on its own.
                let column = sample_matrix.iter().skip(j).step_by(n_cols).copied();
                *auc_j = precompute_column(
                    column,
                    &sample_labels,
                    &stats[j],
                    slot_slab,
                    count_slab,
                    &mut buf,
                );
            });

        // The score is the column past the last feature, and its AUC is the
        // exact all-rows one rather than the sampled sweep's.
        auc[n_features] = all_rows.score_auc;

        // Pass A is linear in rows and pass B is flat, so which one dominates
        // is a property of the run, not a constant. Logged so that stays
        // visible without a profiler.
        tracing::debug!(
            "rescore dashboard precompute: pass A (all {n_rows} rows) {:.3} s, \
             pass B ({n_sampled} sampled) {:.3} s",
            t_pass_a.as_secs_f64(),
            t0.elapsed().as_secs_f64(),
        );

        let n_targets = view.is_target.iter().filter(|&&t| t).count();
        let names: Vec<&str> = view
            .feature_names
            .iter()
            .map(|n| &**n)
            .chain(std::iter::once("discriminant_score"))
            .collect();
        let titles = build_titles(&names);
        let subtitles = build_subtitles(&slots, n_sampled);
        let basis = build_basis(n_sampled, n_rows);
        let cells = build_cells(view, &stats, &auc);

        Ok(Self {
            feature_names: view.feature_names.to_vec(),
            stats,
            auc,
            gain: view.gain.to_vec(),
            counts,
            slots,
            titles,
            subtitles,
            basis,
            q_curves: all_rows
                .qvalue_curves
                .into_iter()
                .zip(Q_ZOOMS)
                .map(|(points, zoom)| QCurve::new(points, zoom))
                .collect(),
            pp_curve: all_rows.pp_curve,
            overview_header: build_overview_header(
                n_targets,
                n_rows - n_targets,
                all_rows.score_auc,
                view.thresholds
                    .get(..OVERVIEW_SHOWN)
                    .unwrap_or(view.thresholds),
            ),
            fdr_rows: build_fdr_rows(view.thresholds),
            cells,
        })
    }

    pub fn n_features(&self) -> usize {
        self.feature_names.len()
    }

    /// What everything on screen was computed over, e.g. `"2.01M rows,
    /// histograms from a 250k sample"`.
    pub fn basis(&self) -> &str {
        &self.basis
    }

    /// How many zoom levels the FDR curve offers.
    pub fn n_q_zooms(&self) -> usize {
        self.q_curves.len()
    }

    /// The FDR curve at zoom level `i`, wrapping so the caller can hold a bare
    /// counter and never bounds-check.
    pub fn q_curve(&self, i: usize) -> &QCurve {
        &self.q_curves[i % self.q_curves.len()]
    }

    /// The discriminant score's column. It is stored and queried exactly like
    /// a feature, one past the last of them.
    pub fn score_column(&self) -> usize {
        self.n_features()
    }

    fn slot_index(&self, column: usize, t: Axis, clip: bool) -> usize {
        (column * N_AXES + t.index()) * N_CLIPS + usize::from(clip)
    }

    /// The stored histogram for a column. Pure indexing — this is the whole of
    /// what a redraw does.
    pub fn hist(&self, column: usize, t: Axis, clip: bool) -> HistView<'_> {
        let slot = self.slot_index(column, t, clip);
        let s = &self.slots[slot];
        let base = slot * BINS_PER_SLOT;
        let (target, decoy) = self.counts[base..base + BINS_PER_SLOT].split_at(N_BINS);
        HistView {
            lo: s.lo,
            hi: s.hi,
            target,
            decoy,
            n_target: s.n_target,
            n_decoy: s.n_decoy,
        }
    }

    pub fn title(&self, column: usize, t: Axis) -> &str {
        &self.titles[column * N_AXES + t.index()]
    }

    pub fn subtitle(&self, column: usize, t: Axis, clip: bool) -> &str {
        &self.subtitles[self.slot_index(column, t, clip)]
    }

    /// A column's AUC: sampled for a feature, exact for the score.
    pub fn auc(&self, column: usize) -> f64 {
        self.auc[column]
    }

    /// Feature `j`'s value under `col`, or `None` for the name column.
    pub fn feature_value(&self, j: usize, col: FeatureColumn) -> Option<f64> {
        column_value(col, &self.stats[j], self.auc[j], self.gain[j])
    }
}

/// Pass A: one row-major sweep per thread over a contiguous row range, then a
/// Chan-combination merge.
///
/// Row-major rather than one strided sweep per column. A strided read of
/// `feature_column(j)` pulls a 64-byte line for every 8-byte value and discards
/// 7/8 of it; this reads every line once and uses all of it, with the whole
/// accumulator block (~11 KB at 131 features) resident the whole time.
fn pass_a(view: &RescoreView<'_>) -> Vec<ColumnStats> {
    let n_features = view.n_features();
    let n_cols = n_features + 1;
    let n_rows = view.n_rows();
    let threads = rayon::current_num_threads().max(1);
    let chunk_rows = n_rows.div_ceil(threads).max(1);

    view.features
        .par_chunks(chunk_rows * n_features)
        .zip(view.is_target.par_chunks(chunk_rows))
        .zip(view.score.par_chunks(chunk_rows))
        .map(|((rows, labels), scores)| {
            let mut acc = vec![ColumnStats::IDENTITY; n_cols];
            for ((row, &is_t), &s) in rows.chunks_exact(n_features).zip(labels).zip(scores) {
                for (a, &v) in acc.iter_mut().zip(row) {
                    a.push(v, is_t);
                }
                acc[n_features].push(s as f64, is_t);
            }
            acc
        })
        .reduce(
            || vec![ColumnStats::IDENTITY; n_cols],
            |mut a, b| {
                for (x, y) in a.iter_mut().zip(&b) {
                    x.merge(y);
                }
                a
            },
        )
}

/// The all-rows tables plus the exact score AUC. Each scans or sorts every row
/// and none depends on the others, so they run concurrently.
struct AllRows {
    qvalue_curves: Vec<Vec<(f64, f64)>>,
    pp_curve: Vec<(f64, f64)>,
    score_auc: f64,
}

fn all_rows_tables(view: &RescoreView<'_>) -> AllRows {
    let ((qvalue_curves, pp_curve), score_auc) = rayon::join(
        || {
            rayon::join(
                || curves::qvalue_curves(view.qvalue, view.is_target, Q_CURVE_POINTS, &Q_ZOOMS),
                || curves::pp_curve(view.score, view.is_target, PP_CURVE_POINTS),
            )
        },
        || stats::auc_exact(view.score.iter().map(|&s| s as f64), view.is_target),
    );
    AllRows {
        qvalue_curves,
        pp_curve,
        score_auc,
    }
}

/// Deterministic mixed LCG. Two runs over the same data must show the same
/// histograms — a sample that moved between runs would be indistinguishable
/// from a bug — and this crate carries no `rand` dependency.
struct Lcg(u64);

impl Lcg {
    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        // An LCG's low bits are weak, so mix before use.
        let mut z = self.0;
        z = (z ^ (z >> 33)).wrapping_mul(0xff51_afd7_ed55_8ccd);
        z = (z ^ (z >> 33)).wrapping_mul(0xc4ce_b9fe_1a85_ec53);
        z ^ (z >> 33)
    }

    /// Uniform in `0..n`, by multiply-shift rather than modulo, so there is no
    /// bias toward low indices.
    fn below(&mut self, n: usize) -> usize {
        ((self.next_u64() as u128 * n as u128) >> 64) as usize
    }
}

/// Gather `min(sample, n_rows)` rows into a compact `m x (n_features + 1)` matrix: the
/// feature row, then the discriminant score as one more column.
///
/// Rows are drawn with replacement by random index. `sample >= n_rows` takes
/// every row in order instead, which makes pass B exact.
fn gather_sample(view: &RescoreView<'_>, sample: usize) -> (Vec<f64>, Vec<bool>) {
    let n_rows = view.n_rows();
    let n_features = view.n_features();
    let m = sample.min(n_rows);
    let mut rng = Lcg(0xC0FF_EE00_1234_5678);
    let mut matrix = Vec::with_capacity(m * (n_features + 1));
    let mut labels = Vec::with_capacity(m);
    for k in 0..m {
        let i = if m == n_rows { k } else { rng.below(n_rows) };
        matrix.extend_from_slice(view.row(i));
        matrix.push(view.score[i] as f64);
        labels.push(view.is_target[i]);
    }
    (matrix, labels)
}

/// Per-column scratch, reused across every transform of one column.
struct ColumnBuf {
    /// Finite `(value, is_target)` pairs, sorted ascending. The one sort.
    pairs: Vec<(f64, bool)>,
    /// Mid-rank percentile of each entry of `pairs`, i.e. the whole
    /// `RankPercentile` transform, for the cost of one linear pass.
    pct: Vec<f64>,
}

impl ColumnBuf {
    fn with_capacity(m: usize) -> Self {
        Self {
            pairs: Vec::with_capacity(m),
            pct: Vec::with_capacity(m),
        }
    }
}

/// Accumulator for one `(transform, clip)` histogram while the sorted column is
/// walked once for all of them.
struct Fill<'a> {
    lo: f64,
    hi: f64,
    span: f64,
    target: &'a mut [u32],
    decoy: &'a mut [u32],
    n_out: u32,
    plottable: bool,
}

impl Fill<'_> {
    #[inline]
    fn push(&mut self, y: f64, is_target: bool) {
        if !self.plottable {
            return;
        }
        match stats::bin_index(y, self.lo, self.hi, self.span) {
            Some(b) => {
                if is_target {
                    self.target[b] += 1;
                } else {
                    self.decoy[b] += 1;
                }
            }
            None => self.n_out += 1,
        }
    }
}

/// Fill one column's `N_AXES * N_CLIPS` histograms and return its
/// Mann-Whitney AUC. One sort, then one walk per transform.
fn precompute_column(
    values: impl Iterator<Item = f64>,
    labels: &[bool],
    exact: &ColumnStats,
    slots: &mut [Slot],
    counts: &mut [u32],
    buf: &mut ColumnBuf,
) -> f64 {
    buf.pairs.clear();
    buf.pairs.extend(
        values
            .zip(labels)
            .filter(|(v, _)| v.is_finite())
            .map(|(v, &t)| (v, t)),
    );
    // `total_cmp` rather than `partial_cmp(..).unwrap_or(Equal)`: NaNs are
    // already filtered out, and a total order lets `sort_unstable_by` keep its
    // guarantees instead of relying on a comparator that could be inconsistent.
    buf.pairs.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    let auc = stats::auc_from_sorted(&buf.pairs);
    stats::mid_rank_percentiles(&buf.pairs, &mut buf.pct);
    let n_non_finite = (labels.len() - buf.pairs.len()) as u32;

    for (axis_i, &t) in Axis::ALL.iter().enumerate() {
        // Survivors of a domain restriction are a suffix of a sorted column.
        let start = buf.pairs.partition_point(|p| !t.accepts(p.0));
        let survivors = &buf.pairs[start..];
        let pct = &buf.pct[start..];
        let ranges = [exact_range(t, exact), clip_range(t, survivors, pct)];

        let slab_base = axis_i * N_CLIPS * BINS_PER_SLOT;
        let (unclipped, clipped) =
            counts[slab_base..slab_base + N_CLIPS * BINS_PER_SLOT].split_at_mut(BINS_PER_SLOT);
        let (ut, ud) = unclipped.split_at_mut(N_BINS);
        let (ct, cd) = clipped.split_at_mut(N_BINS);

        let mut fills = [new_fill(ranges[0], ut, ud), new_fill(ranges[1], ct, cd)];

        let mut dropped = n_non_finite + start as u32;
        for (k, &(v, is_t)) in survivors.iter().enumerate() {
            let y = match t {
                // Sorted position, computed for the whole column above.
                Axis::RankPercentile => pct[k],
                Axis::Value(t) => match t.apply(v) {
                    Some(y) => y,
                    // Only reachable when the map overflows a finite input
                    // (`square` of ~1e154). Counted as a drop, like any other
                    // value the transform cannot place.
                    None => {
                        dropped += 1;
                        continue;
                    }
                },
            };
            for f in &mut fills {
                f.push(y, is_t);
            }
        }

        for (clip_i, f) in fills.iter().enumerate() {
            let slot = axis_i * N_CLIPS + clip_i;
            slots[slot] = Slot {
                lo: f.lo,
                hi: f.hi,
                n_target: f.target.iter().sum(),
                n_decoy: f.decoy.iter().sum(),
                dropped,
                n_out: f.n_out,
                plottable: f.plottable,
            };
        }
    }
    auc
}

fn new_fill<'a>(
    range: Option<(f64, f64)>,
    target: &'a mut [u32],
    decoy: &'a mut [u32],
) -> Fill<'a> {
    match range {
        Some((lo, hi)) => Fill {
            lo,
            hi,
            span: hi - lo,
            target,
            decoy,
            n_out: 0,
            plottable: true,
        },
        None => Fill {
            lo: Slot::EMPTY.lo,
            hi: Slot::EMPTY.hi,
            span: 1.0,
            target,
            decoy,
            n_out: 0,
            plottable: false,
        },
    }
}

/// Clipped axis bounds, from the sample's own percentiles.
///
/// This is why the outlier failure cannot recur: a lone 1e6 row is outside
/// `p99.5` by construction, so it cannot stretch the range until the
/// bulk collapses into bin 0.
fn clip_range(t: Axis, survivors: &[(f64, bool)], pct: &[f64]) -> Option<(f64, f64)> {
    let k = survivors.len();
    if k == 0 {
        return None;
    }
    let rank = |p: f64| (((p / 100.0) * (k - 1) as f64).round() as usize).min(k - 1);
    let (r_lo, r_hi) = (rank(CLIP_LO), rank(CLIP_HI));

    let (lo, hi) = match t {
        Axis::RankPercentile => (pct[r_lo], pct[r_hi]),
        // Monotone: the p-th percentile of `T(x)` is `T` of the p-th
        // percentile of `x`, so two array reads answer it.
        Axis::Value(t) if t.is_monotone() => {
            (t.apply(survivors[r_lo].0)?, t.apply(survivors[r_hi].0)?)
        }
        // `Square` decreases below zero, so its order statistics are those of
        // `|v|`, which one outward walk from the sign change produces without
        // a second sort.
        Axis::Value(t) => {
            let (x_lo, x_hi) = abs_order_stats(survivors, r_lo, r_hi);
            (t.apply(x_lo)?, t.apply(x_hi)?)
        }
    };
    Some(widen(lo, hi))
}

/// Unclipped axis bounds, from pass A's exact all-rows values.
///
/// Deliberately not sample-derived: this is the view that is *supposed* to show
/// the spike, so its endpoints have to be the real ones. Each transform's lower
/// bound is the smallest all-rows value it accepts, which is why pass A tracks
/// the positive and absolute floors rather than only min/max.
fn exact_range(t: Axis, c: &ColumnStats) -> Option<(f64, f64)> {
    // Checked before the `RankPercentile` arm: a wholly non-finite column has
    // no rank axis either, and handing one back would mark the slot plottable
    // when there is nothing to plot.
    if !c.has_finite() {
        return None;
    }
    let t = match t {
        // Bounded by construction, whatever else the column holds.
        Axis::RankPercentile => return Some((0.0, 100.0)),
        Axis::Value(t) => t,
    };
    let (x_lo, x_hi) = match t {
        XTransform::Log10 => (c.min_pos, c.hi),
        XTransform::Sqrt => (if c.n_zero > 0 { 0.0 } else { c.min_pos }, c.hi),
        XTransform::Square => (c.min_abs, c.lo.abs().max(c.hi.abs())),
        _ => (c.lo, c.hi),
    };
    // An infinite floor means the transform accepts nothing in this column.
    if !x_lo.is_finite() || !x_hi.is_finite() || x_lo > x_hi {
        return None;
    }
    Some(widen(t.apply(x_lo)?, t.apply(x_hi)?))
}

/// A degenerate range would collapse the plot onto one edge; give it unit width
/// instead, which is what a single-distinct-value column needs.
fn widen(a: f64, b: f64) -> (f64, f64) {
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    if hi > lo { (lo, hi) } else { (lo, lo + 1.0) }
}

/// The `r_lo`-th and `r_hi`-th smallest `|v|` over an ascending slice.
///
/// One outward walk from the sign change: negatives run left in increasing
/// magnitude, non-negatives run right, and merging the two enumerates `|v|` in
/// ascending order in `O(r_hi)` without allocating or sorting again.
fn abs_order_stats(sorted: &[(f64, bool)], r_lo: usize, r_hi: usize) -> (f64, f64) {
    debug_assert!(r_lo <= r_hi && r_hi < sorted.len());
    let split = sorted.partition_point(|p| p.0 < 0.0);
    let (mut up, mut down) = (split, split);
    let (mut lo, mut hi) = (0.0f64, 0.0f64);
    for r in 0..=r_hi {
        let a = sorted.get(up).map_or(f64::INFINITY, |p| p.0.abs());
        let b = if down > 0 {
            sorted[down - 1].0.abs()
        } else {
            f64::INFINITY
        };
        let v = if a <= b {
            up += 1;
            a
        } else {
            down -= 1;
            b
        };
        if r == r_lo {
            lo = v;
        }
        if r == r_hi {
            hi = v;
        }
    }
    (lo, hi)
}

fn build_titles(names: &[&str]) -> Vec<String> {
    let mut out = Vec::with_capacity(names.len() * N_AXES);
    for name in names {
        for t in Axis::ALL {
            out.push(format!("{name} [{}]", t.label()));
        }
    }
    out
}

/// One diagnostics line per stored histogram, as a fraction of the sample.
///
/// Two different exclusions at two different steps, so they are worded as
/// separate clauses rather than nested in one parenthetical: `dropped` counts
/// values the transform refused (including non-finite ones), while `outside
/// range` counts values that survived the transform and were then excluded at
/// binning time. Neither is a subset of the other.
///
/// Kept short on purpose. The histogram panel is under half the terminal width,
/// so anything longer is truncated by the border it is drawn on. Facts about
/// the run as a whole — the row count, whether the histograms are sampled —
/// belong to [`Dashboard::basis`], which is drawn once across the full width;
/// the non-finite fraction is already the table's `NaN%` column.
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
    shown: &[curves::ThresholdRow],
) -> String {
    let counts = shown
        .iter()
        .map(|r| {
            format!(
                "q<={:.0}%: {} (t{}/d{})",
                r.q * 100.0,
                r.total,
                r.targets,
                r.decoys
            )
        })
        .collect::<Vec<_>>()
        .join("   ");
    format!("{n_targets} targets   {n_decoys} decoys   AUC {score_auc:.4}\n{counts}")
}

fn build_fdr_rows(thresholds: &[curves::ThresholdRow]) -> Vec<[String; 4]> {
    thresholds
        .iter()
        .map(|r| {
            [
                format!("{:.2}", r.q),
                r.total.to_string(),
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

/// The feature table's cells, formatted once. Table rows borrow these `String`s
/// so a frame allocates none.
fn build_cells(
    view: &RescoreView<'_>,
    stats: &[ColumnStats],
    auc: &[f64],
) -> Vec<[String; N_FEATURE_COLUMNS]> {
    (0..view.n_features())
        .map(|j| {
            FeatureColumn::ALL.map(|col| {
                match column_value(col, &stats[j], auc[j], view.gain[j]) {
                    Some(v) => fmt_stat(v),
                    // The name column, the one that is not a number.
                    None => view.feature_names[j].to_string(),
                }
            })
        })
        .collect()
}

/// The one place a [`FeatureColumn`] turns into a number. `None` is the name
/// column, which has no numeric value and sorts as text.
fn column_value(col: FeatureColumn, c: &ColumnStats, auc: f64, gain: f32) -> Option<f64> {
    Some(match col {
        FeatureColumn::Name => return None,
        FeatureColumn::TargetMean => c.mean(true),
        FeatureColumn::DecoyMean => c.mean(false),
        FeatureColumn::Auc => auc,
        FeatureColumn::CohensD => c.cohens_d(),
        FeatureColumn::NanFrac => c.nan_frac(),
        FeatureColumn::Gain => gain as f64,
    })
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
pub(crate) mod tests {
    use super::*;

    /// AUC read off a value-binned histogram — the shape that produced the
    /// bug the fixtures below cover. Test-only, and deliberately not reachable
    /// from anything that ships.
    fn binned_auc(values: &[f64], is_target: &[bool]) -> f64 {
        let (lo, hi) = values
            .iter()
            .copied()
            .filter(|v| v.is_finite())
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), v| {
                (l.min(v), h.max(v))
            });
        let (mut t, mut d) = (vec![0.0f64; N_BINS], vec![0.0f64; N_BINS]);
        for (i, &v) in values.iter().enumerate() {
            if let Some(b) = stats::bin_index(v, lo, hi, hi - lo) {
                let side = if is_target[i] { &mut t } else { &mut d };
                side[b] += 1.0;
            }
        }
        let (nt, nd) = (t.iter().sum::<f64>(), d.iter().sum::<f64>());
        let (mut below, mut acc) = (0.0f64, 0.0f64);
        for i in 0..N_BINS {
            acc += t[i] * (below + 0.5 * d[i]);
            below += d[i];
        }
        acc / (nt * nd)
    }

    /// A synthetic run. THE definition of what test rows look like — the
    /// alternating labels, the perfectly-separating score and the
    /// target/decoy q-values are asserted against from two modules, so they
    /// get one home rather than a copy in each.
    pub(crate) struct Fixture {
        names: Vec<Arc<str>>,
        matrix: Vec<f64>,
        is_target: Vec<bool>,
        score: Vec<f32>,
        qvalue: Vec<f32>,
        gain: Vec<f32>,
        thresholds: Vec<curves::ThresholdRow>,
    }

    impl Fixture {
        pub(crate) fn view(&self) -> RescoreView<'_> {
            RescoreView {
                feature_names: &self.names,
                features: &self.matrix,
                is_target: &self.is_target,
                score: &self.score,
                qvalue: &self.qvalue,
                thresholds: &self.thresholds,
                gain: &self.gain,
            }
        }

        pub(crate) fn build(&self) -> Dashboard {
            Dashboard::build(&self.view(), DEFAULT_SAMPLE).expect("well-formed fixture")
        }
    }

    /// A named column, generated from `(row index, is_target)`.
    type Column<'a> = (&'a str, &'a dyn Fn(usize, bool) -> f64);

    /// `n_rows` rows of `columns.len()` features.
    fn fixture(n_rows: usize, columns: &[Column<'_>]) -> Fixture {
        let is_target: Vec<bool> = (0..n_rows).map(|i| i % 2 == 0).collect();
        let mut matrix = Vec::with_capacity(n_rows * columns.len());
        for (i, &is_t) in is_target.iter().enumerate() {
            for (_, f) in columns {
                matrix.push(f(i, is_t));
            }
        }
        fixture_from_matrix(
            columns.iter().map(|(n, _)| *n),
            matrix,
            is_target,
            columns.len(),
        )
    }

    /// One named column with explicit labels, for the cases that need a
    /// specific target/decoy split rather than the alternating default.
    pub(crate) fn fixture_one_column(name: &str, values: &[f64], is_target: &[bool]) -> Fixture {
        assert_eq!(values.len(), is_target.len());
        fixture_from_matrix(
            std::iter::once(name),
            values.to_vec(),
            is_target.to_vec(),
            1,
        )
    }

    /// The same rows, from columns given literally rather than generated.
    pub(crate) fn fixture_from_columns(names: &[&str], columns: &[&[f64]]) -> Fixture {
        let n_rows = columns[0].len();
        assert!(columns.iter().all(|c| c.len() == n_rows));
        let mut matrix = Vec::with_capacity(n_rows * columns.len());
        for i in 0..n_rows {
            for c in columns {
                matrix.push(c[i]);
            }
        }
        let is_target = (0..n_rows).map(|i| i % 2 == 0).collect();
        fixture_from_matrix(names.iter().copied(), matrix, is_target, columns.len())
    }

    fn fixture_from_matrix<'a>(
        names: impl Iterator<Item = &'a str>,
        matrix: Vec<f64>,
        is_target: Vec<bool>,
        n_features: usize,
    ) -> Fixture {
        let qvalue: Vec<f32> = is_target
            .iter()
            .map(|&t| if t { 0.001 } else { 0.9 })
            .collect();
        let thresholds = [0.01f32, 0.05, 0.1, 0.5, 1.0]
            .into_iter()
            .map(|q| {
                let (targets, decoys) = qvalue
                    .iter()
                    .zip(&is_target)
                    .filter(|(v, _)| **v <= q)
                    .fold(
                        (0, 0),
                        |(t, d), (_, &is_t)| if is_t { (t + 1, d) } else { (t, d + 1) },
                    );
                curves::ThresholdRow {
                    q,
                    total: targets + decoys,
                    targets,
                    decoys,
                }
            })
            .collect();
        Fixture {
            names: names.map(Arc::from).collect(),
            matrix,
            score: is_target
                .iter()
                .enumerate()
                .map(|(i, &t)| if t { i as f32 } else { -(i as f32) })
                .collect(),
            is_target,
            qvalue,
            thresholds,
            // Distinct per column, so a test sorting on gain can tell which
            // feature it landed on.
            gain: (0..n_features).map(|j| j as f32).collect(),
        }
    }

    fn total(h: &HistView<'_>) -> u32 {
        h.n_target + h.n_decoy
    }

    /// Largest gap between two histograms' CDFs, bin by bin.
    ///
    /// KS and not bin-by-bin total variation, which is the wrong instrument
    /// here: a column taking few distinct values plots as a handful of sharp
    /// spikes, and a hair's difference in the axis bounds slides whole spikes
    /// into neighbouring bins. TV scores that as though all the mass moved,
    /// while the two plots are indistinguishable on screen. KS compares
    /// cumulative mass, so a sub-bin shift costs a sub-bin amount.
    fn ks_distance(a: &[u32], b: &[u32], na: f64, nb: f64) -> f64 {
        let (mut ca, mut cb, mut worst) = (0.0f64, 0.0f64, 0.0f64);
        for (&x, &y) in a.iter().zip(b) {
            ca += x as f64 / na;
            cb += y as f64 / nb;
            worst = worst.max((ca - cb).abs());
        }
        worst
    }

    /// What [`DEFAULT_SAMPLE`]'s doc comment promises, checked rather than
    /// asserted: at the shipped sample size, every stored histogram is within
    /// a KS distance of 0.02 of the whole-data one and every feature AUC
    /// within 0.01.
    ///
    /// `#[ignore]`d on cost alone: two million-row builds, 3.7 s in debug
    /// against 0.2 s for the rest of the suite (0.2 s in release). Run with
    /// `cargo test -p rescore_dash --release -- --ignored`.
    #[test]
    #[ignore = "two million-row builds; run explicitly with --ignored"]
    fn sample_size_holds_its_accuracy_claim() {
        const N_ROWS: usize = 1_000_000;
        let mut rng = Lcg(0xDEAD_BEEF);
        // Deliberately mixed shapes: a clean separator, a heavy tail, a
        // mostly-zero column, and a discrete one. The discrete column is the
        // case that made total variation useless as a metric.
        type RandomColumn<'a> = &'a dyn Fn(&mut Lcg, bool) -> f64;
        let columns: [RandomColumn<'_>; 4] = [
            &|r, t| r.below(1_000_000) as f64 / 1e6 + if t { 0.5 } else { 0.0 },
            &|r, _| (r.below(1_000_000) as f64 / 1e6).powi(-3).min(1e9),
            &|r, _| {
                if r.below(10) == 0 {
                    r.below(1000) as f64
                } else {
                    0.0
                }
            },
            &|r, _| r.below(100) as f64,
        ];

        let mut matrix = Vec::with_capacity(N_ROWS * columns.len());
        let mut is_target = Vec::with_capacity(N_ROWS);
        let mut score = Vec::with_capacity(N_ROWS);
        let mut qvalue = Vec::with_capacity(N_ROWS);
        for i in 0..N_ROWS {
            let t = i % 3 != 0;
            for f in columns {
                matrix.push(f(&mut rng, t));
            }
            is_target.push(t);
            score.push((rng.below(1_000_000) as f64 / 1e6 + if t { 0.6 } else { 0.0 }) as f32);
            qvalue.push(rng.below(1_000_000) as f32 / 1e6);
        }
        let names: Vec<Arc<str>> = (0..columns.len())
            .map(|i| Arc::from(format!("f{i}").as_str()))
            .collect();
        let gain = vec![1.0f32; columns.len()];
        let view = RescoreView {
            feature_names: &names,
            features: &matrix,
            is_target: &is_target,
            score: &score,
            qvalue: &qvalue,
            thresholds: &[],
            gain: &gain,
        };

        let exact = Dashboard::build(&view, usize::MAX).expect("well-formed");
        let sampled = Dashboard::build(&view, DEFAULT_SAMPLE).expect("well-formed");

        let (mut worst_ks, mut worst_auc) = (0.0f64, 0.0f64);
        let (mut ks_at, mut auc_at) = (String::new(), String::new());
        for j in 0..columns.len() {
            let (a, b) = (sampled.auc(j), exact.auc(j));
            if a.is_finite() && b.is_finite() && (a - b).abs() > worst_auc {
                worst_auc = (a - b).abs();
                auc_at = format!("f{j}");
            }
            for t in Axis::ALL {
                let (s, e) = (sampled.hist(j, t, true), exact.hist(j, t, true));
                let (ns, ne) = (total(&s) as f64, total(&e) as f64);
                if ns == 0.0 || ne == 0.0 {
                    continue;
                }
                // Both classes pooled: the panel draws them on one axis, so
                // that is the distribution a reader compares.
                let pool = |h: &HistView<'_>| -> Vec<u32> {
                    h.target.iter().zip(h.decoy).map(|(a, b)| a + b).collect()
                };
                let ks = ks_distance(&pool(&s), &pool(&e), ns, ne);
                if ks > worst_ks {
                    worst_ks = ks;
                    ks_at = format!("f{j} {t:?}");
                }
            }
        }

        assert!(
            worst_ks < 0.02,
            "worst KS {worst_ks:.4} at {ks_at} exceeds the sample-size claim"
        );
        assert!(
            worst_auc < 0.01,
            "worst AUC error {worst_auc:.4} at {auc_at} exceeds the sample-size claim"
        );
    }

    /// Zooming has to re-grid, not re-slice. Every zoom keeps a full
    /// `Q_CURVE_POINTS` grid over its own range, and the y axis rescales with
    /// it — a `q <= 0.01` panel still scaled to the q = 1 target count would
    /// draw the part being zoomed into flat against the axis.
    #[test]
    fn zooming_the_fdr_curve_regrids_and_rescales() {
        let mut f = fixture(400, &[("a", &|i, _| i as f64)]);
        // Targets spread over q in (0, 0.2], so every zoom level has a
        // different number of them below its bound.
        for (i, q) in f.qvalue.iter_mut().enumerate() {
            *q = (i as f32 + 1.0) / 400.0 * 0.2;
        }
        let dash = f.build();
        assert_eq!(dash.n_q_zooms(), Q_ZOOMS.len());

        let mut prev_ymax = f64::INFINITY;
        for (i, &zoom) in Q_ZOOMS.iter().enumerate() {
            let c = dash.q_curve(i);
            assert_eq!(c.zoom, zoom);
            assert_eq!(
                c.points.len(),
                Q_CURVE_POINTS,
                "zoom {zoom} must keep a full grid, not a slice of the wide one"
            );
            assert!(
                (c.points.last().unwrap().0 - zoom).abs() < 1e-12,
                "zoom {zoom} must reach its own bound"
            );
            assert!(
                c.ymax <= prev_ymax,
                "zoom {zoom} y bound {} exceeds the wider view's {prev_ymax}",
                c.ymax
            );
            prev_ymax = c.ymax;
        }
        assert!(
            dash.q_curve(Q_ZOOMS.len() - 1).ymax < dash.q_curve(0).ymax,
            "the tightest zoom must actually rescale y, not just x"
        );

        // The point of re-gridding, stated as a number: slicing the widest
        // curve down to the tightest range leaves almost nothing to plot.
        let tightest = *Q_ZOOMS.last().unwrap();
        let sliced = dash
            .q_curve(0)
            .points
            .iter()
            .filter(|&&(x, _)| x <= tightest)
            .count();
        assert!(
            sliced * 10 < dash.q_curve(Q_ZOOMS.len() - 1).points.len(),
            "a slice would give {sliced} points where the zoom gives {}",
            dash.q_curve(Q_ZOOMS.len() - 1).points.len()
        );
    }

    #[test]
    fn pass_a_matches_a_naive_per_column_sweep() {
        let f = fixture(
            500,
            &[
                ("a", &|i, t| if t { i as f64 } else { -(i as f64) }),
                ("b", &|i, _| (i as f64 * 0.31).sin()),
                ("c", &|i, _| if i % 7 == 0 { f64::NAN } else { i as f64 }),
            ],
        );
        let view = f.view();
        let got = pass_a(&view);
        for (j, g) in got.iter().enumerate().take(3) {
            let mut want = ColumnStats::IDENTITY;
            for i in 0..view.n_rows() {
                want.push(view.row(i)[j], view.is_target[i]);
            }
            assert_eq!(g.n_target, want.n_target, "column {j}");
            assert_eq!(g.n_nan, want.n_nan, "column {j}");
            assert!(
                (g.target_mean - want.target_mean).abs() < 1e-9,
                "column {j}"
            );
            assert!(
                (g.target_var() - want.target_var()).abs() < 1e-9,
                "column {j}"
            );
            assert_eq!((g.lo, g.hi), (want.lo, want.hi), "column {j}");
        }
    }

    /// The score rides along as one more column of the same sweep, so it must
    /// agree with a sweep of the score vector on its own.
    #[test]
    fn pass_a_folds_the_score_in_as_the_last_column() {
        let f = fixture(200, &[("a", &|i, _| i as f64)]);
        let view = f.view();
        let got = pass_a(&view);
        let mut want = ColumnStats::IDENTITY;
        for i in 0..view.n_rows() {
            want.push(view.score[i] as f64, view.is_target[i]);
        }
        let s = &got[view.n_features()];
        assert_eq!(s.n_target, want.n_target);
        assert!((s.target_mean - want.target_mean).abs() < 1e-9);
        assert_eq!((s.lo, s.hi), (want.lo, want.hi));
    }

    /// A sample at least as large as the run takes every row in order, which is
    /// what makes small runs — and these tests — exact rather than approximate.
    #[test]
    fn a_sample_larger_than_the_run_takes_every_row_once() {
        let f = fixture(40, &[("a", &|i, _| i as f64)]);
        let view = f.view();
        let (matrix, labels) = gather_sample(&view, 1000);
        assert_eq!(labels.len(), 40);
        assert_eq!(labels, view.is_target);
        let col: Vec<f64> = matrix.iter().step_by(2).copied().collect();
        assert_eq!(col, (0..40).map(|i| i as f64).collect::<Vec<_>>());
    }

    /// The pipeline hands the dashboard rows sorted descending by score, so a
    /// contiguous sample would be all high-scoring targets. Random indices must
    /// reach across the whole run.
    #[test]
    fn a_smaller_sample_draws_from_the_whole_run_not_a_prefix() {
        let f = fixture(10_000, &[("a", &|i, _| i as f64)]);
        let view = f.view();
        let (matrix, labels) = gather_sample(&view, 500);
        assert_eq!(labels.len(), 500);
        let col: Vec<f64> = matrix.iter().step_by(2).copied().collect();
        let max = col.iter().cloned().fold(f64::MIN, f64::max);
        assert!(
            max > 9_000.0,
            "sample never reached the tail of the run, max was {max}"
        );
        assert!(
            col.iter().any(|&v| v < 1_000.0),
            "sample never reached the head of the run"
        );
    }

    #[test]
    fn the_sample_is_drawn_by_random_index_not_off_the_front() {
        let f = fixture(5_000, &[("a", &|i, _| i as f64)]);
        let view = f.view();
        let (matrix, labels) = gather_sample(&view, 300);
        assert_eq!(labels.len(), 300);

        // The column is `i`, so a sampled value IS its row index. The rows
        // arrive sorted by score, so a contiguous or strided draw would be all
        // high-scoring targets — the reason this is by random index at all.
        let rows: Vec<f64> = matrix
            .iter()
            .step_by(view.n_features() + 1)
            .copied()
            .collect();
        let max = rows.iter().copied().fold(0.0f64, f64::max);
        assert!(
            max > 4_000.0,
            "sample never reached the tail of the run; largest row index {max}"
        );
        let ascending = rows.windows(2).all(|w| w[0] <= w[1]);
        assert!(
            !ascending,
            "sample is a prefix or a stride, not a random draw"
        );
    }

    /// The bug this design exists to kill. One decoy row at 1e6 in a perfectly
    /// separated column reported an AUC of 0.5 under value binning, because the
    /// histogram it was read off binned uniformly over min..max.
    #[test]
    fn a_single_outlier_does_not_collapse_the_auc() {
        let mut values: Vec<f64> = Vec::new();
        let mut is_target = Vec::new();
        for i in 0..100 {
            values.push(i as f64);
            is_target.push(false);
        }
        for i in 100..200 {
            values.push(i as f64);
            is_target.push(true);
        }
        values.push(1e6);
        is_target.push(false);

        let dash = fixture_one_column("outlier_feature", &values, &is_target).build();

        assert!(
            dash.auc(0) > 0.98,
            "one outlier must not collapse the AUC, got {}",
            dash.auc(0)
        );
        // And the fixture really does still defeat the old value-binned path,
        // so this keeps covering the bug rather than becoming a tautology.
        assert!(
            binned_auc(&values, &is_target) < 0.6,
            "fixture no longer reproduces the binning failure"
        );
    }

    /// The other half of the same fix: the clipped axis must exclude the
    /// outlier, while the unclipped axis must still stretch to it. Two views of
    /// the same column, and both are wanted.
    #[test]
    fn clipping_excludes_an_outlier_that_the_unclipped_axis_still_shows() {
        let mut values: Vec<f64> = (0..1000).map(|i| i as f64 / 1000.0).collect();
        values.push(1e6);
        let is_target: Vec<bool> = (0..values.len()).map(|i| i % 2 == 0).collect();
        let dash = fixture_one_column("spiky", &values, &is_target).build();

        let clipped = dash.hist(0, Axis::Value(XTransform::Linear), true);
        assert!(
            clipped.hi < 10.0,
            "clipped axis must trim the 1e6 row, got hi = {}",
            clipped.hi
        );
        let unclipped = dash.hist(0, Axis::Value(XTransform::Linear), false);
        assert!(
            unclipped.hi >= 1e6,
            "unclipped axis is supposed to show the spike, got hi = {}",
            unclipped.hi
        );
        // The bulk survives clipping; the unclipped view crushes it into bin 0.
        assert!(clipped.target.iter().filter(|&&c| c > 0).count() > 100);
        assert_eq!(unclipped.target[0] + unclipped.decoy[0], 1000);
    }

    /// The unclipped lower bound is the smallest *all-rows* value the transform
    /// accepts, not the column minimum. For `log10` that is `min_pos`, which is
    /// exactly why pass A tracks it.
    #[test]
    fn the_unclipped_log10_axis_starts_at_the_smallest_positive_value() {
        let values = vec![-5.0, 0.0, 1e-6, 1.0, 100.0, 42.0];
        let dash = fixture_one_column(
            "mixed_sign",
            &values,
            &[true, false, true, false, true, false],
        )
        .build();
        let h = dash.hist(0, Axis::Value(XTransform::Log10), false);
        assert!((h.lo - (-6.0)).abs() < 1e-9, "log10(1e-6), got {}", h.lo);
        assert!((h.hi - 2.0).abs() < 1e-9, "log10(100), got {}", h.hi);
        // The negative and the zero are outside log10's domain: three of six.
        assert_eq!(total(&h), 4);
    }

    /// `square` is the one non-monotone transform, so its order statistics come
    /// from `|v|`. A column straddling zero must bottom out at the value
    /// closest to zero, not at the square of an endpoint.
    #[test]
    fn the_square_axis_bottoms_out_at_the_value_nearest_zero() {
        let values = vec![-8.0, -0.5, 0.25, 3.0];
        let dash =
            fixture_one_column("straddles_zero", &values, &[true, false, true, false]).build();
        let h = dash.hist(0, Axis::Value(XTransform::Square), false);
        assert!((h.lo - 0.0625).abs() < 1e-9, "0.25^2, got {}", h.lo);
        assert!((h.hi - 64.0).abs() < 1e-9, "(-8)^2, got {}", h.hi);
    }

    #[test]
    fn abs_order_stats_enumerates_magnitudes_in_ascending_order() {
        let sorted: Vec<(f64, bool)> = [-8.0, -3.0, -0.5, 0.25, 2.0, 9.0]
            .into_iter()
            .map(|v| (v, true))
            .collect();
        // |v| ascending is 0.25, 0.5, 2, 3, 8, 9.
        let want = [0.25, 0.5, 2.0, 3.0, 8.0, 9.0];
        for (r, &w) in want.iter().enumerate() {
            let (lo, hi) = abs_order_stats(&sorted, r, r);
            assert_eq!(lo, w, "rank {r}");
            assert_eq!(hi, w, "rank {r}");
        }
        assert_eq!(abs_order_stats(&sorted, 1, 4), (0.5, 8.0));
    }

    #[test]
    fn abs_order_stats_handles_an_all_negative_slice() {
        let sorted: Vec<(f64, bool)> = [-9.0, -4.0, -1.0].into_iter().map(|v| (v, true)).collect();
        assert_eq!(abs_order_stats(&sorted, 0, 2), (1.0, 9.0));
    }

    /// Every stored histogram must account for the sample exactly: binned plus
    /// dropped plus out-of-range is the whole sample, with no value counted
    /// twice and none lost.
    #[test]
    fn every_stored_histogram_accounts_for_the_whole_sample() {
        let f = fixture(
            300,
            &[
                ("mixed", &|i, t| if t { i as f64 } else { -(i as f64) }),
                ("heavy_tail", &|i, _| ((i % 50) as f64).exp()),
                ("all_nan", &|_, _| f64::NAN),
                ("constant", &|_, _| 3.0),
                ("with_zeros", &|i, _| {
                    if i % 3 == 0 { 0.0 } else { i as f64 }
                }),
            ],
        );
        let dash = f.build();
        let m = f.view().n_rows().min(DEFAULT_SAMPLE);
        for column in 0..=dash.n_features() {
            for t in Axis::ALL {
                for clip in [false, true] {
                    let slot = dash.slot_index(column, t, clip);
                    let s = &dash.slots[slot];
                    let h = dash.hist(column, t, clip);
                    let binned = total(&h);
                    assert_eq!(
                        binned as usize + s.dropped as usize + s.n_out as usize,
                        m,
                        "column {column} {t:?} clip={clip} lost or double-counted rows"
                    );
                    assert_eq!(
                        binned,
                        h.target.iter().sum::<u32>() + h.decoy.iter().sum::<u32>(),
                        "stored totals disagree with the stored bins"
                    );
                }
            }
        }
    }

    /// The unclipped axis comes from pass A's exact all-rows values, and the
    /// sample is a subset of those rows, so nothing sampled can fall outside
    /// it — for any transform, including the two that are not monotone.
    ///
    /// This is the invariant that says the unclipped view is genuinely exact
    /// rather than sample-derived. A single out-of-range value here would mean
    /// a transform's exact lower bound was computed from the wrong statistic
    /// (`min` instead of `min_pos` for `log10`, or an endpoint instead of
    /// `min_abs` for `square`), which is precisely the mistake the extra pass-A
    /// fields exist to prevent.
    #[test]
    fn nothing_sampled_ever_falls_outside_the_unclipped_axis() {
        let f = fixture(
            2_000,
            &[
                ("straddles_zero", &|i, _| (i as f64 - 1000.0) / 7.0),
                ("heavy_tail", &|i, _| ((i % 60) as f64 / 4.0).exp()),
                ("tiny_positives", &|i, _| (i as f64 + 1.0) * 1e-9),
                ("with_zeros_and_negatives", &|i, _| match i % 3 {
                    0 => 0.0,
                    1 => -(i as f64),
                    _ => i as f64,
                }),
                ("sparse_nan", &|i, _| {
                    if i % 11 == 0 { f64::NAN } else { i as f64 }
                }),
                ("huge", &|i, _| i as f64 * 1e12),
            ],
        );
        let dash = Dashboard::build(&f.view(), 500).expect("well formed");
        for column in 0..=dash.n_features() {
            for t in Axis::ALL {
                let slot = dash.slot_index(column, t, false);
                assert_eq!(
                    dash.slots[slot].n_out, 0,
                    "column {column} {t:?}: {} sampled values outside an axis \
                     that is supposed to bound every row",
                    dash.slots[slot].n_out
                );
            }
        }
    }

    /// An all-NaN column has nothing any transform can plot. It must come back
    /// empty and labelled, not with a fabricated range.
    #[test]
    fn a_wholly_non_finite_column_is_marked_unplottable() {
        let f = fixture(
            50,
            &[("real", &|i, _| i as f64), ("all_nan", &|_, _| f64::NAN)],
        );
        let dash = f.build();
        for t in Axis::ALL {
            for clip in [false, true] {
                let h = dash.hist(1, t, clip);
                assert_eq!(total(&h), 0, "{t:?} clip={clip}");
                assert!(
                    dash.subtitle(1, t, clip).contains("nothing this transform"),
                    "{t:?} clip={clip}: {}",
                    dash.subtitle(1, t, clip)
                );
            }
        }
        assert!(dash.auc(1).is_nan());
        // The neighbouring real column is unaffected.
        assert!(total(&dash.hist(0, Axis::Value(XTransform::Linear), true)) > 0);
    }

    /// `RankPercentile` needs no transform and no range search: sorted position
    /// *is* the percentile, so the histogram spans 0..100 and is close to flat.
    #[test]
    fn rank_percentile_spans_the_percentile_axis_and_is_near_uniform() {
        let f = fixture(2_000, &[("heavy_tail", &|i, _| ((i % 40) as f64).exp())]);
        let dash = f.build();
        let h = dash.hist(0, Axis::RankPercentile, false);
        assert_eq!((h.lo, h.hi), (0.0, 100.0));
        let occupied = h
            .target
            .iter()
            .zip(h.decoy)
            .filter(|&(&a, &b)| a + b > 0)
            .count();
        assert!(
            occupied > 20,
            "a heavy tail must still spread across the rank axis, got {occupied} bins"
        );
    }

    /// A constant column has no range to plot. Both axes get unit width rather
    /// than collapsing onto a single edge, and every row lands in bin 0.
    #[test]
    fn a_constant_column_gets_a_unit_wide_axis() {
        let f = fixture(100, &[("constant", &|_, _| 3.0)]);
        let dash = f.build();
        for clip in [false, true] {
            let h = dash.hist(0, Axis::Value(XTransform::Linear), clip);
            assert_eq!((h.lo, h.hi), (3.0, 4.0), "clip={clip}");
            assert_eq!(h.target[0] + h.decoy[0], 100, "clip={clip}");
        }
    }

    #[test]
    fn the_score_auc_is_exact_and_the_score_gets_its_own_histograms() {
        // The fixture's score separates the classes perfectly, so the exact
        // AUC is 1 by definition. A sampled AUC could only reach it by luck.
        let f = fixture(400, &[("a", &|i, _| i as f64)]);
        let dash = f.build();
        assert_eq!(
            dash.auc(dash.score_column()),
            1.0,
            "the score AUC is not sampled"
        );
        let tiny = Dashboard::build(&f.view(), 20).expect("well-formed");
        assert_eq!(
            tiny.auc(tiny.score_column()),
            1.0,
            "a small sample must not reach the score's AUC"
        );
        assert!(total(&dash.hist(dash.score_column(), Axis::Value(XTransform::Linear), true)) > 0);
        assert!(
            dash.title(dash.score_column(), Axis::Value(XTransform::Linear))
                .contains("discriminant_score")
        );
    }

    #[test]
    fn table_cells_are_formatted_once_per_feature() {
        let f = fixture(
            80,
            &[("alpha", &|i, _| i as f64), ("beta", &|_, _| f64::NAN)],
        );
        let dash = f.build();
        assert_eq!(dash.cells.len(), 2);
        assert_eq!(dash.cells[0][0], "alpha");
        assert_eq!(dash.cells[1][3], "-", "a NaN AUC prints as a dash");
        assert_eq!(dash.cells[1][5], "1.0000", "the NaN fraction is 1");
    }

    /// A reader has to be able to tell an all-rows histogram from a sampled
    /// one. It is a fact about the run, so it lives on the dashboard rather
    /// than being repeated in every panel's subtitle.
    #[test]
    fn the_basis_says_whether_the_histograms_were_sampled() {
        let f = fixture(80, &[("a", &|i, _| i as f64)]);
        let dash = f.build();
        assert_eq!(dash.basis(), "all 80 rows");

        let big = fixture(2_000, &[("a", &|i, _| i as f64)]);
        let dash = Dashboard::build(&big.view(), 500).expect("well formed");
        assert_eq!(dash.basis(), "2000 rows, histograms from a 500 sample");
    }

    /// The subtitle has to fit the histogram panel, which is under half of a
    /// terminal's width. Anything longer is silently cut off by the border it
    /// is drawn on, which is how the run-level facts ended up on the banner.
    #[test]
    fn subtitles_fit_a_narrow_panel() {
        let f = fixture(200, &[("a", &|i, _| i as f64 - 100.0)]);
        let dash = f.build();
        for t in Axis::ALL {
            for clip in [false, true] {
                let s = dash.subtitle(0, t, clip);
                assert!(s.len() <= 60, "{t:?} clip={clip}: {} chars: {s}", s.len());
            }
        }
    }

    /// `log10` rejects the non-positive half, and the subtitle is the only
    /// place that says so — a histogram that silently shrank would read as
    /// "those rows do not exist".
    #[test]
    fn the_subtitle_reports_what_a_transform_refused() {
        let f = fixture(100, &[("half_negative", &|i, _| i as f64 - 50.0)]);
        let dash = f.build();
        let sub = dash.subtitle(0, Axis::Value(XTransform::Log10), true);
        assert!(sub.contains("dropped 51.0% by transform"), "{sub}");
    }

    #[test]
    fn build_rejects_a_malformed_view() {
        let names: Vec<Arc<str>> = vec![Arc::from("a")];
        let view = RescoreView {
            feature_names: &names,
            features: &[1.0, 2.0],
            is_target: &[true, false],
            score: &[0.0],
            qvalue: &[0.0, 0.0],
            thresholds: &[],
            gain: &[0.0],
        };
        assert!(matches!(
            Dashboard::build(&view, DEFAULT_SAMPLE),
            Err(ViewError::RowLen { .. })
        ));
    }

    #[test]
    fn fmt_count_scales_its_unit() {
        assert_eq!(fmt_count(42), "42");
        assert_eq!(fmt_count(9_999), "9999");
        assert_eq!(fmt_count(250_000), "250k");
        assert_eq!(fmt_count(2_013_441), "2.01M");
    }
}
