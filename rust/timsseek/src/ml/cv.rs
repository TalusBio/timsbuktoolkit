pub use super::TargetDecoy;
use crate::scoring::timings::TimedStep;
pub use forust_ml::constraints::{
    Constraint,
    ConstraintMap,
};
pub use forust_ml::errors::ForustError;
pub use forust_ml::gradientbooster::{
    GrowPolicy,
    ImportanceMethod,
    MissingNodeTreatment,
};
pub use forust_ml::metric::{
    EvaluationMetric,
    Metric,
};
pub use forust_ml::objective::ObjectiveType;
pub use forust_ml::sampler::SampleMethod;
pub use forust_ml::{
    GradientBooster,
    Matrix,
};
use serde::Serialize;
pub use std::collections::{
    HashMap,
    HashSet,
};
use std::sync::Arc;

pub struct GBMConfig {
    iterations: usize,
    learning_rate: f32,
    max_depth: usize,
    max_leaves: usize,
    l1: f32,
    l2: f32,
    gamma: f32,
    max_delta_step: f32,
    min_leaf_weight: f32,
    base_score: f64,
    nbins: u16,
    parallel: bool,
    allow_missing_splits: bool,
    monotone_constraints: Option<ConstraintMap>,
    subsample: f32,
    top_rate: f64,
    other_rate: f64,
    colsample_bytree: f64,
    seed: u64,
    missing: f64,
    create_missing_branch: bool,
    sample_method: SampleMethod,
    grow_policy: GrowPolicy,
    evaluation_metric: Option<Metric>,
    early_stopping_rounds: Option<usize>,
    initialize_base_score: bool,
    terminate_missing_features: HashSet<usize>,
    missing_node_treatment: MissingNodeTreatment,
    log_iterations: usize,
    force_children_to_bound_parent: bool,
}

impl Clone for GBMConfig {
    fn clone(&self) -> Self {
        let Self {
            iterations,
            learning_rate,
            max_depth,
            max_leaves,
            l1,
            l2,
            gamma,
            max_delta_step,
            min_leaf_weight,
            base_score,
            nbins,
            parallel,
            allow_missing_splits,
            monotone_constraints,
            subsample,
            top_rate,
            other_rate,
            colsample_bytree,
            seed,
            missing,
            create_missing_branch,
            sample_method,
            grow_policy,
            evaluation_metric,
            early_stopping_rounds,
            initialize_base_score,
            terminate_missing_features,
            missing_node_treatment,
            log_iterations,
            force_children_to_bound_parent,
        } = self;

        Self {
            iterations: *iterations,
            learning_rate: *learning_rate,
            max_depth: *max_depth,
            max_leaves: *max_leaves,
            l1: *l1,
            l2: *l2,
            gamma: *gamma,
            max_delta_step: *max_delta_step,
            min_leaf_weight: *min_leaf_weight,
            base_score: *base_score,
            nbins: *nbins,
            parallel: *parallel,
            allow_missing_splits: *allow_missing_splits,
            monotone_constraints: monotone_constraints.clone(),
            subsample: *subsample,
            top_rate: *top_rate,
            other_rate: *other_rate,
            colsample_bytree: *colsample_bytree,
            seed: *seed,
            missing: *missing,
            create_missing_branch: *create_missing_branch,
            sample_method: match sample_method {
                SampleMethod::None => SampleMethod::None,
                SampleMethod::Random => SampleMethod::Random,
                SampleMethod::Goss => SampleMethod::Goss,
            },
            grow_policy: match grow_policy {
                GrowPolicy::DepthWise => GrowPolicy::DepthWise,
                GrowPolicy::LossGuide => GrowPolicy::LossGuide,
            },
            evaluation_metric: *evaluation_metric,
            early_stopping_rounds: *early_stopping_rounds,
            initialize_base_score: *initialize_base_score,
            terminate_missing_features: terminate_missing_features.clone(),
            missing_node_treatment: *missing_node_treatment,
            log_iterations: *log_iterations,
            force_children_to_bound_parent: *force_children_to_bound_parent,
        }
    }
}

impl Default for GBMConfig {
    fn default() -> Self {
        GBMConfig {
            iterations: 1000,
            learning_rate: 0.1,
            max_depth: 6,
            max_leaves: usize::MAX,
            l1: 0.,
            l2: 1.,
            gamma: 0.,
            max_delta_step: 0.,
            min_leaf_weight: 5.,
            base_score: 0.5,
            nbins: 256,
            parallel: true,
            allow_missing_splits: true,
            monotone_constraints: None,
            subsample: 0.8,
            top_rate: 0.1,
            other_rate: 0.2,
            colsample_bytree: 0.8,
            seed: 0,
            missing: f64::NAN,
            create_missing_branch: false,
            sample_method: SampleMethod::Random,
            grow_policy: GrowPolicy::DepthWise,
            evaluation_metric: Some(Metric::LogLoss),
            // evaluation_metric: None,
            early_stopping_rounds: Some(100),
            initialize_base_score: true,
            terminate_missing_features: HashSet::new(),
            missing_node_treatment: MissingNodeTreatment::AssignToParent,
            log_iterations: 50,
            force_children_to_bound_parent: false,
        }
    }
}

impl GBMConfig {
    fn try_build(&self) -> Result<GradientBooster, ForustError> {
        let Self {
            iterations,
            learning_rate,
            max_depth,
            max_leaves,
            l1,
            l2,
            gamma,
            max_delta_step,
            min_leaf_weight,
            base_score,
            nbins,
            parallel,
            allow_missing_splits,
            monotone_constraints,
            subsample,
            top_rate,
            other_rate,
            colsample_bytree,
            seed,
            missing,
            create_missing_branch,
            sample_method,
            grow_policy,
            evaluation_metric,
            early_stopping_rounds,
            initialize_base_score,
            terminate_missing_features,
            missing_node_treatment,
            log_iterations,
            force_children_to_bound_parent,
        } = self;

        GradientBooster::new(
            ObjectiveType::LogLoss,
            *iterations,
            *learning_rate,
            *max_depth,
            *max_leaves,
            *l1,
            *l2,
            *gamma,
            *max_delta_step,
            *min_leaf_weight,
            *base_score,
            *nbins,
            *parallel,
            *allow_missing_splits,
            monotone_constraints.clone(),
            *subsample,
            *top_rate,
            *other_rate,
            *colsample_bytree,
            *seed,
            *missing,
            *create_missing_branch,
            match sample_method {
                SampleMethod::None => SampleMethod::None,
                SampleMethod::Random => SampleMethod::Random,
                SampleMethod::Goss => SampleMethod::Goss,
            },
            match grow_policy {
                GrowPolicy::DepthWise => GrowPolicy::DepthWise,
                GrowPolicy::LossGuide => GrowPolicy::LossGuide,
            },
            *evaluation_metric,
            *early_stopping_rounds,
            *initialize_base_score,
            terminate_missing_features.clone(),
            *missing_node_treatment,
            *log_iterations,
            *force_children_to_bound_parent,
        )
    }
}

/// Label + score access for a scored record — everything
/// [`CrossValidatedScorer`] needs from `T` once the feature matrix exists.
pub trait FeatureLike {
    fn get_y(&self) -> f64;
    fn assign_score(&mut self, score: f64) -> ();
    fn get_score(&self) -> f64;
}

/// A record that can also walk *itself* into a feature vector. Required only by
/// the self-computing ctor ([`CrossValidatedScorer::new_from_shuffled`]); the
/// production path hands in an externally built matrix
/// (`PrecomputedFeatures::from_row_major`) and needs [`FeatureLike`] alone.
pub trait FeatureVector: FeatureLike {
    /// Note: the returned iterator MUST yield exactly N elements
    /// for every element of this type.
    fn as_feature(&self) -> impl IntoIterator<Item = f64> + '_;
}

/// Read-only, random-access view of the rows being cross-fitted.
///
/// Everything a [`FoldModel`] is allowed to know about the data: the width and
/// names of the lane matrix, THE fold partition, one row's values, and one
/// row's target/decoy label. Deliberately no `ncols` (it is
/// `column_names().len()`) and no `names` argument threaded through the model
/// impls — the dataset already knows both, so there is no second copy to drift.
pub trait FoldDataset {
    fn nrows(&self) -> usize;
    /// Lane column names; also defines the row width.
    fn column_names(&self) -> Vec<Arc<str>>;
    /// THE partition. One definition, one place.
    fn get_fold(&self, i: usize) -> usize;
    /// Writes row `i` into `out` (`out.len() == column_names().len()`).
    fn get_values(&self, i: usize, out: &mut [f64]);
    fn is_decoy(&self, i: usize) -> bool;
}

/// A model [`CrossValidatedScorer`] can cross-fit: fit on a row-index slice,
/// score another, and report per-column importance.
///
/// `fit` receives BOTH the training rows and an early-stopping (`val`) slice.
/// A model without early stopping is expected to ignore `val`; only the
/// scorer's guarantee matters, namely that neither slice is ever scored by the
/// model fitted from it.
pub trait FoldModel: Sized {
    type Config;
    type Error;
    /// Fit one fold's model.
    ///
    /// `fold` is the identity of the fold being fitted, and it is a PARAMETER
    /// rather than config state because the caller is the only thing that knows
    /// it: a config is shared across every fold, and deriving the fold from
    /// `train[0]` breaks on an empty slice and is meaningless for a partition
    /// that trains on several folds at once (`qvalues::crossfit_lda` trains on
    /// all folds but one). Models with a stochastic initialization mix it into
    /// their seed, so that folds differ from each other while a rerun of the
    /// same fold does not. Deterministic models ignore it, the same way a model
    /// without early stopping ignores `val`.
    fn fit<D: FoldDataset>(
        cfg: &Self::Config,
        data: &D,
        fold: usize,
        train: &[usize],
        val: &[usize],
    ) -> Result<Self, Self::Error>;
    fn predict<D: FoldDataset>(&self, data: &D, rows: &[usize]) -> Result<Vec<f64>, Self::Error>;
    /// Lane-indexed, always `column_names().len()` long, so the vector is
    /// positionally comparable across folds and across models by construction.
    ///
    /// # The NAN-vs-finite contract
    /// * `f32::NAN` means **"this model does not report a value for this
    ///   column"** — the column is absent from the sidecar for this fold.
    /// * Any FINITE value, **including `0.0`**, is a reported measurement and
    ///   reaches the sidecar as-is.
    ///
    /// The distinction is load-bearing and the two cases are NOT the same
    /// thing. A tree model that never split on a column has measured nothing
    /// about it (`NAN`); an LDA whose `|coef|` came out exactly `0.0` for a
    /// dead or constant column has measured precisely that, and dropping it
    /// would delete the columns an operator most wants to see. `0.0` used to
    /// carry both meanings, and the `!= 0.0` filter at the [`fold_feature_stats`]
    /// boundary silently applied the first model's semantics to every model.
    fn importance(&self) -> Vec<f32>;
}

/// The one [`FoldDataset`] impl: the already-materialized row-major matrix that
/// [`PrecomputedFeatures`] holds, plus its column names and the fold count that
/// defines the partition. `get_values` is a `copy_from_slice` out of a slab
/// that is already in memory.
pub struct RowMajorDataset {
    features: PrecomputedFeatures,
    names: Vec<Arc<str>>,
    n_folds: usize,
}

impl RowMajorDataset {
    pub(crate) fn new(features: PrecomputedFeatures, names: Vec<Arc<str>>, n_folds: usize) -> Self {
        assert!(n_folds > 0, "n_folds must be positive");
        assert_eq!(
            names.len(),
            features.ncols,
            "column names must be one per matrix column"
        );
        // `is_decoy` is the ONLY label channel the model traits get, so the
        // labels have to be binary for it to be lossless. Both producers
        // (`FeatureVector::get_y` walks and `from_row_major`) feed 0.0/1.0.
        debug_assert!(
            features.responses.iter().all(|&y| y == 0.0 || y == 1.0),
            "FoldDataset labels must be binary 0.0/1.0"
        );
        Self {
            features,
            names,
            n_folds,
        }
    }
}

impl FoldDataset for RowMajorDataset {
    fn nrows(&self) -> usize {
        self.features.responses.len()
    }

    fn column_names(&self) -> Vec<Arc<str>> {
        self.names.clone()
    }

    fn get_fold(&self, i: usize) -> usize {
        i % self.n_folds
    }

    fn get_values(&self, i: usize, out: &mut [f64]) {
        out.copy_from_slice(self.features.row(i));
    }

    fn is_decoy(&self, i: usize) -> bool {
        self.features.responses[i] <= 0.5
    }
}

/// Materialized view into a `DataBuffer`: feature matrix + label slice.
type FoldView<'a> = (Matrix<'a, f64>, &'a [f64]);

/// Column-major scratch space for one `forust` call. `forust` wants
/// feature-major, so gathering a set of rows out of the row-major dataset is a
/// transpose; this is where it happens and the only place it does.
#[derive(Default)]
struct DataBuffer {
    fold_buffer: Vec<f64>,
    response_buffer: Vec<f64>,
    nrows: usize,
    ncols: usize,
}

impl DataBuffer {
    /// Gather `rows` (in the given order) out of `data` into feature-major
    /// layout: `fold_buffer[feature_idx * nrows + sample_idx]`.
    fn fill_from<D: FoldDataset>(&mut self, data: &D, rows: &[usize], ncols: usize) {
        self.fold_buffer.clear();
        self.response_buffer.clear();
        self.nrows = rows.len();
        self.ncols = ncols;
        self.fold_buffer.resize(self.ncols * self.nrows, 0.0);

        let mut row = vec![0.0f64; ncols];
        for (sample_idx, &elem_idx) in rows.iter().enumerate() {
            data.get_values(elem_idx, &mut row);
            for (feature_idx, &val) in row.iter().enumerate() {
                self.fold_buffer[feature_idx * self.nrows + sample_idx] = val;
            }
            self.response_buffer
                .push(if data.is_decoy(elem_idx) { 0.0 } else { 1.0 });
        }
    }

    fn as_matrix(&self) -> FoldView<'_> {
        let mat = Matrix::new(self.fold_buffer.as_slice(), self.nrows, self.ncols);
        assert_eq!(self.fold_buffer.len(), self.nrows * self.ncols);
        assert_eq!(mat.rows, self.nrows);
        assert_eq!(self.response_buffer.len(), self.nrows);
        (mat, self.response_buffer.as_slice())
    }
}

/// Class weights: decoys 1.0, targets 0.5 — the historical GBM weighting,
/// unchanged, just expressed through the dataset's label channel.
///
/// `pub(crate)` because [`crate::ml::mlp_fold`] weights its training rows the
/// same way and used to re-inline the two literals. Two copies of a training
/// objective is not a duplication a test would catch: both models would simply
/// be fitted against different weightings, and nothing downstream compares them.
pub(crate) fn fold_weights<D: FoldDataset>(data: &D, rows: &[usize]) -> Vec<f64> {
    rows.iter()
        .map(|&i| if data.is_decoy(i) { 1.0 } else { 0.5 })
        .collect()
}

/// [`FoldModel`] adapter for `forust`'s [`GradientBooster`]. A newtype only
/// because `GradientBooster` is a foreign type; it also carries the lane width
/// so [`FoldModel::importance`] can return a full-width, lane-indexed vector
/// (forust reports only the columns it split on).
pub struct GbmFoldModel {
    booster: GradientBooster,
    ncols: usize,
}

impl FoldModel for GbmFoldModel {
    type Config = GBMConfig;
    type Error = ForustError;

    /// `fold` is DELIBERATELY IGNORED: forust seeds itself from
    /// `GBMConfig::seed`, so mixing the fold in here would fight the one seed
    /// the booster already has. Folds still differ, because they are fitted on
    /// different rows.
    fn fit<D: FoldDataset>(
        cfg: &GBMConfig,
        data: &D,
        fold: usize,
        train: &[usize],
        val: &[usize],
    ) -> Result<Self, ForustError> {
        let _ = fold;
        let ncols = data.column_names().len();
        let mut train_buffer = DataBuffer::default();
        let mut val_buffer = DataBuffer::default();
        train_buffer.fill_from(data, train, ncols);
        val_buffer.fill_from(data, val, ncols);

        let (matrix, response) = train_buffer.as_matrix();
        let (v_matrix, v_response) = val_buffer.as_matrix();

        let eval_weight = fold_weights(data, val);
        assert_eq!(eval_weight.len(), v_response.len());
        let evaluation_data = Some(vec![(v_matrix, v_response, eval_weight.as_ref())]);

        let sample_weight = fold_weights(data, train);
        assert_eq!(sample_weight.len(), response.len());

        let mut booster = cfg.try_build()?;
        booster.fit(&matrix, response, &sample_weight, evaluation_data)?;
        Ok(Self { booster, ncols })
    }

    fn predict<D: FoldDataset>(&self, data: &D, rows: &[usize]) -> Result<Vec<f64>, ForustError> {
        let mut buffer = DataBuffer::default();
        buffer.fill_from(data, rows, self.ncols);
        let (matrix, _response) = buffer.as_matrix();
        Ok(self.booster.predict(&matrix, true))
    }

    /// forust reports only the columns it actually split on. Every other lane
    /// is `NAN` — "not reported" — rather than `0.0`, because this model has no
    /// gain measurement for a column it never used. See the
    /// [`FoldModel::importance`] contract.
    fn importance(&self) -> Vec<f32> {
        let raw = self
            .booster
            .calculate_feature_importance(ImportanceMethod::Gain, true);
        debug_assert!(
            raw.keys().all(|&idx| idx < self.ncols),
            "forust importance index exceeds the lane width"
        );
        let mut out = vec![f32::NAN; self.ncols];
        for (idx, gain) in raw {
            if let Some(slot) = out.get_mut(idx) {
                // A NAN here would be indistinguishable from "not reported".
                // forust has no reason to emit one; pin that rather than let it
                // silently delete a row from the sidecar.
                debug_assert!(
                    gain.is_finite(),
                    "forust reported a non-finite gain for column {idx}"
                );
                *slot = gain;
            }
        }
        out
    }
}

/// Row-major precomputed feature matrix: features[sample_idx * ncols + feature_idx]
pub(crate) struct PrecomputedFeatures {
    features: Vec<f64>,
    responses: Vec<f64>,
    ncols: usize,
}

impl PrecomputedFeatures {
    fn from_data(data: &[impl FeatureVector]) -> Self {
        let ncols = data
            .first()
            .map_or(0, |d| d.as_feature().into_iter().count());
        let mut features = Vec::with_capacity(data.len() * ncols);
        let mut responses = Vec::with_capacity(data.len());
        for elem in data {
            features.extend(elem.as_feature());
            responses.push(elem.get_y());
        }
        Self {
            features,
            responses,
            ncols,
        }
    }

    /// Build from an already-materialized row-major matrix
    /// (`features[i*ncols + j]`) + responses, instead of walking
    /// [`FeatureVector::as_feature`]. Rows MUST align with the `data` the
    /// scorer is constructed from (same order). This is how the lane-matrix
    /// consumer trains GBM on a prebuilt lane feature set.
    pub(crate) fn from_row_major(features: Vec<f64>, ncols: usize, responses: Vec<f64>) -> Self {
        assert_eq!(
            features.len(),
            ncols.saturating_mul(responses.len()),
            "row-major feature buffer must be nrows*ncols"
        );
        Self {
            features,
            responses,
            ncols,
        }
    }

    fn row(&self, idx: usize) -> &[f64] {
        &self.features[idx * self.ncols..(idx + 1) * self.ncols]
    }
}

/// Per-feature summary within a fold. `mean` is computed over finite
/// values only (NaN/+-Inf skipped). `nan_ratio` is the fraction of
/// non-finite values seen for this feature in the fold (0.0..=1.0).
#[derive(Debug, Serialize)]
pub struct FeatureStat {
    pub name: Arc<str>,
    pub mean: f32,
    pub nan_ratio: f32,
}

/// Per-fold feature statistics. `feature_stats` preserves `feature_names`
/// insertion order. `feature_importance` sorted by gain descending
/// (top features first) so the JSON reads top-down.
#[derive(Debug, Serialize)]
pub struct FoldStats {
    pub fold: u8,
    pub feature_stats: Vec<FeatureStat>,
    pub feature_importance: Vec<(Arc<str>, f32)>,
}

pub type RescoreFeatureStats = Vec<FoldStats>;

/// Per-fold feature means/NaN ratios + model importance, for ANY fold
/// partition.
///
/// Partition-agnostic on purpose: `fold_rows[f]` is simply "the rows summarized
/// under fold `f`" and `models[f]` is "the model whose importance is reported
/// there". The two rescoring partitions in this crate disagree on both — the
/// [`CrossValidatedScorer`] fits on fold `f` and scores the others, while
/// `qvalues::crossfit_lda` fits on everything BUT fold `f` and scores only fold
/// `f` — and both are leak-free. Keeping this function ignorant of which one it
/// is handed is what lets the sidecar have one implementation without the two
/// partitions being forced to converge.
///
/// Column names (and therefore the row width) come from `data`, so the stats
/// align with the matrix the model saw by construction.
pub(crate) fn fold_feature_stats<D: FoldDataset, M: FoldModel>(
    data: &D,
    fold_rows: &[Vec<usize>],
    models: &[Option<&M>],
) -> RescoreFeatureStats {
    let names = data.column_names();
    let mut out: RescoreFeatureStats = Vec::with_capacity(fold_rows.len());
    let mut row_buf = vec![0.0f64; names.len()];
    for (fold, rows) in fold_rows.iter().enumerate() {
        // --- Importance, back to the (name, gain) sidecar shape ---
        // UNREPORTED (`NAN`) columns are dropped; every REPORTED column is
        // emitted, `0.0` included. See the [`FoldModel::importance`] contract
        // for why those are different things. The drop matters because the
        // sidecar and the dashboard's fold-averaged gain both treat a feature
        // as "reported by this fold" simply by being present, so a model's
        // unmeasured columns would otherwise pad the averaging divisor with
        // values it never produced.
        //
        // The filter is `is_finite`, NOT `!= 0.0`: the old rule was written for
        // the tree model, where 0.0 did mean "never split on this", and then
        // silently deleted an LDA's genuine zero-weight columns — exactly the
        // dead/constant features an operator reads the sidecar to find.
        let importance: Vec<(Arc<str>, f32)> = match models.get(fold).copied().flatten() {
            Some(model) => {
                let raw_imp = model.importance();
                debug_assert_eq!(
                    raw_imp.len(),
                    names.len(),
                    "FoldModel::importance must be lane-indexed"
                );
                let mut pairs: Vec<(Arc<str>, f32)> = raw_imp
                    .into_iter()
                    .zip(names.iter())
                    .filter(|(v, _)| v.is_finite())
                    .map(|(v, n)| (n.clone(), v))
                    .collect();
                // Every surviving value is finite, so this comparison is total;
                // the sort is stable, so equal gains (e.g. a run of genuine
                // 0.0s) keep lane order.
                pairs.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
                pairs
            }
            None => Vec::new(),
        };

        // --- Per-feature finite-only means + NaN ratio ---
        let mut sums: Vec<f64> = vec![0.0; names.len()];
        let mut finite_counts: Vec<u32> = vec![0; names.len()];
        let mut nan_counts: Vec<u32> = vec![0; names.len()];
        let mut n: usize = 0;
        for &i in rows.iter() {
            // Read from the dataset (the exact columns the model trained
            // on) rather than re-walking a per-record feature fn, so the
            // stats align with `names` by construction even when features
            // are supplied externally.
            data.get_values(i, &mut row_buf);
            for (j, &v) in row_buf.iter().enumerate() {
                if v.is_finite() {
                    sums[j] += v;
                    finite_counts[j] += 1;
                } else {
                    nan_counts[j] += 1;
                }
            }
            n += 1;
        }
        let feature_stats: Vec<FeatureStat> = if n > 0 {
            names
                .iter()
                .zip(sums.iter())
                .zip(finite_counts.iter())
                .zip(nan_counts.iter())
                .map(|(((name, s), fc), nc)| {
                    let mean = if *fc > 0 { *s / *fc as f64 } else { f64::NAN };
                    FeatureStat {
                        name: name.clone(),
                        mean: mean as f32,
                        nan_ratio: *nc as f32 / n as f32,
                    }
                })
                .collect()
        } else {
            Vec::new()
        };

        out.push(FoldStats {
            fold: fold as u8,
            feature_stats,
            feature_importance: importance,
        });
    }
    out
}

/// This is meant to IN ESSENCE ...
///
/// Provided we have a number of splits k >= 3.
/// We would train k classifiers, for for classifier n,
/// we would use as train data fold n, as early stop data
/// fold n + 1 and thus we can safely use as test/scoring
/// data any other fold.
///
/// So the score for any point in the data is the average of
/// the results for all classifiers that didint use it
/// for either training or early_stopping_rounds.
///
/// `M` is the model being cross-fitted ([`FoldModel`]); the partition above is
/// the scorer's and does not vary with it.
pub struct CrossValidatedScorer<T: FeatureLike, M: FoldModel> {
    n_folds: u8,
    data: Vec<T>,
    dataset: RowMajorDataset,
    /// `fold_rows[f]` = the row indices with `dataset.get_fold(i) == f`, in
    /// ascending order. Derived from `get_fold` so the partition still has
    /// exactly one definition; materialized once because every fit and every
    /// scoring pass needs it.
    fold_rows: Vec<Vec<usize>>,
    fold_models: Vec<Option<M>>,
    config: M::Config,
}

impl<T: FeatureLike, M: FoldModel> CrossValidatedScorer<T, M> {
    /// Create a new CrossValidatedScorer
    ///
    /// NOTE: THIS ASSUMES YOUR DATA IS ALREADY SHUFFLED
    /// FOLDS WILL BE ASSIGNED IN ORDER (0, 1, 2, ..., n_folds-1, 0, 1, ...)
    /// IF YOUR DATA IS ORDERED IN ANY WAY, COULD LEAD TO BIASED RESULTS.
    ///
    /// The self-walking ctor has no name source, so columns are named
    /// positionally (`feature_0`, ...). Callers that have real names use
    /// [`Self::new_from_shuffled_with_precomputed`].
    pub fn new_from_shuffled(n_folds: u8, data: Vec<T>, config: M::Config) -> Self
    where
        T: FeatureVector,
    {
        let precomputed = PrecomputedFeatures::from_data(&data);
        let names: Vec<Arc<str>> = (0..precomputed.ncols)
            .map(|j| Arc::from(format!("feature_{j}").as_str()))
            .collect();
        Self::new(n_folds, data, config, precomputed, names)
    }

    /// Like [`Self::new_from_shuffled`] but with the feature matrix supplied
    /// externally (already row-major, row-aligned with `data`) along with its
    /// column names. Labels + fold assignment are positional, so the caller
    /// MUST build `precomputed` from the SAME (already-shuffled) `data` order.
    pub(crate) fn new_from_shuffled_with_precomputed(
        n_folds: u8,
        data: Vec<T>,
        config: M::Config,
        precomputed: PrecomputedFeatures,
        names: Vec<Arc<str>>,
    ) -> Self {
        Self::new(n_folds, data, config, precomputed, names)
    }

    fn new(
        n_folds: u8,
        data: Vec<T>,
        config: M::Config,
        precomputed: PrecomputedFeatures,
        names: Vec<Arc<str>>,
    ) -> Self {
        assert_eq!(
            precomputed.responses.len(),
            data.len(),
            "precomputed rows must align 1:1 with data"
        );
        let dataset = RowMajorDataset::new(precomputed, names, n_folds as usize);
        let mut fold_rows: Vec<Vec<usize>> = vec![Vec::new(); n_folds as usize];
        for i in 0..dataset.nrows() {
            fold_rows[dataset.get_fold(i)].push(i);
        }

        Self {
            n_folds,
            data,
            dataset,
            fold_rows,
            fold_models: Vec::new(),
            config,
        }
    }

    pub fn fit(&mut self) -> Result<(), M::Error> {
        self.fit_step()?;
        self.assign_scores()?;

        Ok(())
    }

    pub fn fit_step(&mut self) -> Result<(), M::Error> {
        self.fold_models.clear();
        (0..self.n_folds).for_each(|_| self.fold_models.push(None));
        for fold in 0..self.n_folds {
            let step = TimedStep::begin_stderr(format_args!(
                "  Training fold {}/{}...",
                fold + 1,
                self.n_folds
            ));
            self.fit_fold(fold)?;
            step.finish();
        }

        Ok(())
    }

    pub fn get_scores(&self) -> Result<Vec<f64>, M::Error> {
        // Each row is scored by every fold classifier except its own training
        // fold and that fold's early-stopping fold, i.e. `n_folds - 2` times.
        // At n_folds == 2 that averaging divisor is 0 and the scores silently
        // go inf/NaN instead of erroring, so pin the floor here.
        debug_assert!(
            self.n_folds >= 3,
            "get_scores averages over n_folds - 2 classifiers; n_folds must be >= 3, got {}",
            self.n_folds
        );
        let mut scores = vec![0.0; self.data.len()];

        let step = TimedStep::begin_stderr("  Scoring folds...");
        for train_i in 0..self.n_folds {
            let early_stop_i = self.next_fold(train_i);

            for inference_i in 0..self.n_folds {
                if inference_i == train_i {
                    continue;
                };
                if inference_i == early_stop_i {
                    continue;
                };
                let scorer = self.fold_models[train_i as usize].as_ref().unwrap();
                let rows = &self.fold_rows[inference_i as usize];
                let preds = scorer.predict(&self.dataset, rows)?;
                // Now we need to add the predictions to the scores
                for (&row, pred) in rows.iter().zip(preds) {
                    scores[row] += pred;
                }
            }
        }
        step.finish();

        let div_factor = (self.n_folds - 2) as f64;
        scores.iter_mut().for_each(|x| {
            *x /= div_factor;
        });

        Ok(scores)
    }

    fn assign_scores(&mut self) -> Result<(), M::Error> {
        let scores = self.get_scores()?;
        for (v, s) in self.data.iter_mut().zip(scores.iter()) {
            v.assign_score(*s);
        }
        Ok(())
    }

    pub fn score(self) -> Vec<T> {
        self.data
    }

    /// Read-only access to the scored items (e.g. to query feature names).
    pub fn data(&self) -> &[T] {
        &self.data
    }

    /// Compute per-fold feature means + per-fold model importance.
    ///
    /// The stats for fold `f` summarize `fold_rows[f]` — for THIS partition,
    /// the rows fold `f`'s model was TRAINED on. Folds with no fitted model
    /// (shouldn't happen post-fit) produce empty maps. See
    /// [`fold_feature_stats`], which the LDA cross-fit shares with a different
    /// partition.
    pub fn feature_stats(&self) -> RescoreFeatureStats {
        let models: Vec<Option<&M>> = self.fold_models.iter().map(|m| m.as_ref()).collect();
        fold_feature_stats(&self.dataset, &self.fold_rows, &models)
    }

    /// THE fold partition this scorer built: `fold_rows[f]` is the rows with
    /// `get_fold(i) == f`, ascending.
    ///
    /// Test-only, and it exists for one assertion. The hybrid rescorers feed
    /// this scorer a column that was cross-fit under a DIFFERENT train/score
    /// split (`qvalues::crossfit`), and the two splits are only leak-free
    /// together if they agree on the fold ASSIGNMENT — a disagreement is
    /// invisible to this scorer's own cross-validation and shows up solely as an
    /// optimistic FDR. Exposing the partition lets
    /// `qvalues::crossfit_holds_out_exactly_the_rows_the_gbm_scorer_trains_on`
    /// compare the two against each other rather than each against a repeated
    /// `i % n_folds`.
    #[cfg(test)]
    pub(crate) fn fold_rows(&self) -> &[Vec<usize>] {
        &self.fold_rows
    }

    fn next_fold(&self, fold: u8) -> u8 {
        let mut maybe_next = fold + 1;
        if maybe_next >= self.n_folds {
            maybe_next -= self.n_folds;
        }
        assert!(maybe_next < self.n_folds);
        maybe_next
    }

    /// Fit fold `fold`: trained on that fold's rows, early-stopped on the next
    /// fold's. Every other fold is therefore safe to score with it.
    fn fit_fold(&mut self, fold: u8) -> Result<(), M::Error> {
        let val_fold = self.next_fold(fold);
        let model = M::fit(
            &self.config,
            &self.dataset,
            fold as usize,
            &self.fold_rows[fold as usize],
            &self.fold_rows[val_fold as usize],
        )?;
        self.fold_models[fold as usize] = Some(model);
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::distr::{
        Distribution,
        Uniform,
    };

    use crate::ml::cv::{
        FeatureLike,
        FeatureVector,
    };

    struct MyFeature {
        vals: [f64; 5],
        class: f64,
        score: f64,
    }

    impl FeatureVector for MyFeature {
        fn as_feature(&self) -> impl IntoIterator<Item = f64> + '_ {
            self.vals
        }
    }

    impl FeatureLike for MyFeature {
        fn get_y(&self) -> f64 {
            self.class
        }

        fn assign_score(&mut self, score: f64) {
            self.score = score;
        }

        fn get_score(&self) -> f64 {
            self.score
        }
    }

    /// Seeded so the test is reproducible: with `rand::rng()` the thresholds
    /// below had to absorb run-to-run sampling noise, which is exactly what
    /// made them uninformative.
    fn random_data(num_targets: usize, num_decoys: usize, seed: u64) -> Vec<MyFeature> {
        use rand::SeedableRng;
        let between_unch = Uniform::try_from(1.0..10.0).unwrap();
        let between_targ = Uniform::try_from(1.0..20.0).unwrap();
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
        let mut out = Vec::new();

        for _nt in 0..num_targets {
            let tmp = MyFeature {
                vals: [
                    between_unch.sample(&mut rng),
                    between_unch.sample(&mut rng),
                    between_targ.sample(&mut rng),
                    between_targ.sample(&mut rng),
                    between_targ.sample(&mut rng),
                ],
                class: 1.0,
                score: f64::NAN,
            };
            out.push(tmp);
        }
        for _nt in 0..num_decoys {
            let tmp = MyFeature {
                vals: [
                    between_unch.sample(&mut rng),
                    between_unch.sample(&mut rng),
                    between_unch.sample(&mut rng),
                    between_unch.sample(&mut rng),
                    between_unch.sample(&mut rng),
                ],
                class: 0.0,
                score: f64::NAN,
            };
            out.push(tmp);
        }
        out
    }

    #[test]
    fn test_cv_train() {
        let config = GBMConfig::default();
        let data = random_data(500, 500, 0xC0FFEE);
        let data_len = data.len();

        let mut scorer =
            CrossValidatedScorer::<MyFeature, GbmFoldModel>::new_from_shuffled(3, data, config);
        scorer.fit().unwrap();

        // Rows 0..500 are the targets, 500..1000 the decoys.
        let out = scorer.get_scores().unwrap();
        let avg_t: f64 = out[..500].iter().sum::<f64>() / 500.0;
        let avg_d: f64 = out[500..].iter().sum::<f64>() / 500.0;
        let num_t_gt_d = out[..500].iter().filter(|&&x| x > avg_d).count();
        let num_d_gt_t = out[500..].iter().filter(|&&x| x > avg_t).count();
        // 3 of the 5 features are drawn from 2x the decoy uniform range, so a
        // target usually lands outside the decoy support on at least one of
        // them. With the RNG seeded these are exact-run facts rather than
        // distributional guesses: this seed yields 470 and 0. Thresholds are
        // tightened from the old 300 / 100 (which had to absorb unseeded
        // sampling noise) but kept off the observed values so an unrelated GBM
        // tweak still has slack.
        assert!(
            num_t_gt_d > 440,
            "num_t_gt_d: {}, averages={} {}",
            num_t_gt_d,
            avg_t,
            avg_d
        );
        assert!(
            num_d_gt_t < 25,
            "num_d_gt_t: {}, averages {} and {}",
            num_d_gt_t,
            avg_t,
            avg_d
        );
        assert_eq!(out.len(), data_len);
    }

    /// A [`FoldModel`] that reports whatever importance vector it is configured
    /// with, so the sentinel contract can be tested without a real fit.
    struct StubModel(Vec<f32>);

    impl FoldModel for StubModel {
        type Config = Vec<f32>;
        type Error = ();

        fn fit<D: FoldDataset>(
            cfg: &Vec<f32>,
            _data: &D,
            _fold: usize,
            _train: &[usize],
            _val: &[usize],
        ) -> Result<Self, ()> {
            Ok(StubModel(cfg.clone()))
        }

        fn predict<D: FoldDataset>(&self, _data: &D, rows: &[usize]) -> Result<Vec<f64>, ()> {
            Ok(vec![0.0; rows.len()])
        }

        fn importance(&self) -> Vec<f32> {
            self.0.clone()
        }
    }

    /// THE [`FoldModel::importance`] sentinel contract, at the sidecar boundary:
    /// `NAN` means "this model reports nothing for this column" and is dropped;
    /// every FINITE value reaches the sidecar, `0.0` included.
    ///
    /// Regression guard. The boundary used to filter `!= 0.0`, which was right
    /// for a tree model (where 0.0 did mean "never split on") and wrong for
    /// every other model: an LDA's `|coef|` of exactly 0.0 is a measurement of a
    /// dead or constant column, and silently deleting those rows removes
    /// exactly what an operator reads the sidecar to find. Both halves are
    /// asserted — dropping the NAN alone would also pass if 0.0 were dropped
    /// too, so the surviving zero is the load-bearing assertion.
    #[test]
    fn importance_nan_is_unreported_but_zero_is_a_value() {
        let dataset = RowMajorDataset::new(
            PrecomputedFeatures::from_row_major(
                vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
                3,
                vec![1.0, 0.0],
            ),
            vec![Arc::from("unreported"), Arc::from("zero"), Arc::from("big")],
            1,
        );
        let model = StubModel(vec![f32::NAN, 0.0, 2.5]);
        let stats = fold_feature_stats(&dataset, &[vec![0, 1]], &[Some(&model)]);

        let imp = &stats[0].feature_importance;
        // Gain-descending, and the genuine zero SURVIVES.
        assert_eq!(
            imp.iter().map(|(n, v)| (&**n, *v)).collect::<Vec<_>>(),
            vec![("big", 2.5), ("zero", 0.0)],
            "a finite 0.0 is a reported measurement and must reach the sidecar"
        );
        // The NAN does not.
        assert!(
            !imp.iter().any(|(n, _)| &**n == "unreported"),
            "a NAN importance means 'not reported' and must be absent: {imp:?}"
        );

        // The stats half is unaffected by the sentinel: it is always full width.
        assert_eq!(stats[0].feature_stats.len(), 3);
    }

    /// The GBM side of the same contract: forust reports only the columns it
    /// split on, and the ones it did not must come back as `NAN` (absent from
    /// the sidecar), NOT as a 0.0 gain it never measured.
    ///
    /// `feature_4` is constant, so no tree can split on it. The other four are
    /// the usual separable draws, so their presence keeps the absence
    /// assertion from passing vacuously.
    #[test]
    fn gbm_importance_omits_columns_it_never_split_on() {
        let mut data = random_data(300, 300, 0xBEEF);
        for d in data.iter_mut() {
            d.vals[4] = 1.0;
        }

        let mut scorer = CrossValidatedScorer::<MyFeature, GbmFoldModel>::new_from_shuffled(
            3,
            data,
            GBMConfig::default(),
        );
        scorer.fit().unwrap();

        for fold in scorer.feature_stats() {
            let names: Vec<&str> = fold
                .feature_importance
                .iter()
                .map(|(n, _)| &**n)
                .collect::<Vec<_>>();
            assert!(
                !names.contains(&"feature_4"),
                "fold {}: constant column was never split on, so it must be \
                 unreported rather than a fabricated 0.0 gain; got {names:?}",
                fold.fold
            );
            assert!(
                !names.is_empty(),
                "fold {}: no feature reported at all, so the check above is vacuous",
                fold.fold
            );
            assert!(
                fold.feature_importance.iter().all(|(_, g)| g.is_finite()),
                "fold {}: the sidecar must never carry the NAN sentinel itself",
                fold.fold
            );
        }
    }
}
