//! [`FoldModel`] adapter over the dependency-free MLP in [`crate::ml::mlp`].
//!
//! All the numerics live in `mlp.rs`; this module is the wiring between them
//! and the cross-fitting traits in `cv.rs`. Its whole job is to keep two
//! properties true:
//!
//!  1. **Every fitted statistic is train-fold-only.** The [`ColumnTransform`]
//!     (cull set, standardization moments, imputation means) is fitted from the
//!     `train` row slice, stored on the model, and applied UNCHANGED to whatever
//!     rows [`FoldModel::predict`] is later handed. Fitting it over all rows, or
//!     refitting it inside `predict` from the batch being scored, would use
//!     held-out feature values — the cross-fit's leak-freedom is a property of
//!     the caller's partition, and this impl must not undermine it.
//!  2. **The fit is a pure function of the config, the fold index and the
//!     row slices.** The RNG comes from [`MlpConfig::rng_for_fold`], i.e. from
//!     the configured seed and the `fold` argument [`FoldModel::fit`] is handed
//!     and nothing else — no call counter, no clock, no hash-map iteration
//!     order — so folds differ from one another while the same build on the
//!     same input produces bit-identical scores. Early stopping does not weaken
//!     this: the stopping decision is a strict comparison of floats computed by
//!     the same deterministic forward pass, and the inner-validation carve
//!     ([`inner_val_split`]) is a pure function of `train.len()`.
//!  3. **Reporting cannot move the fit.** Two feature-importance metrics are
//!     selectable ([`MlpImportance`]) and the choice reaches no weight, no RNG
//!     draw and no optimizer step: the gradient metric is a DEDICATED
//!     forward/backward pass that runs after training and wipes the parameter
//!     gradients it accumulated, so two runs differing only in
//!     `TIMSSEEK_MLP_IMPORTANCE` produce bit-identical scores and differ only in
//!     the sidecar. That is what makes diffing the two rankings meaningful.

use crate::ml::cv::{
    FoldDataset,
    FoldModel,
};
use crate::ml::mlp::{
    Adam,
    ColumnTransform,
    Mlp,
    MlpBatch,
    MlpConfig,
    MlpImportance,
    MlpRuntimeProfile,
    Tensor,
    TrainOutcome,
};
use std::cell::RefCell;
use std::sync::Arc;
use std::time::{
    Duration,
    Instant,
};

/// The ways fitting or scoring an [`MlpFoldModel`] can fail.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MlpFoldError {
    /// `predict` was handed a dataset whose lane width disagrees with the one
    /// the transform was fitted against. Scoring anyway would silently read
    /// column `j` of a different matrix.
    WidthMismatch { fitted: usize, got: usize },
    /// The cull left nothing to train on: every lane column was entirely
    /// non-finite on the train rows, or measured a std at or below `MIN_STD`
    /// (see [`ColumnTransform`] for what that does and does not guarantee).
    ///
    /// Carries its own diagnostic context because the callers that report it
    /// cannot reconstruct it: [`FoldModel::fit`]'s error travels up through
    /// `CrossValidatedScorer::fit`, which knows nothing about columns, and the
    /// abort in `qvalues::rescore_mlp_with` prints only what the error says.
    ///
    /// `train_rows == 0` is a DIFFERENT failure wearing the same variant: with
    /// no rows every column is vacuously non-finite, so the cull takes all of
    /// them and the operator's actual problem is "this fold got zero rows", not
    /// "the lane is dead". [`Display`](std::fmt::Display) separates the two —
    /// the message used to blame the cull in both cases, which sent the reader
    /// looking at features when the fixture or the fold partition was the
    /// problem.
    NoUsableColumns {
        /// The fold index this fit was for.
        fold: usize,
        /// Lane width, i.e. how many columns were offered and therefore culled.
        ncols: usize,
        /// How many rows the cull was fitted over. Zero means the train slice
        /// was empty.
        train_rows: usize,
    },
    /// `predict` produced a non-finite logit.
    ///
    /// The MLP is the ONLY rescorer whose discriminant is not finite by
    /// construction: the LDA rejects a non-finite `coef` at fit time
    /// (`lda::LdaModel` rejects a non-finite solve), and a forust prediction is a sum of leaf
    /// values. Here the score is whatever the net's output layer computes, and
    /// weights that diverged during training (or overflowed an activation) put a
    /// `NaN` or an infinity in it.
    ///
    /// Caught HERE rather than downstream because downstream diagnoses it wrongly.
    /// A non-finite score reaches `qvalues::finalize` -> `assign_qval`, whose
    /// release-mode `assert!` says "Expecting scores to be sorted in descending
    /// order" — a `NaN` makes every comparison in that window false, so the run
    /// aborts several stages later blaming the SORT for a model failure. Reporting
    /// it as an `MlpFoldError` puts the abort next to the model that produced it
    /// and routes it through the diagnostics `qvalues::abort_standalone_mlp`
    /// already has.
    NonFiniteScore {
        /// Dataset row index of the FIRST non-finite score, so the caller can
        /// look at that row.
        row: usize,
        /// How many of the scored rows came out non-finite. One is a suspicious
        /// row; all of them is a diverged net.
        n_bad: usize,
        /// How many rows were in this batch.
        scored: usize,
    },
}

impl std::fmt::Display for MlpFoldError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MlpFoldError::WidthMismatch { fitted, got } => {
                write!(
                    f,
                    "dataset has {got} lane columns, model was fitted on {fitted}"
                )
            }
            MlpFoldError::NoUsableColumns {
                fold,
                ncols,
                train_rows: 0,
            } => write!(
                f,
                "fold {fold}: the train slice is EMPTY, so all {ncols} lane columns are \
                 vacuously dead — this is a zero-row fold, not a dead feature set"
            ),
            MlpFoldError::NoUsableColumns {
                fold,
                ncols,
                train_rows,
            } => write!(
                f,
                "fold {fold}: all {ncols} lane columns were culled (each one non-finite or \
                 constant across the {train_rows} train rows)"
            ),
            MlpFoldError::NonFiniteScore { row, n_bad, scored } => write!(
                f,
                "the net produced a NON-FINITE logit for {n_bad} of {scored} scored rows (first \
                 at dataset row {row}) — the weights or an activation diverged during training"
            ),
        }
    }
}

impl std::error::Error for MlpFoldError {}

#[cfg(test)]
fn gather<D: FoldDataset>(data: &D, rows: &[usize], ncols: usize) -> Vec<f64> {
    let mut feat = vec![0.0f64; rows.len() * ncols];
    for (k, &i) in rows.iter().enumerate() {
        data.get_values(i, &mut feat[k * ncols..(k + 1) * ncols]);
    }
    feat
}

/// A trained MLP plus THE input transform it was trained through.
///
/// The two are one unit on purpose: the net's weights are only meaningful on
/// rows standardized by this exact transform, so they are stored together and
/// `predict` has no way to reach a differently-fitted one.
pub struct MlpFoldModel {
    /// [`Mlp::forward`] writes into the model's own activation buffers, so it
    /// needs `&mut`, while [`FoldModel::predict`] takes `&self`. A `RefCell`
    /// rather than making the buffers external: `Mlp` is already `!Sync` (it
    /// holds `Box<dyn Layer>`), so this gives up nothing.
    net: RefCell<Mlp>,
    /// Fitted on the TRAIN rows at `fit` time and never refitted. This is the
    /// field the whole leak argument rests on.
    transform: ColumnTransform,
    /// Mean training loss of the epoch whose weights the net holds. Diagnostics
    /// only — `fit` logs it, and [`MlpFoldModel::final_loss`] hands it to a
    /// failing test's message.
    #[cfg_attr(not(test), allow(dead_code))]
    final_loss: f32,
    /// Rows that actually reached the optimizer: `train.len()` minus whatever
    /// the inner-validation carve took. The ONE externally checkable statement
    /// that the carve is excluded from training — see
    /// `inner_validation_rows_are_withheld_from_the_optimizer`.
    #[cfg_attr(not(test), allow(dead_code))]
    n_train_rows: usize,
    /// What the training run did — epochs run, epoch kept, whether it rolled
    /// back. `fit` logs it; the tests assert on it.
    #[cfg_attr(not(test), allow(dead_code))]
    outcome: TrainOutcome,
    /// Which metric [`FoldModel::importance`] reports. Copied off the config at
    /// `fit` time because `importance(&self)` has no config in scope, and because
    /// the `Grad` metric has to be MEASURED during `fit` — the train rows are gone
    /// by the time anyone asks.
    importance_kind: MlpImportance,
    /// Mean `|dL/dx|` per TRANSFORMED input over the train rows, from the
    /// dedicated post-training pass in [`Mlp::input_grad_importance`]. `None` when
    /// nothing asked for it: the pass costs about one epoch, so it runs only when
    /// [`MlpConfig::importance`] selects it or when the `debug` comparison log is
    /// enabled (the same `tracing::enabled!` gate the held-out trace uses).
    grad_importance: Option<Vec<f32>>,
    /// Bound inference storage to the same row count used for train chunks.
    inference_batch_size: usize,
}

// Both accessors are TEST SURFACE and are scoped to it. Nothing in the live path
// reads either one, but the leak property this module exists to keep
// (`transform` is fitted on the train rows and never refitted) is only checkable
// from outside through `transform()`, and `final_loss()` is what makes a failing
// convergence assertion diagnosable. `pub` would have advertised them as API.
#[cfg(test)]
impl MlpFoldModel {
    /// The transform fitted at `fit` time — the one `predict` uses.
    pub(crate) fn transform(&self) -> &ColumnTransform {
        &self.transform
    }

    /// Mean loss over the final training epoch.
    pub(crate) fn final_loss(&self) -> f32 {
        self.final_loss
    }

    /// Rows the optimizer actually saw.
    pub(crate) fn n_train_rows(&self) -> usize {
        self.n_train_rows
    }

    /// What the training run did.
    pub(crate) fn outcome(&self) -> TrainOutcome {
        self.outcome
    }
}

/// Map a per-TRANSFORMED-INPUT vector onto lane columns, which is what
/// [`FoldModel::importance`] has to return.
///
/// Shared by both importance metrics on purpose: the two differ in what they
/// measure, never in how a measurement is placed, and a second copy of this loop
/// is how the `_isna` folding rule would come to disagree with itself.
///
/// **Culled lanes come back `0.0`, not `NAN`.** The output starts zeroed and a
/// culled lane is never written, which is the [`FoldModel::importance`] contract's
/// "a culled column genuinely contributed nothing" case: a finite `0.0` reaches
/// the sidecar and tells the operator the column was culled, where `NAN` would be
/// filtered out at the `fold_feature_stats` boundary and delete the cull report.
///
/// **`_isna` companions have no lane of their own**, so each one's value is summed
/// into the column it flags. For either metric the resulting number answers one
/// question — how much this column matters, through its value AND through its
/// missingness — which is the only reading that makes the two comparable.
fn fold_inputs_to_lanes(transform: &ColumnTransform, raw: &[f32]) -> Vec<f32> {
    // Every transformed input is either a standardized column or a companion, so
    // with the right length the loop below places all of them; a short/long
    // vector is the only way one could go missing.
    debug_assert_eq!(
        raw.len(),
        transform.width(),
        "an importance vector must report one value per transformed input"
    );
    let mut out = vec![0.0f32; transform.ncols_lane()];
    for (k, &v) in raw.iter().enumerate() {
        let lane = transform
            .lane_of_input(k)
            .or_else(|| transform.isna_lane_of_input(k));
        if let Some(lane) = lane {
            out[lane] += v;
        }
    }
    debug_assert!(
        out.iter().all(|v| v.is_finite()),
        "a non-finite importance would read as the NAN 'unreported' sentinel"
    );
    debug_assert!(
        {
            let (a, b) = (out.iter().sum::<f32>(), raw.iter().sum::<f32>());
            (a - b).abs() <= 1e-3 * b.abs().max(1.0)
        },
        "folding companions into their parent lane must not drop or double any value"
    );
    out
}

/// How many lane columns the per-fold `debug` comparison line names per metric.
///
/// The line exists so the two rankings can be eyeballed WITHOUT a second search
/// (requirement: diffable in one run). It is a top-N and not the whole vector
/// because the vector is 101 or 128 columns wide and three folds of that at
/// `debug` is a wall of text; the authoritative diff is still two runs, one per
/// `TIMSSEEK_MLP_IMPORTANCE` value, whose sidecars are directly comparable.
const IMPORTANCE_LOG_TOP_N: usize = 10;

/// `name=value` for the top `IMPORTANCE_LOG_TOP_N` lanes, descending.
///
/// Ties keep lane order (the sort is stable), matching how `fold_feature_stats`
/// orders the sidecar, so the log and the file agree on which of two equal
/// columns comes first.
fn top_lanes(lane_imp: &[f32], names: &[Arc<str>], n: usize) -> String {
    let mut idx: Vec<usize> = (0..lane_imp.len()).collect();
    idx.sort_by(|&a, &b| {
        lane_imp[b]
            .partial_cmp(&lane_imp[a])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    idx.into_iter()
        .take(n)
        .map(|j| format!("{}={:.6}", names[j], lane_imp[j]))
        .collect::<Vec<_>>()
        .join(" ")
}

/// One row in every `INNER_VAL_STRIDE` of `train` is carved off as the inner
/// validation slice — 20%.
///
/// Why a STRIDE and not a tail cut or a random sample:
///
///  * **Deterministic with no RNG.** It is a pure function of `train.len()`, so
///    it cannot consume a draw and shift the weight init or the shuffle
///    sequence, and two runs of the same build carve the same rows by
///    construction. A random sample would have to be seeded from somewhere, and
///    every seed source available here is already spoken for.
///  * **Class-proportional under any row ORDER.** The rescore path shuffles
///    before building the matrix, so a tail cut would usually be fine there —
///    but `FoldModel::fit` is handed whatever slice the caller has, and the
///    obvious unshuffled layout is targets-then-decoys, where a tail cut is a
///    single-class validation set and the stopping decision becomes nonsense. A
///    stride takes 20% of every contiguous run of rows, whatever it holds.
const INNER_VAL_STRIDE: usize = 5;

/// Below this many train rows, no inner slice is carved and no early stopping
/// happens on the `val = &[]` path — the full epoch budget runs, as it did
/// before early stopping existed.
///
/// `160` is `MIN_INNER_VAL_ROWS * INNER_VAL_STRIDE`: the floor is really on the
/// SIZE OF THE VALIDATION SET, 32 rows. A held-out loss averaged over fewer
/// rows than that is dominated by which rows landed in it, and a patience rule
/// reading it stops on sampling noise rather than on the overfitting turn —
/// worse than not stopping at all. The 20% also costs training rows, which a
/// small fold cannot spare.
///
/// This is a floor chosen for the reason above, not a tuned value. Production
/// folds are 10^4-10^5 rows and are nowhere near it; the fixtures in this crate
/// straddle it, so some of them early-stop and some do not.
const MIN_TRAIN_ROWS_FOR_INNER_VAL: usize = MIN_INNER_VAL_ROWS * INNER_VAL_STRIDE;

/// The floor on the number of rows the stopping decision is averaged over.
///
/// It applies to EITHER source of those rows — the inner carve (via
/// [`MIN_TRAIN_ROWS_FOR_INNER_VAL`], which is this times the stride) and the
/// outer `val` slice a caller passes (via [`usable_val`]). The argument above is
/// about the SIZE of the held-out set and says nothing about where it came from:
/// a 30-row outer slice makes exactly the same noisy stopping decision a 30-row
/// carve would, and an outer slice is not privileged just because the caller
/// chose it. Under the floor the knob goes inert and the full epoch budget runs,
/// which is what the pre-early-stopping behaviour was.
const MIN_INNER_VAL_ROWS: usize = 32;

/// The outer `val` slice if it is big enough to make a stopping decision on,
/// otherwise nothing — see [`MIN_INNER_VAL_ROWS`].
///
/// Returning an empty slice rather than a flag so the rest of `fit` has ONE
/// notion of "is there a held-out set": an under-floor slice is indistinguishable
/// from the `&[]` the `crossfit` driver passes, and therefore takes the same
/// path (an inner carve if `train` can afford one, no early stopping if it
/// cannot).
fn usable_val(val: &[usize]) -> &[usize] {
    if val.len() < MIN_INNER_VAL_ROWS {
        &[]
    } else {
        val
    }
}

/// Split `0..len` into (rows to TRAIN on, rows to make the stopping decision
/// on), or `(all, none)` when `len` is under [`MIN_TRAIN_ROWS_FOR_INNER_VAL`].
///
/// Positions, not dataset row indices: both halves index the caller's `train`
/// slice, whose position `k` names dataset row `train[k]`.
///
/// Offset `INNER_VAL_STRIDE - 1` rather than `0` so position 0 trains; nothing
/// depends on which offset it is, only on it being fixed.
fn inner_val_split(len: usize) -> (Vec<usize>, Vec<usize>) {
    if len < MIN_TRAIN_ROWS_FOR_INNER_VAL {
        return ((0..len).collect(), Vec::new());
    }
    (0..len).partition(|k| k % INNER_VAL_STRIDE != INNER_VAL_STRIDE - 1)
}

fn batch_count(nrows: usize, batch_size: usize) -> usize {
    nrows.div_ceil(batch_size.max(1))
}

#[derive(Default)]
struct BatchVisitTimings {
    prepare_cpu: Duration,
    consumer_wait: Duration,
    batches: u64,
}

impl std::ops::AddAssign for BatchVisitTimings {
    fn add_assign(&mut self, rhs: Self) {
        self.prepare_cpu += rhs.prepare_cpu;
        self.consumer_wait += rhs.consumer_wait;
        self.batches += rhs.batches;
    }
}

fn runtime_profile_enabled() -> bool {
    std::env::var("TIMSSEEK_MLP_PROFILE").is_ok_and(|raw| {
        !matches!(
            raw.trim().to_ascii_lowercase().as_str(),
            "" | "0" | "false" | "off"
        )
    })
}

/// Visit transformed batches in `chunk_order`, with exactly two reusable buffers.
/// A scoped producer projects the next chunk while the consumer performs neural
/// work on the current one. No transformed row survives the visit.
fn visit_prefetched_batches<D: FoldDataset>(
    data: &D,
    rows: &[usize],
    transform: &ColumnTransform,
    batch_size: usize,
    chunk_order: &[usize],
    consume: &mut dyn FnMut(&MlpBatch),
    profile: bool,
) -> BatchVisitTimings {
    let ncols = data.column_names().len();
    let width = transform.width();
    let batch_size = batch_size.max(1);
    debug_assert!(
        chunk_order
            .iter()
            .all(|&chunk| chunk < batch_count(rows.len(), batch_size))
    );

    std::thread::scope(|scope| {
        let (empty_tx, empty_rx) = std::sync::mpsc::sync_channel(2);
        let (ready_tx, ready_rx) = std::sync::mpsc::sync_channel(2);
        empty_tx.send(MlpBatch::buffer(batch_size, width)).unwrap();
        empty_tx.send(MlpBatch::buffer(batch_size, width)).unwrap();

        let producer = scope.spawn(move || {
            let mut raw = vec![0.0f64; ncols];
            let mut prepare_cpu = Duration::ZERO;
            for &chunk_id in chunk_order {
                let mut batch = empty_rx.recv().expect("MLP batch consumer disappeared");
                let started = profile.then(Instant::now);
                let start = chunk_id * batch_size;
                let end = (start + batch_size).min(rows.len());
                let source = &rows[start..end];
                batch.x.reshape(source.len(), width);
                batch.y.clear();
                batch.w.clear();
                for (k, &i) in source.iter().enumerate() {
                    data.get_values(i, &mut raw);
                    transform.apply(&raw, batch.x.row_mut(k));
                    let decoy = data.is_decoy(i);
                    batch.y.push(if decoy { 0.0 } else { 1.0 });
                    batch.w.push(if decoy { 1.0 } else { 0.5 });
                }
                if let Some(started) = started {
                    prepare_cpu += started.elapsed();
                }
                ready_tx
                    .send(batch)
                    .expect("MLP batch consumer disappeared");
            }
            prepare_cpu
        });

        let mut consumer_wait = Duration::ZERO;
        for _ in chunk_order {
            let started = profile.then(Instant::now);
            let batch = ready_rx.recv().expect("MLP batch producer disappeared");
            if let Some(started) = started {
                consumer_wait += started.elapsed();
            }
            consume(&batch);
            let _ = empty_tx.send(batch);
        }
        BatchVisitTimings {
            prepare_cpu: producer.join().expect("MLP batch producer panicked"),
            consumer_wait,
            batches: chunk_order.len() as u64,
        }
    })
}

impl FoldModel for MlpFoldModel {
    type Config = MlpConfig;
    type Error = MlpFoldError;

    fn fit_progress_units(cfg: &MlpConfig) -> u64 {
        cfg.epochs as u64
    }

    fn fit<D: FoldDataset>(
        cfg: &MlpConfig,
        data: &D,
        fold: usize,
        train: &[usize],
        val: &[usize],
    ) -> Result<Self, MlpFoldError> {
        Self::fit_with_progress(cfg, data, fold, train, val, &|| {})
    }

    /// `val` is the EARLY-STOPPING SLICE, and it is used as one whenever
    /// [`MlpConfig::early_stopping_patience`] is set. It reaches
    /// [`ColumnTransform::apply`] and [`Mlp::eval_loss`] and nothing else — not
    /// the transform FIT, not the optimizer, not the RNG — so the only thing it
    /// can move is which epoch's weights are kept, which is what early stopping
    /// means.
    ///
    /// # The `val = &[]` callers
    /// A caller that has no fold to spare may pass an empty slice. Rather than
    /// silently training the full budget there, an inner slice is carved out of `train` — see
    /// [`inner_val_split`] for the rule and [`MIN_TRAIN_ROWS_FOR_INNER_VAL`] for
    /// when it declines to. Those rows are excluded from the optimizer and used
    /// only for the stopping decision.
    ///
    /// An outer `val` slice SHORTER THAN [`MIN_INNER_VAL_ROWS`] takes that same
    /// path, i.e. it is ignored: the floor is on the size of the held-out set,
    /// and a 30-row outer slice makes the same sampling-noise stopping decision a
    /// 30-row carve would. See [`usable_val`].
    ///
    /// **The [`ColumnTransform`] is still fitted over ALL of `train`,
    /// inner-validation rows included.** Deliberate, and the choice worth
    /// stating: those rows are train-fold rows — the outer partition does not
    /// hold them out from this model, this function does, for a decision that
    /// never leaves it — so including them breaks no leak property. What it buys
    /// is that the cull set, the standardization moments and the imputation
    /// means do not move when the patience knob does: a config knob that changed
    /// WHICH LANE COLUMNS ARE ALIVE would change the importance sidecar and the
    /// [`MlpFoldError::NoUsableColumns`] failure mode along with it. The cost is
    /// that the inner-validation loss is very slightly optimistic, through
    /// label-free per-column moments that 80% of the same rows already
    /// determine — second-order next to the epoch it is being used to choose.
    ///
    /// `fold` is NOT ignored: it is the only thing that distinguishes one
    /// fold's initialization from another's, via
    /// [`MlpConfig::rng_for_fold`]. Passing the same `fold` for every fold
    /// would give every fold the same weight init and the same shuffle order —
    /// still leak-free, but the folds would be correlated in a way nothing
    /// downstream would notice.
    ///
    /// Everything fitted here — the cull set, the standardization moments, the
    /// imputation means, the weights — comes from `train` and only `train`.
    fn fit_with_progress<D: FoldDataset>(
        cfg: &MlpConfig,
        data: &D,
        fold: usize,
        train: &[usize],
        val: &[usize],
        epoch_finished: &(dyn Fn() + Sync),
    ) -> Result<Self, MlpFoldError> {
        let profile_enabled = runtime_profile_enabled();
        let fit_started = profile_enabled.then(Instant::now);
        let names = data.column_names();
        let ncols = names.len();
        let transform_started = profile_enabled.then(Instant::now);
        let transform = ColumnTransform::fit_streaming(ncols, train.iter().copied(), |i, out| {
            data.get_values(i, out)
        });
        let transform_wall = transform_started.map_or(Duration::ZERO, |t| t.elapsed());

        if !transform.culled().is_empty() {
            let culled: Vec<&str> = transform.culled().iter().map(|&j| &*names[j]).collect();
            tracing::warn!(
                "MLP fold {}: culled {}/{} lane columns (all non-finite or constant on the \
                 {} train rows): {:?}",
                fold,
                culled.len(),
                ncols,
                train.len(),
                culled,
            );
        }

        let width = transform.width();
        if width == 0 {
            return Err(MlpFoldError::NoUsableColumns {
                fold,
                ncols,
                train_rows: train.len(),
            });
        }

        // Where the stopping decision gets measured. A big enough outer `val`
        // slice is preferred whenever there is one — it costs no training rows —
        // and the inner carve is the fallback for the callers that pass `&[]`
        // (and for an outer slice under [`MIN_INNER_VAL_ROWS`], which is not a
        // set worth stopping on; see [`usable_val`]). With early stopping off
        // there is no decision to make, so no carve either and every train row
        // trains, exactly as before.
        let early_stopping = cfg.early_stopping_patience.is_some();
        let val_rows = usable_val(val);
        if !val.is_empty() && val_rows.is_empty() {
            tracing::debug!(
                "MLP fold {}: the outer validation slice has {} rows, under the {} a held-out \
                 measurement needs — ignoring it, so this fit is exactly the `val = &[]` fit. \
                 That falls back to the inner carve ONLY with early stopping on and enough \
                 train rows to carve from (too few rows logs just below); with early stopping \
                 off there is no carve and no stopping decision to make",
                fold,
                val.len(),
                MIN_INNER_VAL_ROWS,
            );
        }
        let (fit_pos, inner_pos) = if early_stopping && val_rows.is_empty() {
            inner_val_split(train.len())
        } else {
            ((0..train.len()).collect(), Vec::new())
        };
        if early_stopping && val_rows.is_empty() && inner_pos.is_empty() {
            tracing::debug!(
                "MLP fold {}: {} train rows is under the {} needed to carve an inner \
                 validation slice — training the full {} epoch budget with no early stopping",
                fold,
                train.len(),
                MIN_TRAIN_ROWS_FOR_INNER_VAL,
                cfg.epochs,
            );
        }

        // Labels are target=1 / decoy=0 and weights are target=0.5 / decoy=1.0,
        // matching the GBM objective. Only row indices are retained here; feature
        // values are regenerated into the two-buffer loader on every pass.
        let fit_rows: Vec<usize> = fit_pos.iter().map(|&p| train[p]).collect();

        // The held-out set, transformed through the TRAIN-fitted transform. It
        // reaches `Mlp::eval_loss` and nothing else — no optimizer step, no
        // transform fit, no RNG draw.
        //
        // Streamed whenever early stopping is on: it is the stopping decision,
        // not just a log line. With early stopping off it is visited only for the
        // trace, as before.
        let want_val = early_stopping || tracing::enabled!(tracing::Level::DEBUG);
        let inner_rows: Vec<usize> = inner_pos.iter().map(|&p| train[p]).collect();
        let validation_rows = if !want_val {
            None
        } else if !val_rows.is_empty() {
            Some(val_rows)
        } else if !inner_pos.is_empty() {
            Some(inner_rows.as_slice())
        } else {
            None
        };

        // Seeded from (config seed, fold index) alone — see the module docs.
        let mut rng = cfg.rng_for_fold(fold);
        let mut net = Mlp::feedforward(width, &cfg.hidden, cfg.activation, &mut rng);
        let mut opt = Adam::new(cfg.lr).with_weight_decay(cfg.weight_decay);
        let tag = format!("MLP fold {fold}: ");
        let n_train_batches = batch_count(fit_rows.len(), cfg.batch_size);
        let mut train_loader = BatchVisitTimings::default();
        let mut visit_train = |order: &[usize], consume: &mut dyn FnMut(&MlpBatch)| {
            train_loader += visit_prefetched_batches(
                data,
                &fit_rows,
                &transform,
                cfg.batch_size,
                order,
                consume,
                profile_enabled,
            );
        };
        let val_order: Vec<usize> = validation_rows
            .map(|rows| (0..batch_count(rows.len(), cfg.batch_size)).collect())
            .unwrap_or_default();
        let mut validation_loader = BatchVisitTimings::default();
        let mut visit_val = |consume: &mut dyn FnMut(&MlpBatch)| {
            validation_loader += visit_prefetched_batches(
                data,
                validation_rows.expect("validation callback requires rows"),
                &transform,
                cfg.batch_size,
                &val_order,
                consume,
                profile_enabled,
            );
        };
        #[allow(clippy::type_complexity)]
        let val_callback: Option<&mut dyn FnMut(&mut dyn FnMut(&MlpBatch))> = validation_rows
            .is_some()
            .then_some(&mut visit_val as &mut dyn FnMut(&mut dyn FnMut(&MlpBatch)));
        let mut runtime_profile = MlpRuntimeProfile::default();
        let outcome = net.train_reporting_streaming(
            cfg,
            n_train_batches,
            &mut visit_train,
            val_callback,
            &mut opt,
            &mut rng,
            &tag,
            epoch_finished,
            profile_enabled.then_some(&mut runtime_profile),
        );

        // --- Feature importance, measured now because the rows are only here ---
        //
        // A DEDICATED PASS AFTER TRAINING, over the TRAIN rows only. Three
        // properties, all load-bearing:
        //
        //  * It measures THE FITTED NET. `train_reporting` has already rolled the
        //    weights back to the best held-out epoch, so a gradient accumulated
        //    during the last epoch would describe parameters that were thrown
        //    away. See `Mlp::input_grad_importance`.
        //  * It reads THE ROWS THE OPTIMIZER SAW. Attributing over the
        //    held-out/val rows would make the sidecar's importance a function of
        //    rows this model is about to SCORE — a reporting-side leak of exactly
        //    the kind the rest of this module is built to avoid. (The inner
        //    early-stopping carve is train-fold data and would be defensible too,
        //    but it is excluded: `fit_rows` is precisely the index set that
        //    produced these weights.)
        //  * It touches no weight, no RNG draw and no optimizer step, so the fit
        //    is bit-identical whether or not this runs.
        //
        // Gated so it costs nothing when nobody asked: roughly one epoch of
        // forward/backward on top of the ~20-30 a fold runs. Same
        // `tracing::enabled!` shape as the held-out trace above — a diagnostic
        // nobody asked to see should not cost a pass.
        let traced = tracing::enabled!(tracing::Level::DEBUG);
        let want_grad = cfg.importance == MlpImportance::Grad || traced;
        let importance_started = (profile_enabled && want_grad).then(Instant::now);
        let mut importance_loader = BatchVisitTimings::default();
        let grad_importance = want_grad.then(|| {
            let mut total = vec![0.0f64; width];
            let mut nrows = 0usize;
            let order: Vec<usize> = (0..batch_count(fit_rows.len(), cfg.batch_size)).collect();
            let mut consume = |batch: &MlpBatch| {
                let values = net.input_grad_importance(
                    &batch.x,
                    &batch.y,
                    &batch.w,
                    &cfg.loss,
                    batch.rows(),
                );
                for (dst, value) in total.iter_mut().zip(values) {
                    *dst += value as f64 * batch.rows() as f64;
                }
                nrows += batch.rows();
            };
            importance_loader += visit_prefetched_batches(
                data,
                &fit_rows,
                &transform,
                cfg.batch_size,
                &order,
                &mut consume,
                profile_enabled,
            );
            total
                .into_iter()
                .map(|v| {
                    if nrows > 0 {
                        (v / nrows as f64) as f32
                    } else {
                        0.0
                    }
                })
                .collect::<Vec<f32>>()
        });
        let importance_wall = importance_started.map_or(Duration::ZERO, |t| t.elapsed());

        // The per-fold COMPARISON line: both rankings, one fold, no second run
        // needed to see whether they disagree. At `debug` because it is 2 lines
        // per fold of ~10 columns each, which is trace-volume next to the `info`
        // summary below. The authoritative diff is still two runs (one per
        // `TIMSSEEK_MLP_IMPORTANCE` value), whose sidecars carry every column.
        if traced {
            let w1 = fold_inputs_to_lanes(&transform, &net.input_importance());
            tracing::debug!(
                "{tag}importance top {} by w1:   {}",
                IMPORTANCE_LOG_TOP_N,
                top_lanes(&w1, &names, IMPORTANCE_LOG_TOP_N),
            );
            if let Some(g) = grad_importance.as_ref() {
                let grad = fold_inputs_to_lanes(&transform, g);
                tracing::debug!(
                    "{tag}importance top {} by grad: {}",
                    IMPORTANCE_LOG_TOP_N,
                    top_lanes(&grad, &names, IMPORTANCE_LOG_TOP_N),
                );
            }
        }

        // ONE grep-able line per fold, at `info` — the summary a sweep reads.
        //
        // `info` and not `debug` because the per-epoch trace is already at
        // `debug`: turning that on to find the stopping epoch means parsing
        // `epochs` lines per fold to recover 5 numbers, which is exactly the
        // thing this line exists to avoid. The cost is bounded and tiny —
        // `N_RESCORE_FOLDS` lines per run, 3 today — and the CLI's default
        // filter is `info`, so a sweep's run log carries them with no extra
        // flags. Extract with:
        //   grep -o 'MLP fold summary:.*' run.log
        //
        // 1-based epochs, matching the per-epoch trace. `kept_epoch` /
        // `best_held_out_loss` are `none` when there was no held-out set (early
        // stopping off and the trace off) or when every measurement was NaN.
        let kept_epoch = match outcome.best_epoch {
            Some(b) => (b + 1).to_string(),
            None => "none".to_string(),
        };
        let best_held_out = match outcome.best_val_loss {
            Some(v) => format!("{v:.6}"),
            None => "none".to_string(),
        };
        tracing::info!(
            "MLP fold summary: fold={} epochs_run={} epoch_budget={} kept_epoch={} \
             best_held_out_loss={} final_train_loss={:.6} restored={} train_rows={} \
             fit_rows={} inputs={} lane_columns={} importance={:?}",
            fold,
            outcome.epochs_run,
            cfg.epochs,
            kept_epoch,
            best_held_out,
            outcome.final_train_loss,
            outcome.restored,
            train.len(),
            fit_rows.len(),
            width,
            ncols,
            cfg.importance,
        );

        if let Some(fit_started) = fit_started {
            let total = fit_started.elapsed();
            let accounted = transform_wall
                + runtime_profile.train_wall
                + runtime_profile.validation_wall
                + importance_wall;
            let other = total.saturating_sub(accounted);
            tracing::info!(
                "MLP runtime profile: fold={} total_s={:.3} transform_s={:.3} train_wall_s={:.3} \
                 validation_wall_s={:.3} importance_wall_s={:.3} other_s={:.3} \
                 train_forward_s={:.3} train_loss_s={:.3} train_backward_s={:.3} adam_s={:.3} \
                 validation_forward_s={:.3} validation_loss_s={:.3} \
                 train_prepare_cpu_s={:.3} train_wait_for_batch_s={:.3} \
                 validation_prepare_cpu_s={:.3} validation_wait_for_batch_s={:.3} \
                 importance_prepare_cpu_s={:.3} importance_wait_for_batch_s={:.3} \
                 train_batches={} validation_batches={} importance_batches={}",
                fold,
                total.as_secs_f64(),
                transform_wall.as_secs_f64(),
                runtime_profile.train_wall.as_secs_f64(),
                runtime_profile.validation_wall.as_secs_f64(),
                importance_wall.as_secs_f64(),
                other.as_secs_f64(),
                runtime_profile.train_forward.as_secs_f64(),
                runtime_profile.train_loss.as_secs_f64(),
                runtime_profile.train_backward.as_secs_f64(),
                runtime_profile.adam.as_secs_f64(),
                runtime_profile.validation_forward.as_secs_f64(),
                runtime_profile.validation_loss.as_secs_f64(),
                train_loader.prepare_cpu.as_secs_f64(),
                train_loader.consumer_wait.as_secs_f64(),
                validation_loader.prepare_cpu.as_secs_f64(),
                validation_loader.consumer_wait.as_secs_f64(),
                importance_loader.prepare_cpu.as_secs_f64(),
                importance_loader.consumer_wait.as_secs_f64(),
                runtime_profile.train_batches,
                runtime_profile.validation_batches,
                importance_loader.batches,
            );
        }

        Ok(MlpFoldModel {
            net: RefCell::new(net),
            transform,
            final_loss: outcome.final_train_loss,
            n_train_rows: fit_rows.len(),
            outcome,
            importance_kind: cfg.importance,
            grad_importance,
            inference_batch_size: cfg.batch_size.max(1),
        })
    }

    /// Score `rows` through the SAME transform that was fitted at `fit` time.
    ///
    /// Nothing is refitted here — not the cull set, not the moments, not the
    /// imputation means. A refit would make a row's score depend on the other
    /// rows it happened to be batched with, which is both a leak (the batch may
    /// contain rows the caller is holding out) and non-reproducible.
    ///
    /// The returned score is the raw logit, not a probability: `assign_qval`
    /// needs a monotone ranking and the sigmoid is a monotone no-op on it.
    ///
    /// Every returned value is FINITE, or this is an
    /// [`MlpFoldError::NonFiniteScore`] — see that variant for why the check lives
    /// here and not downstream.
    fn predict<D: FoldDataset>(&self, data: &D, rows: &[usize]) -> Result<Vec<f64>, MlpFoldError> {
        let names = data.column_names();
        let ncols = names.len();
        if ncols != self.transform.ncols_lane() {
            return Err(MlpFoldError::WidthMismatch {
                fitted: self.transform.ncols_lane(),
                got: ncols,
            });
        }
        if rows.is_empty() {
            return Ok(Vec::new());
        }

        // A column the TRAIN rows never saw a non-finite value in has no
        // `_isna` companion, so a NaN arriving here is imputed to the train
        // mean and becomes indistinguishable from an average row. That is a
        // train/score distribution mismatch worth knowing about, and it is
        // reported at runtime rather than asserted: `debug_assert!` compiles
        // out in release, which is exactly where a production run would hit it.
        let mut dirty_flags = vec![false; ncols];
        let mut row = vec![0.0f64; ncols];
        let mut scores = Vec::with_capacity(rows.len());
        let mut net = self.net.borrow_mut();
        let mut x = Tensor::default();
        for chunk in rows.chunks(self.inference_batch_size) {
            x.reshape(chunk.len(), self.transform.width());
            for (k, &i) in chunk.iter().enumerate() {
                data.get_values(i, &mut row);
                self.transform.mark_dirty_clean(&row, &mut dirty_flags);
                self.transform.apply(&row, x.row_mut(k));
            }
            let out = net.forward(&x, false);
            scores.extend((0..chunk.len()).map(|i| out.row(i)[0] as f64));
        }
        let dirty: Vec<usize> = dirty_flags
            .iter()
            .enumerate()
            .filter_map(|(j, &hit)| hit.then_some(j))
            .collect();
        if !dirty.is_empty() {
            let dirty_names: Vec<&str> = dirty.iter().map(|&j| &*names[j]).collect();
            tracing::warn!(
                "MLP predict: {} scored rows carried non-finite values in {} column(s) that were \
                 clean on the train rows, so they impute to the train mean with no _isna flag: \
                 {:?}",
                rows.len(),
                dirty.len(),
                dirty_names,
            );
        }

        // A hard error, not a `debug_assert!`: release is where a diverged net
        // would actually be met, and the alternative is aborting in `assign_qval`
        // with a message about the SORT. See `MlpFoldError::NonFiniteScore`.
        let n_bad = scores.iter().filter(|s| !s.is_finite()).count();
        if n_bad > 0 {
            let first = scores.iter().position(|s| !s.is_finite()).unwrap();
            return Err(MlpFoldError::NonFiniteScore {
                row: rows[first],
                n_bad,
                scored: rows.len(),
            });
        }

        Ok(scores)
    }

    /// Whichever of the TWO metrics [`MlpConfig::importance`] selected, mapped
    /// back to lane indices by [`fold_inputs_to_lanes`].
    ///
    ///  * [`MlpImportance::W1`] (THE DEFAULT, and the content the externally
    ///    consumed sidecar has always carried): `sum_o |W1[j, o]|` over the first
    ///    layer.
    ///  * [`MlpImportance::Grad`]: mean `|dL/dx_j|` over the train rows, measured
    ///    at `fit` time by [`Mlp::input_grad_importance`].
    ///
    /// See [`MlpImportance`] for what each one can and cannot see, and
    /// [`MlpConfig::from_env`] for the `TIMSSEEK_MLP_IMPORTANCE` selector.
    ///
    /// NEVER `NAN`, i.e. never "unreported" under the [`FoldModel::importance`]
    /// contract, for EITHER metric: the first layer has a weight for every input
    /// and the backward pass writes a gradient for every input, so every lane
    /// column that survived the cull has a measurement, and the culled ones report
    /// a finite `0.0` — see [`fold_inputs_to_lanes`], which is where both of those
    /// properties live for both metrics.
    fn importance(&self) -> Vec<f32> {
        let raw = match self.importance_kind {
            MlpImportance::W1 => self.net.borrow().input_importance(),
            // `fit` computes this whenever `importance_kind` is `Grad`, so the
            // fallback is unreachable. It is a clone-or-empty rather than an
            // `expect` because an importance report is a diagnostic: aborting a
            // finished search over a sidecar column would be a worse failure than
            // an empty one, and the `debug_assert` names the bug in test builds.
            MlpImportance::Grad => {
                debug_assert!(
                    self.grad_importance.is_some(),
                    "fit must measure the gradient attribution when it is the selected metric"
                );
                self.grad_importance.clone().unwrap_or_default()
            }
        };
        fold_inputs_to_lanes(&self.transform, &raw)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::ml::cv::{
        PrecomputedFeatures,
        RowMajorDataset,
    };
    use rand::SeedableRng;
    use rand::distr::{
        Distribution,
        Uniform,
    };

    /// Test config: same architecture family as the default, shrunk so the
    /// seed sweeps stay fast.
    fn cfg(seed: u64) -> MlpConfig {
        MlpConfig {
            hidden: vec![16, 8],
            lr: 1e-3,
            epochs: 150,
            batch_size: 64,
            seed,
            ..MlpConfig::default()
        }
    }

    fn dataset(feat: Vec<f64>, ncols: usize, responses: Vec<f64>) -> RowMajorDataset {
        let names: Vec<Arc<str>> = (0..ncols)
            .map(|j| Arc::from(format!("feature_{j}").as_str()))
            .collect();
        RowMajorDataset::new(
            PrecomputedFeatures::from_row_major(feat, ncols, responses),
            names,
            2,
        )
    }

    #[test]
    fn prefetcher_retains_only_two_buffers_and_preserves_chunk_order() {
        let nrows = 10usize;
        let feat: Vec<f64> = (0..nrows)
            .flat_map(|i| [i as f64, (i * i) as f64])
            .collect();
        let data = dataset(feat, 2, vec![0.0; nrows]);
        let rows: Vec<usize> = (0..nrows).collect();
        let transform = ColumnTransform::fit_streaming(2, rows.iter().copied(), |i, out| {
            data.get_values(i, out)
        });
        let order = [2usize, 0, 1];
        let mut seen_rows = Vec::new();
        let mut buffers = std::collections::HashSet::new();
        visit_prefetched_batches(
            &data,
            &rows,
            &transform,
            4,
            &order,
            &mut |batch| {
                buffers.insert(batch.x.row(0).as_ptr() as usize);
                for i in 0..batch.rows() {
                    seen_rows.extend_from_slice(batch.x.row(i));
                }
            },
            false,
        );

        let source_order = [8usize, 9, 0, 1, 2, 3, 4, 5, 6, 7];
        let mut expected = Vec::new();
        let mut raw = vec![0.0; 2];
        let mut transformed = vec![0.0; transform.width()];
        for i in source_order {
            data.get_values(i, &mut raw);
            transform.apply(&raw, &mut transformed);
            expected.extend_from_slice(&transformed);
        }
        assert_eq!(seen_rows, expected);
        assert!(
            buffers.len() <= 2,
            "the prefetcher retained more than its two reusable buffers: {buffers:?}"
        );
    }

    /// `n_each` targets then `n_each` decoys, 3 columns. Columns 0 and 1 are
    /// drawn from twice the decoy range for targets (separable); column 2 is
    /// pure noise in both classes. `shift` moves every value, so two calls with
    /// the same seed and different shifts differ only by a translation.
    fn toy(n_each: usize, seed: u64, shift: f64) -> (Vec<f64>, Vec<f64>) {
        let noise = Uniform::try_from(1.0..10.0).unwrap();
        let signal = Uniform::try_from(1.0..20.0).unwrap();
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
        let mut feat = Vec::with_capacity(2 * n_each * 3);
        let mut y = Vec::with_capacity(2 * n_each);
        for _ in 0..n_each {
            feat.push(signal.sample(&mut rng) + shift);
            feat.push(signal.sample(&mut rng) + shift);
            feat.push(noise.sample(&mut rng) + shift);
            y.push(1.0);
        }
        for _ in 0..n_each {
            feat.push(noise.sample(&mut rng) + shift);
            feat.push(noise.sample(&mut rng) + shift);
            feat.push(noise.sample(&mut rng) + shift);
            y.push(0.0);
        }
        (feat, y)
    }

    /// Fraction of (target, decoy) pairs the scores order correctly — an AUC,
    /// which needs no threshold and is invariant to the logit's offset.
    fn pair_auc(scores: &[f64], is_target: &[bool]) -> f64 {
        let t: Vec<f64> = scores
            .iter()
            .zip(is_target)
            .filter(|&(_, &b)| b)
            .map(|(&s, _)| s)
            .collect();
        let d: Vec<f64> = scores
            .iter()
            .zip(is_target)
            .filter(|&(_, &b)| !b)
            .map(|(&s, _)| s)
            .collect();
        let mut ok = 0usize;
        for a in &t {
            for b in &d {
                if a > b {
                    ok += 1;
                }
            }
        }
        ok as f64 / (t.len() * d.len()) as f64
    }

    /// Rows are laid out targets-then-decoys, so a train/held split has to
    /// interleave rather than cut: take every other row.
    fn split(nrows: usize) -> (Vec<usize>, Vec<usize>) {
        (
            (0..nrows).step_by(2).collect(),
            (1..nrows).step_by(2).collect(),
        )
    }

    /// The learning test. Several seeds on purpose: building `mlp.rs` turned up
    /// two distinct init-dependent training traps (dead ReLUs, and full-batch
    /// descent settling into one), each of which a single-seed test would have
    /// scored as a clean pass or a backprop bug purely on luck.
    ///
    /// The bar is 0.85 rather than something near 1.0 because this toy is not
    /// perfectly separable and cannot be: a target's informative columns are
    /// drawn from twice the decoy range, so `(9/19)^2 ~ 22%` of targets land
    /// entirely inside the decoy support and are indistinguishable from one.
    /// The Bayes-optimal AUC is therefore only ~0.89, and the fits below come
    /// in at 0.87-0.89.
    ///
    /// BOTH bars are asserted and the ABSOLUTE one is the tighter of the two on
    /// this fixture: the best single raw column comes in at ~0.763, so the
    /// relative bar sits at ~0.813 and the absolute 0.85 subsumes it. The
    /// relative bar is kept anyway, and is the one to keep if the fixture
    /// changes — it says the model COMBINED the two informative features rather
    /// than latching onto one, which is a claim about the fit; 0.85 is a number
    /// that happens to be right for these draws.
    #[test]
    fn fitted_model_separates_a_separable_toy_on_held_out_rows() {
        for seed in [7u64, 13, 42, 1234] {
            let (feat, y) = toy(200, seed, 0.0);
            let data = dataset(feat.clone(), 3, y.clone());
            let (train, held) = split(400);

            let model = MlpFoldModel::fit(&cfg(seed), &data, 0, &train, &[]).unwrap();
            let scores = model.predict(&data, &held).unwrap();
            let is_target: Vec<bool> = held.iter().map(|&i| y[i] > 0.5).collect();

            let auc = pair_auc(&scores, &is_target);
            let best_raw = (0..3)
                .map(|j| {
                    let raw: Vec<f64> = held.iter().map(|&i| feat[i * 3 + j]).collect();
                    pair_auc(&raw, &is_target)
                })
                .fold(0.0f64, f64::max);
            assert!(
                auc > 0.85,
                "seed {seed}: held-out AUC {auc} (final loss {})",
                model.final_loss()
            );
            assert!(
                auc > best_raw + 0.05,
                "seed {seed}: held-out AUC {auc} barely beats the best raw column \
                 ({best_raw}) — the fit added nothing"
            );
        }
    }

    /// THE property this module exists to keep, half one: the transform is
    /// fitted from the TRAIN rows, so rows that are neither trained on nor
    /// scored cannot move a score.
    ///
    /// Two datasets share their train rows AND the held rows being scored, and
    /// differ only in a block of extra rows with a wildly different mean and
    /// scale. If the fit reached past `train` — over all rows, say — those extra
    /// rows would move the standardization and the scores would drift.
    /// Bit-identical is the assertion, not approximately-equal: nothing about a
    /// seeded fit is allowed to be fuzzy.
    #[test]
    fn fit_ignores_rows_outside_the_train_slice() {
        for seed in [7u64, 13, 42, 1234] {
            let (base, y) = toy(150, seed, 0.0);
            let (train, held) = split(300);

            // Same 300 rows in both, then 100 extra rows: mild in A, extreme
            // in B. Extra rows are labelled so the matrix stays well-formed;
            // no fit or predict call below ever names their indices.
            let build = |extra_shift: f64, extra_scale: f64| {
                let (extra, extra_y) = toy(50, seed ^ 0xABCD, 0.0);
                let mut feat = base.clone();
                feat.extend(extra.iter().map(|v| v * extra_scale + extra_shift));
                let mut resp = y.clone();
                resp.extend(extra_y);
                dataset(feat, 3, resp)
            };
            let a = build(0.0, 1.0);
            let b = build(1.0e6, 1.0e3);

            let sa = MlpFoldModel::fit(&cfg(seed), &a, 0, &train, &[])
                .unwrap()
                .predict(&a, &held)
                .unwrap();
            let sb = MlpFoldModel::fit(&cfg(seed), &b, 0, &train, &[])
                .unwrap()
                .predict(&b, &held)
                .unwrap();
            assert_eq!(
                sa, sb,
                "seed {seed}: rows outside the train slice moved the scores"
            );

            // Non-vacuity: the extra rows are not inert data the fit would
            // ignore anyway — training ON them changes the answer.
            let with_extra: Vec<usize> = train.iter().copied().chain(300..400).collect();
            let sc = MlpFoldModel::fit(&cfg(seed), &b, 0, &with_extra, &[])
                .unwrap()
                .predict(&b, &held)
                .unwrap();
            assert!(
                sc != sa,
                "seed {seed}: the extra rows changed nothing even when trained on, \
                 so the assertion above proves nothing"
            );
        }
    }

    /// Half two: `predict` applies the STORED transform, so a row's score does
    /// not depend on which other rows it was batched with. A refit inside
    /// `predict` — from the rows being scored — would make it depend on exactly
    /// that, and would also be a leak whenever the batch spans held-out rows.
    ///
    /// The held rows here are shifted far away from the train rows, so a refit
    /// would land on visibly different moments rather than coincidentally
    /// similar ones.
    #[test]
    fn predict_scores_a_row_the_same_alone_as_in_a_batch() {
        for seed in [7u64, 13, 42, 1234] {
            let (mut feat, y) = toy(150, seed, 0.0);
            // Push the HELD rows far off, so the scored batch has a different
            // mean and spread from the train rows. `split` takes even rows as
            // train and odd rows as held, and the layout is targets-then-decoys,
            // so "every odd row" is half of each class — the held half.
            for i in (0..300).filter(|i| i % 2 == 1) {
                for j in 0..3 {
                    feat[i * 3 + j] = feat[i * 3 + j] * 50.0 + 500.0;
                }
            }
            let data = dataset(feat, 3, y);
            let (train, held) = split(300);

            let model = MlpFoldModel::fit(&cfg(seed), &data, 0, &train, &[]).unwrap();
            let batched = model.predict(&data, &held).unwrap();
            for (k, &row) in held.iter().enumerate() {
                let alone = model.predict(&data, &[row]).unwrap();
                assert_eq!(
                    alone[0], batched[k],
                    "seed {seed}: row {row} scored differently alone than in a batch, \
                     so `predict` is refitting from the rows it is given"
                );
            }
        }
    }

    /// The MLP is the only rescorer whose discriminant is not finite by
    /// construction, and a non-finite one must surface HERE rather than as
    /// `assign_qval`'s "Expecting scores to be sorted in descending order" several
    /// stages downstream. See [`MlpFoldError::NonFiniteScore`].
    ///
    /// Provoked by poisoning the fitted net's output weights, because no INPUT can
    /// reach this: the transform imputes every non-finite feature value before the
    /// forward pass, so a diverged fit is the only route and this test manufactures
    /// its end state directly rather than trying to find a config that diverges.
    ///
    /// The clean predict on the same model is the non-vacuity half — it says the
    /// error came from the poison and not from the fixture.
    #[test]
    fn predict_rejects_a_non_finite_logit_instead_of_returning_it() {
        let (feat, y) = toy(60, 7, 0.0);
        let data = dataset(feat, 3, y);
        let (train, held) = split(120);

        let model = MlpFoldModel::fit(&cfg(7), &data, 0, &train, &[]).unwrap();
        assert!(
            model
                .predict(&data, &held)
                .unwrap()
                .iter()
                .all(|s| s.is_finite()),
            "the fixture must score finite before the poison, or the assertion below \
             would pass for the wrong reason"
        );

        // NaN into every parameter of the net. The last `Linear` writes the logit,
        // so this is guaranteed to reach the output regardless of layer layout.
        for (p, _grad, _decay) in model.net.borrow_mut().params_and_grads() {
            p.fill(f32::NAN);
        }

        let err = model
            .predict(&data, &held)
            .expect_err("a NaN logit must not be returned as a score");
        assert_eq!(
            err,
            MlpFoldError::NonFiniteScore {
                row: held[0],
                n_bad: held.len(),
                scored: held.len(),
            },
            "the error has to name the first offending DATASET row (not its index in \
             the batch) and how many of the batch went bad"
        );
        let msg = err.to_string();
        assert!(
            msg.contains("NON-FINITE") && msg.contains("diverged"),
            "the message must say what happened and where to look: {msg}"
        );
    }

    /// The [`FoldModel::importance`] contract: culled columns report a FINITE
    /// `0.0` (which reaches the sidecar and tells the operator they were
    /// culled), never the `NAN` "unreported" sentinel (which would be dropped).
    #[test]
    fn importance_reports_zero_for_culled_columns() {
        for seed in [7u64, 13, 42, 1234] {
            let (base, y) = toy(150, seed, 0.0);
            // Widen to 5 columns: 0,1 informative, 2 noise, 3 all-NaN, 4 constant.
            let mut feat = Vec::with_capacity(300 * 5);
            for i in 0..300 {
                feat.extend_from_slice(&base[i * 3..i * 3 + 3]);
                feat.push(f64::NAN);
                feat.push(7.0);
            }
            let data = dataset(feat, 5, y);
            let (train, _) = split(300);

            let model = MlpFoldModel::fit(&cfg(seed), &data, 0, &train, &[]).unwrap();
            assert_eq!(model.transform().culled(), &[3, 4]);

            let imp = model.importance();
            assert_eq!(imp.len(), 5, "importance must be lane-indexed, full width");
            assert!(
                imp.iter().all(|v| v.is_finite()),
                "seed {seed}: the NAN sentinel would delete the cull report: {imp:?}"
            );
            assert_eq!(imp[3], 0.0, "an all-NaN culled column reports 0.0");
            assert_eq!(imp[4], 0.0, "a constant culled column reports 0.0");
            assert!(
                imp[0] > 0.0 && imp[1] > 0.0 && imp[2] > 0.0,
                "seed {seed}: surviving columns must carry weight, else the zeros \
                 above are vacuous: {imp:?}"
            );
        }
    }

    /// A missable column's `_isna` companion has no lane of its own, so its
    /// weight is folded into the column it flags rather than dropped. The
    /// no-weight-lost half is a `debug_assert` inside `importance`, which this
    /// test (a debug build) executes; asserted here is that the companion exists
    /// at all, so that assert is not vacuous, AND that it maps back to the column
    /// it flags — without which this test would pass for an
    /// `isna_lane_of_input` that always answered lane 0, i.e. for an
    /// implementation that folded every companion's weight into the wrong
    /// feature.
    #[test]
    fn missable_columns_fold_their_companion_into_the_parent_lane() {
        for seed in [7u64, 13, 42] {
            let (mut feat, y) = toy(150, seed, 0.0);
            // Punch NaNs into column 2 on every 7th row, train and held alike.
            for i in (0..300).step_by(7) {
                feat[i * 3 + 2] = f64::NAN;
            }
            let data = dataset(feat, 3, y);
            let (train, held) = split(300);

            let model = MlpFoldModel::fit(&cfg(seed), &data, 0, &train, &[]).unwrap();
            assert!(model.transform().culled().is_empty());
            assert_eq!(
                model.transform().width(),
                4,
                "three survivors + one _isna companion"
            );
            // The companion is transformed input 3 (it follows the three
            // standardized columns) and it flags LANE 2 — the column the NaNs
            // were punched into. This is the property in this test's name.
            assert_eq!(model.transform().isna_lane_of_input(3), Some(2));
            assert_eq!(
                model.transform().isna_lane_of_input(2),
                None,
                "inputs 0..3 are standardized columns, not companions"
            );

            let imp = model.importance();
            assert_eq!(imp.len(), 3);
            assert!(imp.iter().all(|v| v.is_finite() && *v > 0.0), "{imp:?}");

            // Scoring rows that carry the NaN must still work end to end.
            let scores = model.predict(&data, &held).unwrap();
            assert!(scores.iter().all(|s| s.is_finite()), "seed {seed}");
        }
    }

    // ------------------------------------------- input-gradient importance

    /// [`cfg`] with the gradient metric selected. Nothing else differs, which is
    /// the premise of every test below.
    fn grad_cfg(seed: u64) -> MlpConfig {
        MlpConfig {
            importance: MlpImportance::Grad,
            ..cfg(seed)
        }
    }

    /// The 5-column fixture the two importance tests share: lanes 0 and 1
    /// informative, 2 pure noise, 3 all-NaN (culled), 4 constant (culled).
    fn importance_fixture(seed: u64) -> (RowMajorDataset, Vec<usize>) {
        let (base, y) = toy(150, seed, 0.0);
        let mut feat = Vec::with_capacity(300 * 5);
        for i in 0..300 {
            feat.extend_from_slice(&base[i * 3..i * 3 + 3]);
            feat.push(f64::NAN);
            feat.push(7.0);
        }
        (dataset(feat, 5, y), split(300).0)
    }

    /// THE PROPERTY THAT MAKES THE TWO METRICS COMPARABLE: selecting one changes
    /// the sidecar and NOTHING ELSE. Two runs of the same search, one with
    /// `TIMSSEEK_MLP_IMPORTANCE=grad`, must produce the same scores down to the
    /// bit, or the diff of the two rankings is a diff of two different fits.
    ///
    /// The gradient metric runs a real forward/BACKWARD pass over the train rows
    /// after training, so this is not a formality: it is the observable
    /// consequence of that pass taking no optimizer step, drawing no RNG and
    /// leaving no stale gradient (see
    /// `mlp::test::input_grad_importance_leaves_the_fit_untouched`, which asserts
    /// the mechanism; this asserts the effect at the fold boundary).
    ///
    /// The second half is non-vacuity: the two metrics must actually DISAGREE
    /// somewhere, or "both are reportable" would be an expensive way to report one
    /// number twice.
    #[test]
    fn the_importance_metric_choice_does_not_move_the_scores() {
        let mut ranking_differs = false;
        for seed in [7u64, 13, 42, 1234] {
            let (data, train) = importance_fixture(seed);
            let held: Vec<usize> = (1..300).step_by(2).collect();

            let w1 = MlpFoldModel::fit(&cfg(seed), &data, 0, &train, &[]).unwrap();
            let grad = MlpFoldModel::fit(&grad_cfg(seed), &data, 0, &train, &[]).unwrap();

            assert_eq!(
                w1.predict(&data, &held).unwrap(),
                grad.predict(&data, &held).unwrap(),
                "seed {seed}: selecting an importance metric moved the SCORES"
            );
            assert_eq!(
                w1.outcome(),
                grad.outcome(),
                "seed {seed}: the stopping decision moved"
            );

            let (iw, ig) = (w1.importance(), grad.importance());
            assert_ne!(
                iw, ig,
                "seed {seed}: the two metrics reported the same vector"
            );

            // The `debug` comparison line builds its `w1` row from the FITTED NET,
            // because a `Grad`-configured model's `importance()` returns the other
            // metric. Same numbers, which is what makes the two log rows comparable
            // to the two runs' sidecars. This is also the only coverage of that
            // branch's body: the branch itself is gated on
            // `tracing::enabled!(DEBUG)` and no test here installs a subscriber.
            assert_eq!(
                fold_inputs_to_lanes(grad.transform(), &grad.net.borrow().input_importance()),
                iw,
                "seed {seed}: the log's w1 row must be the number the sidecar would carry"
            );

            let top = |v: &[f32]| {
                (0..v.len())
                    .max_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap())
                    .unwrap()
            };
            ranking_differs |= top(&iw) != top(&ig);
        }
        // MEASURED, and the entire justification for a second metric: on this
        // fixture `W1` ranks the PURE-NOISE lane (2) top on seeds 42 and 1234,
        // where the gradient ranks it last on all four. If this ever fails, the
        // two metrics have started agreeing on this fixture and the fixture — not
        // the assertion — is what needs replacing.
        assert!(
            ranking_differs,
            "the two metrics agreed on the top lane for every seed, so nothing here \
             demonstrates that comparing them is worth the pass"
        );
    }

    /// The [`FoldModel::importance`] contract, for the SECOND metric: culled
    /// columns report a FINITE `0.0` (which reaches the sidecar and tells the
    /// operator they were culled), never the `NAN` "unreported" sentinel (which
    /// `fold_feature_stats` filters out).
    ///
    /// The mirror of [`importance_reports_zero_for_culled_columns`], asserted
    /// separately rather than parameterized over the metric: the two arrive at the
    /// zero by different routes — `W1` has no weight row for a culled column,
    /// `Grad` has no gradient for one — and a regression in either could leave the
    /// other passing.
    #[test]
    fn grad_importance_reports_zero_for_culled_columns() {
        for seed in [7u64, 13, 42, 1234] {
            let (data, train) = importance_fixture(seed);
            let model = MlpFoldModel::fit(&grad_cfg(seed), &data, 0, &train, &[]).unwrap();
            assert_eq!(model.transform().culled(), &[3, 4]);

            let imp = model.importance();
            assert_eq!(imp.len(), 5, "importance must be lane-indexed, full width");
            assert!(
                imp.iter().all(|v| v.is_finite()),
                "seed {seed}: the NAN sentinel would delete the cull report: {imp:?}"
            );
            assert_eq!(imp[3], 0.0, "an all-NaN culled column reports 0.0");
            assert_eq!(imp[4], 0.0, "a constant culled column reports 0.0");
            assert!(
                imp[0] > 0.0 && imp[1] > 0.0 && imp[2] > 0.0,
                "seed {seed}: surviving columns must carry attribution, else the zeros \
                 above are vacuous: {imp:?}"
            );

            // The reported measurement survives to the sidecar shape: `0.0` is
            // present, and it is the culled lanes that carry it.
            let stats = crate::ml::cv::fold_feature_stats(
                &data,
                std::slice::from_ref(&train),
                &[Some(&model)],
            );
            let imp_rows = &stats[0].feature_importance;
            assert_eq!(
                imp_rows.len(),
                5,
                "seed {seed}: every reported column must reach the sidecar: {imp_rows:?}"
            );
            let zeros: Vec<&str> = imp_rows
                .iter()
                .filter(|(_, v)| *v == 0.0)
                .map(|(n, _)| &**n)
                .collect();
            assert_eq!(
                zeros,
                vec!["feature_3", "feature_4"],
                "seed {seed}: the culled lanes are the ones that must show up as 0.0"
            );
        }
    }

    /// What the gradient metric buys over the weight norm, at the fold boundary:
    /// it ranks the pure-noise lane BELOW both informative ones.
    ///
    /// Seed-swept, because it is a claim about a trained net.
    ///
    /// MEASURED on this fixture (lane 2 is the noise column):
    ///
    /// | seed | grad[0] | grad[1] | grad[2] | w1[0] | w1[1] | w1[2] |
    /// |------|---------|---------|---------|-------|-------|-------|
    /// | 7    | 0.228   | 0.252   | 0.102   | 11.66 | 15.30 | 13.17 |
    /// | 13   | 0.184   | 0.220   | 0.049   | 12.79 | 14.00 |  9.44 |
    /// | 42   | 0.241   | 0.228   | 0.068   | 10.38 | 13.60 | 14.57 |
    /// | 1234 | 0.170   | 0.233   | 0.088   |  9.80 |  7.70 | 11.76 |
    ///
    /// `W1` ranks the noise lane FIRST on seeds 42 and 1234 and second on seed 7.
    /// Chunk-order training changes the exact trajectory, but the gradient still
    /// ranks it last on all four seeds.
    #[test]
    fn grad_importance_ranks_the_noise_lane_below_the_informative_ones() {
        for seed in [7u64, 13, 42, 1234] {
            let (data, train) = importance_fixture(seed);
            let model = MlpFoldModel::fit(&grad_cfg(seed), &data, 0, &train, &[]).unwrap();
            let g = model.importance();
            assert!(
                g[2] < g[0].min(g[1]),
                "seed {seed}: the noise lane must rank below both informative lanes: {g:?}"
            );
        }
    }

    /// Same config + same fold + same data => bit-identical gradient attribution.
    ///
    /// The sidecar exists to be diffed across runs, so this is the property that
    /// makes the diff readable at all. It also covers the `_isna` companion route:
    /// lane 2 carries NaNs here, so its reported value is the sum of its own
    /// gradient and its companion's, and the `debug_assert` inside
    /// [`fold_inputs_to_lanes`] executes on this (debug) build.
    #[test]
    fn grad_importance_is_bit_reproducible() {
        for seed in [7u64, 13, 42] {
            let (mut feat, y) = toy(150, seed, 0.0);
            // NaNs into lane 2 on every 7th row, so it becomes `Missable` and gets
            // an `_isna` companion whose attribution folds back into it.
            for i in (0..300).step_by(7) {
                feat[i * 3 + 2] = f64::NAN;
            }
            let data = dataset(feat, 3, y);
            let (train, _) = split(300);

            let a = MlpFoldModel::fit(&grad_cfg(seed), &data, 0, &train, &[]).unwrap();
            let b = MlpFoldModel::fit(&grad_cfg(seed), &data, 0, &train, &[]).unwrap();
            assert_eq!(
                a.transform().width(),
                4,
                "seed {seed}: three survivors + one _isna companion, or the companion \
                 folding below is not exercised"
            );
            assert_eq!(a.transform().isna_lane_of_input(3), Some(2));

            let (ia, ib) = (a.importance(), b.importance());
            assert_eq!(
                ia.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                ib.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                "seed {seed}: the attribution moved between two identical fits"
            );
            assert_eq!(ia.len(), 3, "lane-indexed, companions folded in");
            assert!(
                ia.iter().all(|v| v.is_finite() && *v > 0.0),
                "seed {seed}: {ia:?}"
            );
        }
    }

    /// [`top_lanes`] is what the per-fold `debug` comparison line is made of, and
    /// it is not otherwise executed by this suite: the branch that calls it is
    /// gated on `tracing::enabled!(DEBUG)`, and no test here installs a subscriber,
    /// so with no direct test the ordering could be reversed and nothing would say
    /// so.
    ///
    /// Descending, truncated to `n`, and TIES KEEP LANE ORDER — the same rule
    /// `fold_feature_stats` sorts the sidecar by, so the log and the file agree on
    /// which of two equal columns comes first.
    #[test]
    fn the_comparison_log_names_the_strongest_lanes_first() {
        let names: Vec<Arc<str>> = ["a", "b", "c", "d", "e"]
            .iter()
            .map(|s| Arc::from(*s))
            .collect();
        let imp = [0.5f32, 0.0, 2.0, 0.0, 1.0];

        let line = top_lanes(&imp, &names, 3);
        assert_eq!(line, "c=2.000000 e=1.000000 a=0.500000");

        // Ties in lane order, and `n` past the end is not an error.
        let all = top_lanes(&imp, &names, 99);
        assert_eq!(
            all, "c=2.000000 e=1.000000 a=0.500000 b=0.000000 d=0.000000",
            "equal values must keep lane order, matching the sidecar's stable sort"
        );
    }

    /// A dataset of the wrong width would silently read column `j` of a
    /// different matrix, so it is an error rather than a best effort.
    #[test]
    fn predict_rejects_a_dataset_of_a_different_width() {
        let (feat, y) = toy(50, 7, 0.0);
        let narrow = dataset(feat.clone(), 3, y.clone());
        let (train, held) = split(100);
        let model = MlpFoldModel::fit(&cfg(7), &narrow, 0, &train, &[]).unwrap();

        // Same rows, one extra column.
        let mut wide_feat = Vec::with_capacity(100 * 4);
        for i in 0..100 {
            wide_feat.extend_from_slice(&feat[i * 3..i * 3 + 3]);
            wide_feat.push(i as f64);
        }
        let wide = dataset(wide_feat, 4, y);

        assert_eq!(
            model.predict(&wide, &held),
            Err(MlpFoldError::WidthMismatch { fitted: 3, got: 4 })
        );
        // The narrow dataset it WAS fitted on still scores, so the check is not
        // rejecting everything.
        assert!(model.predict(&narrow, &held).is_ok());
    }

    /// Same config + same fold + same data => bit-identical scores, run to run.
    /// The whole q-value pipeline downstream is a ranking of these.
    ///
    /// The STOPPING DECISION is part of what has to reproduce: an early stop is
    /// a branch on a measured float, and a run that stopped at a different epoch
    /// would produce a different (still finite, still plausible) score vector.
    /// `TrainOutcome` carries the epoch it ran to and the epoch it kept, and
    /// both are asserted equal here.
    #[test]
    fn refitting_the_same_fold_reproduces_the_scores_exactly() {
        // 200 train rows, i.e. over `MIN_TRAIN_ROWS_FOR_INNER_VAL`, so this
        // fixture exercises the carve and the stop rather than a fixed budget.
        let (feat, y) = toy(200, 42, 0.0);
        let data = dataset(feat, 3, y);
        let (train, held) = split(400);

        let a = MlpFoldModel::fit(&cfg(42), &data, 0, &train, &[]).unwrap();
        let b = MlpFoldModel::fit(&cfg(42), &data, 0, &train, &[]).unwrap();
        assert_eq!(
            a.predict(&data, &held).unwrap(),
            b.predict(&data, &held).unwrap()
        );
        assert_eq!(a.importance(), b.importance());
        assert_eq!(a.outcome(), b.outcome(), "the stopping decision moved");
        // Non-vacuity for the line above: this fixture goes down the
        // early-stopping path — an inner slice was carved and a best epoch was
        // chosen — so "same outcome" is a comparison of two real decisions and
        // not of two inert budgets. Whether rollback happens depends on whether
        // the final chunk-order epoch is the best one; the stop itself is covered
        // by `mlp::test::early_stopping_stops_early_and_does_not_hurt_the_held_out_loss`.
        assert!(
            a.outcome().best_epoch.is_some(),
            "this train slice ({} rows) is supposed to carve an inner validation \
             slice and choose a best epoch: {:?}",
            train.len(),
            a.outcome()
        );

        // ...and the `fold` ARGUMENT really does reach the init: same config,
        // same rows, different fold => different weights. Without this the
        // parameter could be dropped on the floor and every fold would train
        // from the same initialization, which is still leak-free and so would
        // not fail any other test here.
        let other = MlpFoldModel::fit(&cfg(42), &data, 1, &train, &[]).unwrap();
        assert!(other.predict(&data, &held).unwrap() != a.predict(&data, &held).unwrap());
    }

    /// An all-dead matrix has nothing to train on; that is an error rather than
    /// a net with zero inputs.
    ///
    /// The error's payload is asserted, not just its variant: it is the ONLY
    /// diagnostic the abort in `qvalues::rescore_mlp_with` can print, so a fit
    /// that reported the wrong fold, the wrong lane width or the wrong row count
    /// would send an operator looking in the wrong place.
    #[test]
    fn fit_rejects_a_fully_culled_matrix() {
        let feat = vec![1.0f64; 20 * 2];
        let y: Vec<f64> = (0..20).map(|i| (i % 2) as f64).collect();
        let data = dataset(feat, 2, y);
        let train: Vec<usize> = (0..20).collect();
        let err = MlpFoldModel::fit(&cfg(7), &data, 0, &train, &[]).err();
        assert_eq!(
            err,
            Some(MlpFoldError::NoUsableColumns {
                fold: 0,
                ncols: 2,
                train_rows: 20,
            })
        );
        let msg = err.unwrap().to_string();
        assert!(
            msg.contains("all 2 lane columns were culled") && msg.contains("20 train rows"),
            "the culled-lane message must name the lane width and the row count: {msg}"
        );
    }

    /// A ZERO-ROW fold reaches [`MlpFoldError::NoUsableColumns`] too — every
    /// column is vacuously dead when there is nothing to measure it on — and the
    /// message must say THAT rather than blaming the cull.
    ///
    /// This is the realistic way the variant fires: a genuinely constant column
    /// of any magnitude keeps a hair of floating-point variance in
    /// `sumsq/n - mean^2` at most row counts, so it survives `MIN_STD` and the
    /// "every column is dead" reading is the rarer one. The message used to
    /// assert it in both cases, which is a wrong-diagnosis bug, not a wording
    /// nit: it points the reader at the feature set when the fold partition or
    /// the input size is what broke.
    #[test]
    fn an_empty_train_slice_reports_zero_rows_rather_than_blaming_the_cull() {
        let (feat, y) = toy(20, 7, 2.0);
        let data = dataset(feat, 3, y);

        let err = MlpFoldModel::fit(&cfg(7), &data, 2, &[], &[]).err();
        assert_eq!(
            err,
            Some(MlpFoldError::NoUsableColumns {
                fold: 2,
                ncols: 3,
                train_rows: 0,
            })
        );
        let msg = err.unwrap().to_string();
        assert!(
            msg.contains("train slice is EMPTY") && msg.contains("fold 2"),
            "an empty train slice must be reported as such, not as a culled lane: {msg}"
        );

        // NON-VACUITY: the same dataset with real train rows FITS, so the error
        // above is a property of the empty slice and not of a dead fixture.
        // `split(40).0` rather than `(0..20)`: `toy` lays out targets then
        // decoys, so the first 20 rows are 20 targets and no decoys, and a
        // single-class fit is a second thing that could be going right or wrong
        // here. The even-row slice is class-balanced and costs nothing.
        let (train, _) = split(40);
        assert!(MlpFoldModel::fit(&cfg(7), &data, 2, &train, &[]).is_ok());
    }

    // --------------------------------------------------- inner validation

    /// The carve is a PURE FUNCTION OF `train.len()` — no RNG, no data, no
    /// clock — which is the whole reason it can sit inside a fit that has to
    /// reproduce bit-for-bit.
    #[test]
    fn the_inner_validation_carve_is_a_deterministic_partition() {
        // Under the floor: everything trains, nothing is carved.
        for len in [0usize, 1, 31, MIN_TRAIN_ROWS_FOR_INNER_VAL - 1] {
            let (fit, val) = inner_val_split(len);
            assert_eq!(fit, (0..len).collect::<Vec<_>>(), "len {len}");
            assert!(
                val.is_empty(),
                "len {len}: nothing may be carved under the floor"
            );
        }

        for len in [
            MIN_TRAIN_ROWS_FOR_INNER_VAL,
            MIN_TRAIN_ROWS_FOR_INNER_VAL + 1,
            201,
            1000,
            76_543,
        ] {
            let (fit, val) = inner_val_split(len);

            // Same call, same answer — and the SECOND call is what a rerun of
            // the same fit makes.
            assert_eq!(
                (fit.clone(), val.clone()),
                inner_val_split(len),
                "len {len}"
            );

            // A partition: disjoint, and together exactly `0..len`.
            let mut all: Vec<usize> = fit.iter().chain(val.iter()).copied().collect();
            all.sort_unstable();
            assert_eq!(all, (0..len).collect::<Vec<_>>(), "len {len}");

            // The stride, stated as the rule and not as a row count: every
            // position congruent to `STRIDE - 1` is held out and no other is.
            assert!(
                val.iter()
                    .all(|k| k % INNER_VAL_STRIDE == INNER_VAL_STRIDE - 1),
                "len {len}: a carved position is off the stride"
            );
            assert!(
                fit.iter()
                    .all(|k| k % INNER_VAL_STRIDE != INNER_VAL_STRIDE - 1),
                "len {len}: a stride position was left in the training set"
            );
            assert_eq!(val.len(), len / INNER_VAL_STRIDE, "len {len}: 20%, exactly");
            assert!(
                val.len() >= MIN_INNER_VAL_ROWS,
                "len {len}: the floor is on the VALIDATION SIZE, and it was not met"
            );
            assert!(fit.contains(&0), "position 0 must train");
        }
    }

    /// The carved rows never reach the optimizer.
    ///
    /// Asserted through the row count the fit actually trained on, because it is
    /// the only externally visible statement of the property. The tempting
    /// observational test — perturb the carved rows and check the scores do not
    /// move — cannot be written here: their FEATURE values feed the
    /// [`ColumnTransform`] by design (see [`FoldModel::fit`]), and perturbing
    /// their LABELS legitimately moves the stopping decision and hence the
    /// scores.
    #[test]
    fn inner_validation_rows_are_withheld_from_the_optimizer() {
        let (feat, y) = toy(200, 7, 0.0);
        let data = dataset(feat, 3, y);

        // 200 train rows, comfortably over the floor: 40 carved, 160 trained.
        let (train, _) = split(400);
        assert_eq!(train.len(), 200);
        let m = MlpFoldModel::fit(&cfg(7), &data, 0, &train, &[]).unwrap();
        assert_eq!(m.n_train_rows(), 160, "20% of `train` must be withheld");
        assert!(m.outcome().best_val_loss.is_some());

        // ...but an OUTER validation slice is used as-is and costs no training
        // rows, which is the reason it is preferred over the carve.
        let held: Vec<usize> = (1..400).step_by(2).collect();
        let outer = MlpFoldModel::fit(&cfg(7), &data, 0, &train, &held).unwrap();
        assert_eq!(outer.n_train_rows(), 200);
        assert!(outer.outcome().best_val_loss.is_some());

        // Under the floor, no carve and no early stopping: the full budget runs
        // on every row, exactly as it did before early stopping existed.
        let small: Vec<usize> = (0..100).step_by(2).collect();
        let tiny = MlpFoldModel::fit(&cfg(7), &data, 0, &small, &[]).unwrap();
        assert_eq!(tiny.n_train_rows(), small.len());
        assert_eq!(tiny.outcome().epochs_run, cfg(7).epochs);
        assert_eq!(tiny.outcome().best_epoch, None);

        // And with the knob off there is no carve at any size.
        let off = MlpConfig {
            early_stopping_patience: None,
            ..cfg(7)
        };
        let m_off = MlpFoldModel::fit(&off, &data, 0, &train, &[]).unwrap();
        assert_eq!(m_off.n_train_rows(), 200);
        assert_eq!(m_off.outcome().epochs_run, off.epochs);
        assert!(!m_off.outcome().restored);
    }

    #[test]
    fn progress_reports_each_completed_epoch_including_early_stopping() {
        use std::sync::atomic::{
            AtomicUsize,
            Ordering,
        };

        let (feat, y) = toy(200, 17, 0.0);
        let data = dataset(feat, 3, y);
        let (train, held) = split(400);
        let ticks = AtomicUsize::new(0);
        let advance = || {
            ticks.fetch_add(1, Ordering::Relaxed);
        };

        let model =
            MlpFoldModel::fit_with_progress(&cfg(17), &data, 0, &train, &held, &advance).unwrap();

        assert_eq!(ticks.load(Ordering::Relaxed), model.outcome().epochs_run);
        assert!(model.outcome().epochs_run <= cfg(17).epochs);
    }

    /// [`ColumnTransform::check_clean`]'s report must actually fire, i.e. the
    /// warn branch in [`FoldModel::predict`] must be reachable.
    ///
    /// It needs a column that is CLEAN across the train rows and non-finite only
    /// on a row being scored. That is an awkward fixture to hit by accident, and
    /// no other test in this file does: every one of them either NaNs the column
    /// in the train rows too — which makes it `Missable`, and `check_clean` skips
    /// those on purpose, since a missable column has an `_isna` companion to carry
    /// the missingness — or leaves the column entirely dead, which culls it. So
    /// `dirty` was empty everywhere and this path was never executed.
    ///
    /// Why it matters: on such a column the NaN imputes to the TRAIN MEAN with no
    /// flag, so the row is scored as if it were average on that feature. That is a
    /// train/score distribution mismatch worth telling an operator about, which is
    /// why the report is a returned value rather than a `debug_assert!` — the
    /// production build is exactly where it would happen.
    #[test]
    fn a_nan_only_in_a_scored_row_is_reported_rather_than_silently_imputed() {
        for seed in [7u64, 13, 42] {
            let (mut feat, y) = toy(150, seed, 0.0);
            let (train, held) = split(300);

            // Column 2, on TWO HELD rows only — one per class, so the fixture is
            // not also testing single-class behaviour.
            let dirty_rows = [held[0], held[1]];
            for &i in &dirty_rows {
                feat[i * 3 + 2] = f64::NAN;
            }
            let data = dataset(feat, 3, y);

            let model = MlpFoldModel::fit(&cfg(seed), &data, 0, &train, &[]).unwrap();
            // THE PREMISE: column 2 survived the cull and the fit saw it as
            // CLEAN, so it has no `_isna` companion to carry the missingness.
            assert!(
                model.transform().culled().is_empty(),
                "seed {seed}: a culled column has no clean/dirty distinction to report"
            );
            assert_eq!(
                model.transform().width(),
                3,
                "seed {seed}: three survivors and NO companion — a companion would mean the \
                 fit saw the column as missable, which is the case `check_clean` skips"
            );

            // THE ASSERTION: the branch's condition is true, so the warn in
            // `predict` runs — and it names the offending column.
            let dirty = model.transform().check_clean(
                &gather(&data, &held, 3),
                &(0..held.len()).collect::<Vec<_>>(),
            );
            assert_eq!(
                dirty,
                &[2],
                "seed {seed}: a NaN that only ever appears in a scored row must be reported"
            );

            // ...and scoring still completes, with finite values: the report is a
            // diagnostic, not a failure.
            let scores = model.predict(&data, &held).unwrap();
            assert_eq!(scores.len(), held.len());
            assert!(
                scores.iter().all(|s| s.is_finite()),
                "seed {seed}: the imputed rows must still score finitely"
            );

            // NON-VACUITY for the report: the SAME transform over rows without the
            // NaN reports nothing, so `dirty` above is about those two rows.
            let clean_rows: Vec<usize> = held
                .iter()
                .copied()
                .filter(|i| !dirty_rows.contains(i))
                .collect();
            assert!(
                model
                    .transform()
                    .check_clean(
                        &gather(&data, &clean_rows, 3),
                        &(0..clean_rows.len()).collect::<Vec<_>>()
                    )
                    .is_empty(),
                "seed {seed}: rows with no NaN must not be reported, else the assertion \
                 above holds for any input"
            );
        }
    }

    /// An outer `val` slice under [`MIN_INNER_VAL_ROWS`] is IGNORED, exactly as
    /// an under-floor carve would be.
    ///
    /// The floor is on the SIZE of the held-out set — see [`MIN_INNER_VAL_ROWS`]
    /// — so it cannot be a property of only one of the two ways the rows arrive.
    /// The concrete caller: `CrossValidatedScorer` hands `fold_rows[f + 1]`, which
    /// on this crate's 90-row fixtures is 30 rows, and a patience rule reading a
    /// 30-row loss stops on which rows landed in it.
    ///
    /// Asserted as BIT-IDENTITY with the `val = &[]` fit rather than as a row
    /// count alone: "ignored" means the fit is the one it would have been with no
    /// outer slice at all, down to the epoch kept.
    ///
    /// Several seeds because the fits being compared are real training runs, and
    /// a single init could agree by luck on a fixture where the two paths in fact
    /// diverge — the same reason every other learning test here sweeps.
    #[test]
    fn an_outer_validation_slice_under_the_floor_is_ignored() {
        for seed in [7u64, 13, 42] {
            let (feat, y) = toy(200, seed, 0.0);
            let data = dataset(feat, 3, y);
            let (train, held) = split(400);
            assert_eq!(train.len(), 200);

            // `held` is `toy`'s layout thinned: targets first, decoys second, so
            // its LEADING rows are all one class. Stride over it instead of
            // taking a prefix, or the slice under test is a target-only set no
            // stopping rule could learn from — the assertions here are about
            // whether the slice is USED, and they would hold vacuously on a
            // degenerate one.
            let stride = held.len() / MIN_INNER_VAL_ROWS;
            let spread: Vec<usize> = held
                .iter()
                .copied()
                .step_by(stride)
                .take(MIN_INNER_VAL_ROWS)
                .collect();
            assert_eq!(spread.len(), MIN_INNER_VAL_ROWS);
            assert!(
                spread.iter().any(|&i| data.is_decoy(i))
                    && spread.iter().any(|&i| !data.is_decoy(i)),
                "seed {seed}: the validation slice must carry both classes"
            );

            let none = MlpFoldModel::fit(&cfg(seed), &data, 0, &train, &[]).unwrap();
            let under: Vec<usize> = spread
                .iter()
                .copied()
                .take(MIN_INNER_VAL_ROWS - 1)
                .collect();
            let floored = MlpFoldModel::fit(&cfg(seed), &data, 0, &train, &under).unwrap();
            assert_eq!(
                floored.n_train_rows(),
                160,
                "seed {seed}: an under-floor outer slice must fall through to the inner carve"
            );
            assert_eq!(
                floored.predict(&data, &held).unwrap(),
                none.predict(&data, &held).unwrap(),
                "seed {seed}: an ignored outer slice must give the `val = &[]` fit exactly"
            );
            assert_eq!(floored.outcome(), none.outcome(), "seed {seed}");

            // NON-VACUITY: one more row and the slice IS used — the floor is at
            // `MIN_INNER_VAL_ROWS`, not "small slices are ignored".
            let used = MlpFoldModel::fit(&cfg(seed), &data, 0, &train, &spread).unwrap();
            assert_eq!(
                used.n_train_rows(),
                200,
                "seed {seed}: a slice AT the floor is usable, so no carve and every row trains"
            );
            assert!(used.outcome().best_val_loss.is_some(), "seed {seed}");
        }
    }
}
