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
//!     data.** The RNG comes from [`MlpConfig::rng_for_fold`], i.e. from the
//!     configured seed and the `fold` argument [`FoldModel::fit`] is handed and
//!     nothing else — no call counter, no clock, no hash-map iteration order —
//!     so folds differ from one another while the same build on the same input
//!     produces bit-identical scores.

use crate::ml::cv::{
    FoldDataset,
    FoldModel,
    fold_weights,
};
use crate::ml::mlp::{
    Adam,
    ColumnTransform,
    Mlp,
    MlpConfig,
    Tensor,
};
use std::cell::RefCell;

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
    /// abort in `qvalues::rescore_mlp_lane` prints only what the error says.
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
        }
    }
}

impl std::error::Error for MlpFoldError {}

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
    /// Mean loss of the final training epoch. Diagnostics only — `fit` logs it,
    /// and [`MlpFoldModel::final_loss`] hands it to a failing test's message.
    #[cfg_attr(not(test), allow(dead_code))]
    final_loss: f32,
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
}

/// Gather `rows` out of `data` into a compact row-major slab, in the given
/// order. The slab's row `k` is `rows[k]`, so slab-local indices are `0..len`
/// — which is what [`ColumnTransform::fit`] and [`ColumnTransform::check_clean`]
/// want, since both index `feat` by their `rows` argument.
///
/// Materializing rather than streaming is deliberate: `FoldDataset` only offers
/// row-at-a-time access and the transform needs two passes. The copy is a fixed
/// cost per fold and correctness comes first.
fn gather<D: FoldDataset>(data: &D, rows: &[usize], ncols: usize) -> Vec<f64> {
    let mut feat = vec![0.0f64; rows.len() * ncols];
    for (k, &i) in rows.iter().enumerate() {
        data.get_values(i, &mut feat[k * ncols..(k + 1) * ncols]);
    }
    feat
}

impl FoldModel for MlpFoldModel {
    type Config = MlpConfig;
    type Error = MlpFoldError;

    /// `val` is DELIBERATELY IGNORED. Training runs a fixed `cfg.epochs` with
    /// no early stopping, and [`FoldModel`] documents that a model without
    /// early stopping ignores the validation slice. Callers may pass `&[]`.
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
    fn fit<D: FoldDataset>(
        cfg: &MlpConfig,
        data: &D,
        fold: usize,
        train: &[usize],
        val: &[usize],
    ) -> Result<Self, MlpFoldError> {
        let _ = val;
        let names = data.column_names();
        let ncols = names.len();
        let feat = gather(data, train, ncols);
        // Slab-local indices: `feat` holds the train rows and nothing else, so
        // no held-out value is even reachable from here.
        let local: Vec<usize> = (0..train.len()).collect();
        let transform = ColumnTransform::fit(&feat, ncols, &local);

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

        // Labels: y = 1.0 target / 0.0 decoy, the GBM's convention.
        let mut x = Tensor::new(train.len(), width);
        let mut y = Vec::with_capacity(train.len());
        for (k, &i) in train.iter().enumerate() {
            transform.apply(&feat[k * ncols..(k + 1) * ncols], x.row_mut(k));
            y.push(if data.is_decoy(i) { 0.0 } else { 1.0 });
        }
        // Sample weights come from [`fold_weights`], the ONE definition of the
        // decoy-1.0 / target-0.5 convention. Re-inlining it here is how the MLP
        // and the GBM would come to weight their classes differently without any
        // test noticing — the two models would just be trained on different
        // objectives.
        let w: Vec<f32> = fold_weights(data, train)
            .into_iter()
            .map(|v| v as f32)
            .collect();

        // Seeded from (config seed, fold index) alone — see the module docs.
        let mut rng = cfg.rng_for_fold(fold);
        let mut net = Mlp::feedforward(width, &cfg.hidden, &mut rng);
        let mut opt = Adam::new(cfg.lr).with_weight_decay(cfg.weight_decay);
        let final_loss = net.train(cfg, &x, &y, &w, &mut opt, &mut rng);
        tracing::debug!(
            "MLP fold {}: trained on {} rows x {} inputs ({} lane columns), final epoch loss {}",
            fold,
            train.len(),
            width,
            ncols,
            final_loss,
        );

        Ok(MlpFoldModel {
            net: RefCell::new(net),
            transform,
            final_loss,
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

        let feat = gather(data, rows, ncols);
        let local: Vec<usize> = (0..rows.len()).collect();

        // A column the TRAIN rows never saw a non-finite value in has no
        // `_isna` companion, so a NaN arriving here is imputed to the train
        // mean and becomes indistinguishable from an average row. That is a
        // train/score distribution mismatch worth knowing about, and it is
        // reported at runtime rather than asserted: `debug_assert!` compiles
        // out in release, which is exactly where a production run would hit it.
        let dirty = self.transform.check_clean(&feat, &local);
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

        let mut x = Tensor::new(rows.len(), self.transform.width());
        for k in 0..rows.len() {
            self.transform
                .apply(&feat[k * ncols..(k + 1) * ncols], x.row_mut(k));
        }

        let mut net = self.net.borrow_mut();
        let out = net.forward(&x, false);
        Ok((0..rows.len()).map(|i| out.row(i)[0] as f64).collect())
    }

    /// `sum_o |W1[j, o]|` over the first layer, mapped back to lane indices.
    ///
    /// NEVER `NAN`, i.e. never "unreported" under the [`FoldModel::importance`]
    /// contract: the first layer has a weight for every input, so every lane
    /// column that survived the cull has a measurement.
    ///
    /// **Culled columns report `0.0`, not `NAN`.** A culled column genuinely
    /// contributed nothing — it has no input to weigh — and "this column was
    /// culled" is exactly what an operator reads the sidecar to find; it is the
    /// machine-readable half of the warn `fit` emits. `NAN` would drop those
    /// rows at the `fold_feature_stats` boundary and delete the cull report.
    ///
    /// `_isna` companions have no lane of their own, so each one's weight is
    /// summed into the column it flags: the importance of a missable column is
    /// the weight the net puts on its value plus the weight it puts on its
    /// missingness, which is the one number that answers "how much does this
    /// column matter".
    fn importance(&self) -> Vec<f32> {
        let raw = self.net.borrow().input_importance();
        // Every transformed input is either a standardized column or a
        // companion, so with the right length the loop below places all of
        // them; a short/long vector is the only way one could go missing.
        debug_assert_eq!(
            raw.len(),
            self.transform.width(),
            "input_importance must report one value per transformed input"
        );
        let mut out = vec![0.0f32; self.transform.ncols_lane()];
        for (k, &v) in raw.iter().enumerate() {
            let lane = self
                .transform
                .lane_of_input(k)
                .or_else(|| self.transform.isna_lane_of_input(k));
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
            "folding companions into their parent lane must not drop or double any weight"
        );
        out
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
    use std::sync::Arc;

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
    #[test]
    fn refitting_the_same_fold_reproduces_the_scores_exactly() {
        let (feat, y) = toy(100, 42, 0.0);
        let data = dataset(feat, 3, y);
        let (train, held) = split(200);

        let a = MlpFoldModel::fit(&cfg(42), &data, 0, &train, &[]).unwrap();
        let b = MlpFoldModel::fit(&cfg(42), &data, 0, &train, &[]).unwrap();
        assert_eq!(
            a.predict(&data, &held).unwrap(),
            b.predict(&data, &held).unwrap()
        );
        assert_eq!(a.importance(), b.importance());

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
    /// diagnostic the abort in `qvalues::rescore_mlp_lane` can print, so a fit
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
}
