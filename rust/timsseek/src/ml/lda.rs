//! Sage-style two-class Fisher LDA (github.com/lazear/sage) with ridge
//! shrinkage: `w = (Sw + lambda*I)^-1 (mu_t - mu_d)`, targets projecting high.
//!
//! Fits the LINEAR lane only (`LINEAR_NCOLS`, 101 columns today). That lane is
//! NOT the GBM's input — the GBM sees the ALL lane (`ALL_NCOLS`, 128), including
//! the 22-wide `sequence_counts` block, which has no linear lane and never
//! reaches LDA.
//!
//! Features are z-standardized before the fit (Fisher LDA is scale-invariant,
//! but the `lambda * I` term is not), and non-finite values are imputed to the
//! column mean, i.e. 0 post-standardization — LDA cannot take missing natively.

use crate::ml::cv::{
    FoldDataset,
    FoldModel,
};
use crate::utils::maybe_par::chunked_fold_reduce;
use std::time::Instant;

/// Fraction of `mean(diag(Sw))` added to the diagonal as ridge shrinkage.
/// Small enough to barely perturb a well-conditioned problem, large enough to
/// rescue a rank-deficient one.
pub const DEFAULT_SHRINKAGE: f64 = 1e-2;

/// Row-chunk size for the parallel reductions. Fixed so chunk boundaries — and
/// therefore the summation order of the partial accumulators — are identical
/// across runs, keeping the fit bitwise-deterministic despite parallelism.
const CHUNK_ROWS: usize = 65_536;

pub struct LdaModel {
    /// Discriminant direction in standardized feature space.
    coef: Vec<f64>,
    /// Per-feature standardization mean (finite-only).
    mean: Vec<f64>,
    /// Per-feature standardization std (finite-only, floored so near-constant
    /// features standardize to ~0 rather than exploding).
    inv_std: Vec<f64>,
    ncols: usize,
}

impl LdaModel {
    /// Fit LDA from a flat row-major feature matrix (`feat[i*ncols + j]`) and
    /// boolean decoy labels (`is_decoy[i]`).
    ///
    /// Returns `None` if either class is empty or the linear solve is singular
    /// even after shrinkage.
    pub fn fit(
        feat: &[f64],
        nrows: usize,
        ncols: usize,
        is_decoy: &[bool],
        shrinkage: f64,
    ) -> Option<LdaModel> {
        assert_eq!(is_decoy.len(), nrows);
        assert_eq!(feat.len(), nrows * ncols);
        if nrows == 0 || ncols == 0 {
            return None;
        }

        // --- Standardization stats (finite values only), parallel reduce ---
        let t = Instant::now();
        let (sum, sumsq, cnt) = col_finite_moments(feat, ncols);
        let mut mean = vec![0.0f64; ncols];
        let mut inv_std = vec![1.0f64; ncols];
        for j in 0..ncols {
            if cnt[j] > 0 {
                let n = cnt[j] as f64;
                let m = sum[j] / n;
                let var = (sumsq[j] / n - m * m).max(0.0);
                let std = var.sqrt();
                mean[j] = m;
                inv_std[j] = if std > 1e-12 { 1.0 / std } else { 0.0 };
            }
        }

        // --- Pass 1: per-class means in standardized space (parallel) ---
        let (class_sum, class_cnt) = class_std_sums(feat, nrows, ncols, is_decoy, &mean, &inv_std);
        if class_cnt[0] == 0 || class_cnt[1] == 0 {
            return None;
        }
        let class_mean: [Vec<f64>; 2] = std::array::from_fn(|c| {
            let n = class_cnt[c] as f64;
            (0..ncols).map(|j| class_sum[c][j] / n).collect()
        });
        eprintln!(
            "  LDA: standardized + class means ({nrows} rows x {ncols} feats) in {:.2?}",
            t.elapsed()
        );

        // --- Pass 2: pooled within-class scatter Sw (parallel, D x D) ---
        // Sw = sum_c (1/n_c) sum_{i in c} (z_i - mu_c)(z_i - mu_c)^T
        let t = Instant::now();
        let mut sw =
            within_class_scatter(feat, nrows, ncols, is_decoy, &mean, &inv_std, &class_mean);
        for (c, &cnt) in class_cnt.iter().enumerate() {
            let inv_nc = 1.0 / cnt as f64;
            let off = c * ncols * ncols;
            for e in &mut sw[off..off + ncols * ncols] {
                *e *= inv_nc;
            }
        }
        // Fold the two per-class scatters into one within-class matrix.
        let mut sw_within = vec![0.0f64; ncols * ncols];
        for e in 0..ncols * ncols {
            sw_within[e] = sw[e] + sw[ncols * ncols + e];
        }
        eprintln!("  LDA: within-class scatter in {:.2?}", t.elapsed());

        // --- Ridge shrinkage: Sw += lambda_eff * I ---
        let mut diag_sum = 0.0f64;
        for j in 0..ncols {
            diag_sum += sw_within[j * ncols + j];
        }
        let lambda_eff = shrinkage * (diag_sum / ncols as f64).max(1e-12);
        for j in 0..ncols {
            sw_within[j * ncols + j] += lambda_eff;
        }

        // --- Solve Sw * w = (mu_t - mu_d) ---
        let mu_diff: Vec<f64> = (0..ncols)
            .map(|j| class_mean[1][j] - class_mean[0][j])
            .collect();
        let coef = solve_gauss(sw_within, mu_diff, ncols)?;
        if !coef.iter().all(|c| c.is_finite()) {
            return None;
        }

        Some(LdaModel {
            coef,
            mean,
            inv_std,
            ncols,
        })
    }

    /// Project one raw (un-standardized) feature row onto the discriminant.
    pub fn score(&self, row: &[f64]) -> f64 {
        debug_assert_eq!(row.len(), self.ncols);
        let mut acc = 0.0;
        // `self.ncols` (not `row.len()`) is the authoritative width across the four
        // parallel arrays `row`/`mean`/`inv_std`/`coef`; a zip would need a 4-way nest.
        #[allow(clippy::needless_range_loop)]
        for j in 0..self.ncols {
            let v = row[j];
            let z = if v.is_finite() {
                (v - self.mean[j]) * self.inv_std[j]
            } else {
                0.0
            };
            acc += self.coef[j] * z;
        }
        acc
    }
}

/// Everything [`LdaModel`] needs that is not data. One field today; a struct
/// rather than a bare `f64` so [`FoldModel::Config`] has somewhere to grow.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LdaConfig {
    /// Ridge shrinkage, see [`DEFAULT_SHRINKAGE`].
    pub shrinkage: f64,
}

impl Default for LdaConfig {
    fn default() -> Self {
        Self {
            shrinkage: DEFAULT_SHRINKAGE,
        }
    }
}

/// The only way an LDA fit fails. [`LdaModel::fit`] returns `None` for both
/// causes without distinguishing them, so this enum has one variant rather
/// than inventing a distinction the fit does not actually report.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LdaError {
    /// Either class had no rows, or the within-class scatter stayed singular
    /// even after shrinkage.
    SingularOrEmptyClass,
}

impl std::fmt::Display for LdaError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LdaError::SingularOrEmptyClass => f.write_str("singular or empty class"),
        }
    }
}

impl std::error::Error for LdaError {}

/// [`FoldModel`] adapter over the closed-form fit above.
///
/// Gathers its train rows out of the dataset into the flat row-major buffer
/// [`LdaModel::fit`] wants; `predict` reads scored rows the same way. Neither
/// touches a row outside the index slice it was handed, so leak-freedom is
/// entirely the caller's partition and never this impl's.
impl FoldModel for LdaModel {
    type Config = LdaConfig;
    type Error = LdaError;

    /// `val` is DELIBERATELY IGNORED, and that is not an oversight: this LDA is
    /// closed-form (one linear solve), so there is no iteration to early-stop
    /// and no use for a validation slice. [`FoldModel`] documents that a model
    /// without early stopping is expected to ignore it. Callers may pass an
    /// empty slice.
    ///
    /// `fold` is ignored for the same kind of reason: the fit is a
    /// deterministic linear solve with no initialization to seed, so the fold's
    /// identity would have nothing to change. Two folds differ because they are
    /// fitted on different rows, and a rerun of one is bit-identical.
    fn fit<D: FoldDataset>(
        cfg: &LdaConfig,
        data: &D,
        fold: usize,
        train: &[usize],
        val: &[usize],
    ) -> Result<Self, LdaError> {
        let (_, _) = (fold, val);
        let ncols = data.column_names().len();
        let mut feat = Vec::with_capacity(train.len() * ncols);
        let mut is_decoy = Vec::with_capacity(train.len());
        let mut row = vec![0.0f64; ncols];
        for &i in train {
            data.get_values(i, &mut row);
            feat.extend_from_slice(&row);
            is_decoy.push(data.is_decoy(i));
        }
        LdaModel::fit(&feat, train.len(), ncols, &is_decoy, cfg.shrinkage)
            .ok_or(LdaError::SingularOrEmptyClass)
    }

    fn predict<D: FoldDataset>(&self, data: &D, rows: &[usize]) -> Result<Vec<f64>, LdaError> {
        let mut row = vec![0.0f64; self.ncols];
        Ok(rows
            .iter()
            .map(|&i| {
                data.get_values(i, &mut row);
                self.score(&row)
            })
            .collect())
    }

    /// `|coef|` — the discriminant weights live in standardized space, so their
    /// magnitudes are directly comparable across columns.
    ///
    /// NEVER `NAN`, i.e. never "unreported" under the
    /// [`FoldModel::importance`] contract: the solve looks at every column, so
    /// every column has a measurement. `LdaModel::fit` rejects a non-finite
    /// solution outright, so `coef` is finite by construction and this cannot
    /// emit the sentinel by accident.
    ///
    /// A column that is constant, or has no finite value at all, standardizes
    /// to all-zeros and therefore solves to exactly `0.0`. That zero is a
    /// RESULT, not a gap, and it reaches the sidecar — "this feature is dead in
    /// this fold" is one of the more useful things the sidecar can say.
    fn importance(&self) -> Vec<f32> {
        let out: Vec<f32> = self.coef.iter().map(|c| c.abs() as f32).collect();
        debug_assert!(
            out.iter().all(|v| v.is_finite()),
            "LdaModel::fit guarantees finite coefficients"
        );
        out
    }
}

/// Per-column finite sum / sum-of-squares / count, reduced over fixed row
/// chunks (deterministic summation order).
fn col_finite_moments(feat: &[f64], ncols: usize) -> (Vec<f64>, Vec<f64>, Vec<u64>) {
    chunked_fold_reduce(
        feat,
        CHUNK_ROWS * ncols,
        || (vec![0.0f64; ncols], vec![0.0f64; ncols], vec![0u64; ncols]),
        |acc, _ci, block| {
            for row in block.chunks_exact(ncols) {
                for (j, &v) in row.iter().enumerate() {
                    if v.is_finite() {
                        acc.0[j] += v;
                        acc.1[j] += v * v;
                        acc.2[j] += 1;
                    }
                }
            }
        },
        |out, (s, sq, c)| {
            for j in 0..ncols {
                out.0[j] += s[j];
                out.1[j] += sq[j];
                out.2[j] += c[j];
            }
        },
    )
}

/// Per-class standardized-feature sums + counts (parallel, deterministic).
fn class_std_sums(
    feat: &[f64],
    nrows: usize,
    ncols: usize,
    is_decoy: &[bool],
    mean: &[f64],
    inv_std: &[f64],
) -> ([Vec<f64>; 2], [u64; 2]) {
    debug_assert_eq!(feat.len(), nrows * ncols);
    let std_at = |row: &[f64], j: usize| {
        let v = row[j];
        if v.is_finite() {
            (v - mean[j]) * inv_std[j]
        } else {
            0.0
        }
    };
    chunked_fold_reduce(
        feat,
        CHUNK_ROWS * ncols,
        || ([vec![0.0f64; ncols], vec![0.0f64; ncols]], [0u64; 2]),
        |(csum, ccnt), ci, block| {
            // The block holds rows `ci * CHUNK_ROWS ..` of the matrix; the
            // labels are indexed globally, so recover the offset.
            let r0 = ci * CHUNK_ROWS;
            for (local, row) in block.chunks_exact(ncols).enumerate() {
                let c = if is_decoy[r0 + local] { 0 } else { 1 };
                for (j, s) in csum[c].iter_mut().enumerate() {
                    *s += std_at(row, j);
                }
                ccnt[c] += 1;
            }
        },
        |(class_sum, class_cnt), (csum, ccnt)| {
            for c in 0..2 {
                for j in 0..ncols {
                    class_sum[c][j] += csum[c][j];
                }
                class_cnt[c] += ccnt[c];
            }
        },
    )
}

/// Per-class within-class scatter (un-normalized), returned as two stacked
/// row-major `D x D` blocks: `[class0 (D*D), class1 (D*D)]`. Reduced over fixed
/// row chunks for determinism. This is the `O(n * d^2)` hot loop.
fn within_class_scatter(
    feat: &[f64],
    nrows: usize,
    ncols: usize,
    is_decoy: &[bool],
    mean: &[f64],
    inv_std: &[f64],
    class_mean: &[Vec<f64>; 2],
) -> Vec<f64> {
    debug_assert_eq!(feat.len(), nrows * ncols);
    let dd = ncols * ncols;
    chunked_fold_reduce(
        feat,
        CHUNK_ROWS * ncols,
        || vec![0.0f64; 2 * dd],
        |out, ci, block| {
            // Global row offset of this block, for the label lookup.
            let r0 = ci * CHUNK_ROWS;
            let mut centered = vec![0.0f64; ncols];
            for (local, row) in block.chunks_exact(ncols).enumerate() {
                let c = if is_decoy[r0 + local] { 0 } else { 1 };
                let mu = &class_mean[c];
                for j in 0..ncols {
                    let v = row[j];
                    let z = if v.is_finite() {
                        (v - mean[j]) * inv_std[j]
                    } else {
                        0.0
                    };
                    centered[j] = z - mu[j];
                }
                let base = c * dd;
                for j in 0..ncols {
                    let cj = centered[j];
                    if cj == 0.0 {
                        continue;
                    }
                    let rowbase = base + j * ncols;
                    for k in 0..ncols {
                        out[rowbase + k] += cj * centered[k];
                    }
                }
            }
        },
        |sw, p| {
            for e in 0..2 * dd {
                sw[e] += p[e];
            }
        },
    )
}

/// Solve `A x = b` for a row-major `n x n` matrix `A` via Gaussian elimination
/// with partial pivoting. Returns `None` if `A` is singular. `A` and `b` are
/// consumed (used as scratch).
fn solve_gauss(mut a: Vec<f64>, mut b: Vec<f64>, n: usize) -> Option<Vec<f64>> {
    for col in 0..n {
        // Partial pivot: largest magnitude in this column at/below the diagonal.
        let mut pivot = col;
        let mut best = a[col * n + col].abs();
        for r in (col + 1)..n {
            let v = a[r * n + col].abs();
            if v > best {
                best = v;
                pivot = r;
            }
        }
        if best < 1e-12 {
            return None;
        }
        if pivot != col {
            for k in 0..n {
                a.swap(col * n + k, pivot * n + k);
            }
            b.swap(col, pivot);
        }
        // Eliminate below.
        let diag = a[col * n + col];
        for r in (col + 1)..n {
            let factor = a[r * n + col] / diag;
            if factor == 0.0 {
                continue;
            }
            for k in col..n {
                a[r * n + k] -= factor * a[col * n + k];
            }
            b[r] -= factor * b[col];
        }
    }
    // Back-substitution.
    let mut x = vec![0.0f64; n];
    for col in (0..n).rev() {
        let mut acc = b[col];
        for k in (col + 1)..n {
            acc -= a[col * n + k] * x[k];
        }
        x[col] = acc / a[col * n + col];
    }
    Some(x)
}

#[cfg(test)]
mod test {
    use super::*;

    /// Row-major flat helper for tests.
    fn flat(rows: &[Vec<f64>]) -> (Vec<f64>, usize, usize) {
        let nrows = rows.len();
        let ncols = rows.first().map(|r| r.len()).unwrap_or(0);
        let mut v = Vec::with_capacity(nrows * ncols);
        for r in rows {
            v.extend_from_slice(r);
        }
        (v, nrows, ncols)
    }

    #[test]
    fn solve_gauss_matches_hand_solved_2x2() {
        // [2 1; 1 3] x = [3; 5]  =>  x = [0.8; 1.4] (det 5).
        let a = vec![2.0, 1.0, 1.0, 3.0];
        let b = vec![3.0, 5.0];
        let x = solve_gauss(a, b, 2).unwrap();
        assert!((x[0] - 0.8).abs() < 1e-9, "{x:?}");
        assert!((x[1] - 1.4).abs() < 1e-9, "{x:?}");
    }

    #[test]
    fn singular_returns_none() {
        let a = vec![1.0, 2.0, 2.0, 4.0];
        let b = vec![1.0, 2.0];
        assert!(solve_gauss(a, b, 2).is_none());
    }

    #[test]
    fn separable_two_class() {
        let mut rows = Vec::new();
        let mut labels = Vec::new();
        for i in 0..200 {
            let f = (i % 10) as f64 * 0.1;
            rows.push(vec![3.0 + f, 3.0 - f]);
            labels.push(false); // target
            rows.push(vec![0.0 + f, 0.0 - f]);
            labels.push(true); // decoy
        }
        let (feat, nrows, ncols) = flat(&rows);
        let lda = LdaModel::fit(&feat, nrows, ncols, &labels, DEFAULT_SHRINKAGE).unwrap();
        let mut t_sum = 0.0;
        let mut d_sum = 0.0;
        for (i, &d) in labels.iter().enumerate() {
            let s = lda.score(&feat[i * ncols..(i + 1) * ncols]);
            if d {
                d_sum += s;
            } else {
                t_sum += s;
            }
        }
        assert!(
            t_sum > d_sum,
            "targets {t_sum} should exceed decoys {d_sum}"
        );
    }

    #[test]
    fn nan_imputed_to_column_mean() {
        // Col 1's only finite values are 1.0 and 2.0, so its imputation target
        // (the finite-only column mean) is exactly 1.5.
        let rows = vec![
            vec![1.0, f64::NAN],
            vec![2.0, 1.0],
            vec![3.0, 2.0],
            vec![0.0, f64::NAN],
        ];
        let labels = vec![false, false, true, true];
        let (feat, nrows, ncols) = flat(&rows);
        let lda = LdaModel::fit(&feat, nrows, ncols, &labels, DEFAULT_SHRINKAGE).unwrap();

        for i in 0..nrows {
            assert!(lda.score(&feat[i * ncols..(i + 1) * ncols]).is_finite());
        }

        // The contract is not merely "finite": a NaN must score as if it held
        // the column mean. `is_finite` alone would also pass for impute-to-zero.
        const COL1_MEAN: f64 = 1.5;
        for row in [vec![1.0, f64::NAN], vec![0.0, f64::NAN]] {
            let imputed = vec![row[0], COL1_MEAN];
            assert_eq!(
                lda.score(&row),
                lda.score(&imputed),
                "NaN row {row:?} must score as its column-mean-imputed twin {imputed:?}"
            );
        }
        // Guard the premise: a DIFFERENT col-1 value must move the score, else
        // the equality above would be vacuous (e.g. a zeroed coefficient).
        assert_ne!(lda.score(&[1.0, COL1_MEAN]), lda.score(&[1.0, 5.0]));
    }
}
