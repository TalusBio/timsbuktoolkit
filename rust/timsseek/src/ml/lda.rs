//! Shrinkage Linear Discriminant Analysis rescorer.
//!
//! Sage-style two-class Fisher LDA (github.com/lazear/sage), adapted for the
//! timsseek feature set. Differences from Sage:
//!   * Features are z-standardized (mean/std over the training rows) before the
//!     fit. Sage skips this because exact Fisher LDA is scale-invariant; we add
//!     ridge shrinkage below, which is NOT scale-invariant, so standardization
//!     is required to make the `lambda * I` term comparable across features.
//!   * The within-class scatter is ridge-regularized (`Sw + lambda*I`). The
//!     timsseek feature set (~108 dims) is heavily collinear (main_score family,
//!     the 7 per-ion error dims + their mean, 22 AA counts, ...), which makes
//!     the raw `Sw` near-singular. Shrinkage keeps the linear solve stable
//!     without hand-pruning the feature set, so LDA sees the exact same inputs
//!     as the GBM path (clean apples-to-apples bench).
//!   * NaN/non-finite feature values are imputed to the column mean (i.e. 0
//!     after standardization). The GBM handles missing natively; LDA cannot.
//!
//! For two-class LDA the between-class scatter is rank-1 along `(mu_t - mu_d)`,
//! so the discriminant direction is `w = (Sw + lambda*I)^-1 (mu_t - mu_d)`.
//! Target rows project higher than decoy rows by construction.
//!
//! # Layout & scaling
//! The feature matrix is a single flat row-major `Vec<f64>` (`feat[i*d + j]`),
//! never a `Vec<Vec<f64>>` — at tens of millions of candidates the per-row
//! allocation and pointer-chasing dominate. The heavy cost is the within-class
//! scatter — `O(n * d^2)` (`d^2` FMAs/row) — **linear in the row count but with
//! a large constant**, parallelized (deterministically) over rayon. Nothing
//! here is super-linear in `n`. Skew-taming is done at feature-emit time by the
//! block grammar (log2/ln1p transforms); the linear lane only admits fields
//! that are approx-Gaussian after their declared transform, so no rank-based
//! gaussianization step runs here.

use std::time::Instant;

#[cfg(feature = "rayon")]
use rayon::prelude::*;

/// Fraction of `mean(diag(Sw))` added to the diagonal as ridge shrinkage.
/// Small enough to barely perturb a well-conditioned problem, large enough to
/// rescue a rank-deficient one.
const DEFAULT_SHRINKAGE: f64 = 1e-2;

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
        let (sum, sumsq, cnt) = col_finite_moments(feat, nrows, ncols);
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
        for c in 0..2 {
            let inv_nc = 1.0 / class_cnt[c] as f64;
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

    /// Score every row of a flat row-major matrix (parallel), into `out`.
    pub fn score_all(&self, feat: &[f64], nrows: usize, out: &mut [f64]) {
        let d = self.ncols;
        assert_eq!(feat.len(), nrows * d);
        assert_eq!(out.len(), nrows);
        #[cfg(feature = "rayon")]
        {
            out.par_iter_mut()
                .enumerate()
                .for_each(|(i, o)| *o = self.score(&feat[i * d..(i + 1) * d]));
        }
        #[cfg(not(feature = "rayon"))]
        {
            for (i, o) in out.iter_mut().enumerate() {
                *o = self.score(&feat[i * d..(i + 1) * d]);
            }
        }
    }

    /// Discriminant weights in standardized space, `|coef|`-interpretable as
    /// feature importance.
    pub fn coef(&self) -> &[f64] {
        &self.coef
    }
}

/// Per-column finite sum / sum-of-squares / count, reduced over fixed row
/// chunks (deterministic summation order).
fn col_finite_moments(feat: &[f64], nrows: usize, ncols: usize) -> (Vec<f64>, Vec<f64>, Vec<u64>) {
    let chunk = CHUNK_ROWS * ncols;
    let fold = |acc: &mut (Vec<f64>, Vec<f64>, Vec<u64>), block: &[f64]| {
        for row in block.chunks_exact(ncols) {
            for j in 0..ncols {
                let v = row[j];
                if v.is_finite() {
                    acc.0[j] += v;
                    acc.1[j] += v * v;
                    acc.2[j] += 1;
                }
            }
        }
    };
    let zero = || (vec![0.0f64; ncols], vec![0.0f64; ncols], vec![0u64; ncols]);

    #[cfg(feature = "rayon")]
    let partials: Vec<(Vec<f64>, Vec<f64>, Vec<u64>)> = feat
        .par_chunks(chunk)
        .map(|block| {
            let mut acc = zero();
            fold(&mut acc, block);
            acc
        })
        .collect();
    #[cfg(not(feature = "rayon"))]
    let partials: Vec<(Vec<f64>, Vec<f64>, Vec<u64>)> = feat
        .chunks(chunk)
        .map(|block| {
            let mut acc = zero();
            fold(&mut acc, block);
            acc
        })
        .collect();

    let mut out = zero();
    for (s, sq, c) in &partials {
        for j in 0..ncols {
            out.0[j] += s[j];
            out.1[j] += sq[j];
            out.2[j] += c[j];
        }
    }
    let _ = nrows;
    out
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
    let chunk_rows = CHUNK_ROWS;
    let std_at = |row: &[f64], j: usize| {
        let v = row[j];
        if v.is_finite() {
            (v - mean[j]) * inv_std[j]
        } else {
            0.0
        }
    };
    let per_chunk = |ci: usize| {
        let r0 = ci * chunk_rows;
        let r1 = ((ci + 1) * chunk_rows).min(nrows);
        let mut csum = [vec![0.0f64; ncols], vec![0.0f64; ncols]];
        let mut ccnt = [0u64; 2];
        for i in r0..r1 {
            let row = &feat[i * ncols..(i + 1) * ncols];
            let c = if is_decoy[i] { 0 } else { 1 };
            for j in 0..ncols {
                csum[c][j] += std_at(row, j);
            }
            ccnt[c] += 1;
        }
        (csum, ccnt)
    };
    let nchunks = nrows.div_ceil(chunk_rows);

    #[cfg(feature = "rayon")]
    let partials: Vec<([Vec<f64>; 2], [u64; 2])> =
        (0..nchunks).into_par_iter().map(per_chunk).collect();
    #[cfg(not(feature = "rayon"))]
    let partials: Vec<([Vec<f64>; 2], [u64; 2])> = (0..nchunks).map(per_chunk).collect();

    let mut class_sum = [vec![0.0f64; ncols], vec![0.0f64; ncols]];
    let mut class_cnt = [0u64; 2];
    for (csum, ccnt) in &partials {
        for c in 0..2 {
            for j in 0..ncols {
                class_sum[c][j] += csum[c][j];
            }
            class_cnt[c] += ccnt[c];
        }
    }
    (class_sum, class_cnt)
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
    let chunk_rows = CHUNK_ROWS;
    let dd = ncols * ncols;
    let per_chunk = |ci: usize| -> Vec<f64> {
        let r0 = ci * chunk_rows;
        let r1 = ((ci + 1) * chunk_rows).min(nrows);
        let mut out = vec![0.0f64; 2 * dd];
        let mut centered = vec![0.0f64; ncols];
        for i in r0..r1 {
            let row = &feat[i * ncols..(i + 1) * ncols];
            let c = if is_decoy[i] { 0 } else { 1 };
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
        out
    };
    let nchunks = nrows.div_ceil(chunk_rows);

    #[cfg(feature = "rayon")]
    let partials: Vec<Vec<f64>> = (0..nchunks).into_par_iter().map(per_chunk).collect();
    #[cfg(not(feature = "rayon"))]
    let partials: Vec<Vec<f64>> = (0..nchunks).map(per_chunk).collect();

    let mut sw = vec![0.0f64; 2 * dd];
    for p in &partials {
        for e in 0..2 * dd {
            sw[e] += p[e];
        }
    }
    sw
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
    fn solve_identity() {
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
    fn nan_imputed_not_propagated() {
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
    }
}

pub const DEFAULT_LDA_SHRINKAGE: f64 = DEFAULT_SHRINKAGE;
