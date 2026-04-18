//! Apex-local scoring functions for the METHODS.md Phase 3 scoring pipeline.
//!
//! These functions are computed once at specific apex locations, not per-cycle.
//! They operate on windows of the raw chromatogram data and implement the
//! 11 feature functions plus helper structs described in METHODS.md Sections 3.1-3.5.

use array2d::Array2D;
use timsquery::models::MzMajorIntensityArray;
use timsquery::traits::KeyLike;

use crate::models::query_item::linear_get;

/// Reusable scratch buffers for `coelution_gradient`. Owned by the scorer so
/// allocations happen once per rayon worker split and are reused across every
/// peptide that worker processes. `new()` pre-sizes capacity for a typical
/// peptide so the first hot-path call does no reallocation.
#[derive(Debug)]
pub struct CoelutionScratch {
    /// Rows: active fragments. Cols: coelution window length.
    coel_rows: Array2D<f32>,
    /// Rows: active fragments. Cols: gradient window length - 1.
    grad_rows: Array2D<f32>,
    /// One weight per active fragment in `coel_rows`.
    weights: Vec<f32>,
    /// First-differences buffer reused across fragments inside one call.
    diffs: Vec<f32>,
}

impl CoelutionScratch {
    const TYP_COEL_LEN: usize = 41;
    const TYP_DIFF_LEN: usize = 20;
    /// Typical peptide has up to ~16 fragments; coel window = 2*20+1 = 41 cycles,
    /// gradient diffs length = 2*10 = 20. Pre-sized to these so the first call
    /// in a worker does not reallocate.
    const TYP_FRAGS: usize = 16;

    pub fn new() -> Self {
        Self {
            coel_rows: Array2D::with_capacity(Self::TYP_COEL_LEN, Self::TYP_FRAGS),
            grad_rows: Array2D::with_capacity(Self::TYP_DIFF_LEN, Self::TYP_FRAGS),
            weights: Vec::with_capacity(Self::TYP_FRAGS),
            diffs: Vec::with_capacity(Self::TYP_DIFF_LEN),
        }
    }
}

impl Default for CoelutionScratch {
    fn default() -> Self {
        Self::new()
    }
}

/// Floor value for Scribe score when no signal is observed.
pub const SCRIBE_FLOOR: f32 = -100.0;

// ---------------------------------------------------------------------------
// Structs
// ---------------------------------------------------------------------------

/// The 11 apex-local features described in METHODS.md Section 3.4.
#[derive(Debug, Clone, Copy)]
pub struct ApexFeatures {
    pub peak_shape: f32,
    pub ratio_cv: f32,
    pub centered_apex: f32,
    pub precursor_coelution: f32,
    pub fragment_coverage: f32,
    pub precursor_apex_match: f32,
    pub xic_quality: f32,
    pub fragment_apex_agreement: f32,
    pub isotope_correlation: f32,
    pub gaussian_correlation: f32,
    pub per_frag_gaussian_corr: f32,
}

/// Result of the area-uniqueness calculation (METHODS.md Section 3.1).
#[derive(Debug, Clone, Copy)]
pub struct AreaUniquenessResult {
    pub au_score: f32,
}

/// Result of the coelution-gradient quality calculation (METHODS.md Section 3.2).
#[derive(Debug, Clone, Copy)]
pub struct CoelutionGradientResult {
    pub weighted_coelution: f32,
    pub gradient_consistency: f32,
    pub combined: f32,
}

/// Split product score computed from independent cosine and scribe apexes (METHODS.md Section 3.1).
#[derive(Debug, Clone, Copy)]
pub struct SplitProductScore {
    pub cosine_au: f32,
    pub cosine_cg: f32,
    pub scribe_au: f32,
    pub scribe_cg: f32,
    pub base_score: f32,
    pub cosine_weighted_coelution: f32,
    pub cosine_gradient_consistency: f32,
    pub scribe_weighted_coelution: f32,
    pub scribe_gradient_consistency: f32,
}

// ---------------------------------------------------------------------------
// Scoring weights: (offset, scale) pairs for each feature.
// Final score = base * product(offset + scale * feature_k)
// See METHODS.md Section 3.5.
// ---------------------------------------------------------------------------

pub const SCORING_WEIGHTS: [(f32, f32); 11] = [
    (1.0, 3.5),  // peak_shape
    (1.0, 3.0),  // ratio_cv
    (1.0, 4.5),  // centered_apex
    (1.0, 2.3),  // precursor_coelution
    (0.32, 1.0), // fragment_coverage
    (1.0, 8.0),  // precursor_apex_match
    (1.0, 7.8),  // xic_quality
    (0.43, 1.0), // fragment_apex_agreement
    (1.0, 2.6),  // isotope_correlation
    (0.27, 1.0), // gaussian_correlation
    (0.65, 1.0), // per_frag_gaussian_corr
];

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

/// Pearson correlation between two equal-length slices.
/// Returns 0.0 if either has zero variance.
fn pearson_correlation(a: &[f32], b: &[f32]) -> f32 {
    debug_assert_eq!(a.len(), b.len());
    let n = a.len();
    if n == 0 {
        return 0.0;
    }
    let mean_a: f32 = a.iter().sum::<f32>() / n as f32;
    let mean_b: f32 = b.iter().sum::<f32>() / n as f32;

    let mut cov = 0.0f32;
    let mut var_a = 0.0f32;
    let mut var_b = 0.0f32;
    for i in 0..n {
        let da = a[i] - mean_a;
        let db = b[i] - mean_b;
        cov += da * db;
        var_a += da * da;
        var_b += db * db;
    }
    let denom = (var_a * var_b).sqrt();
    if denom < 1e-12 {
        return 0.0;
    }
    cov / denom
}

/// Center a window (subtract mean) and normalize to unit length.
/// Writes into `out` (must be same length as `window`).
/// Returns `false` if the norm is zero (all-constant input), leaving `out` zeroed.
fn center_normalize(window: &[f32], out: &mut [f32]) -> bool {
    debug_assert_eq!(window.len(), out.len());
    let n = window.len();
    if n == 0 {
        return false;
    }
    let mean: f32 = window.iter().sum::<f32>() / n as f32;
    let mut norm_sq = 0.0f32;
    for i in 0..n {
        let v = window[i] - mean;
        out[i] = v;
        norm_sq += v * v;
    }
    let norm = norm_sq.sqrt();
    if norm < 1e-12 {
        for v in out.iter_mut() {
            *v = 0.0;
        }
        return false;
    }
    for v in out.iter_mut() {
        *v /= norm;
    }
    true
}

/// Estimate FWHM of a profile around the apex by descending from the peak
/// to find the half-maximum width on both sides.
fn estimate_fwhm(profile: &[f32], apex: usize) -> f32 {
    if profile.is_empty() {
        return 2.0;
    }
    let peak_val = profile[apex];
    if peak_val <= 0.0 {
        return 2.0;
    }
    let half_max = peak_val * 0.5;

    // Descend left
    let mut left_dist = 0.0f32;
    for i in (0..apex).rev() {
        if profile[i] <= half_max {
            // Linear interpolation between i and i+1
            let above = profile[i + 1];
            let below = profile[i];
            let frac = if (above - below).abs() > 1e-12 {
                (above - half_max) / (above - below)
            } else {
                0.5
            };
            left_dist = (apex - i - 1) as f32 + frac;
            break;
        }
        if i == 0 {
            left_dist = apex as f32;
        }
    }

    // Descend right
    let mut right_dist = 0.0f32;
    for i in (apex + 1)..profile.len() {
        if profile[i] <= half_max {
            let above = profile[i - 1];
            let below = profile[i];
            let frac = if (above - below).abs() > 1e-12 {
                (above - half_max) / (above - below)
            } else {
                0.5
            };
            right_dist = (i - apex - 1) as f32 + frac;
            break;
        }
        if i == profile.len() - 1 {
            right_dist = (profile.len() - 1 - apex) as f32;
        }
    }

    let fwhm = (left_dist + right_dist).max(1.0);
    fwhm
}

// ---------------------------------------------------------------------------
// Public functions
// ---------------------------------------------------------------------------

/// Area-uniqueness score around an apex (METHODS.md Section 3.1).
///
/// `signal` is the per-cycle profile (e.g. cosine_profile or scribe_profile).
/// `apex` is the argmax index within `signal`.
/// `hw` is the half-width of the peak window (typically 5).
///
/// AU = peak_area * (1 + 200 * peak_area / total_area)
pub fn area_uniqueness(signal: &[f32], apex: usize, hw: usize) -> AreaUniquenessResult {
    let n = signal.len();
    if n == 0 {
        return AreaUniquenessResult { au_score: 0.0 };
    }
    let lo = apex.saturating_sub(hw);
    let hi = (apex + hw + 1).min(n);

    let peak_area: f32 = signal[lo..hi].iter().sum();
    let total: f32 = signal.iter().sum();

    if total <= 0.0 {
        return AreaUniquenessResult { au_score: 0.0 };
    }

    let uniqueness = peak_area / total;
    let au_score = peak_area * (1.0 + 200.0 * uniqueness);

    AreaUniquenessResult { au_score }
}

/// Coelution-gradient quality at a given apex (METHODS.md Section 3.2).
///
/// Measures whether the fragment XICs co-elute and change in the same direction.
/// Weights use raw predicted intensity fractions (NOT sqrt-transformed).
///
/// Returns `combined = 1.0` if fewer than 2 active fragments in the window.
pub fn coelution_gradient<T: KeyLike>(
    fragments: &MzMajorIntensityArray<T, f32>,
    expected: &[(T, f32)],
    apex: usize,
    coelution_hw: usize,
    gradient_hw: usize,
    scratch: &mut CoelutionScratch,
) -> CoelutionGradientResult {
    let n_cycles = fragments.num_cycles();
    if n_cycles == 0 {
        return CoelutionGradientResult {
            weighted_coelution: 0.0,
            gradient_consistency: 0.0,
            combined: 1.0,
        };
    }

    let coel_lo = apex.saturating_sub(coelution_hw);
    let coel_hi = (apex + coelution_hw + 1).min(n_cycles);
    let coel_len = coel_hi - coel_lo;

    let grad_lo = apex.saturating_sub(gradient_hw);
    let grad_hi = (apex + gradient_hw + 1).min(n_cycles);
    let grad_len = grad_hi - grad_lo;
    let diff_len = grad_len.saturating_sub(1);

    let weight_sum: f32 = expected.iter().map(|(_, v)| *v).sum();
    if weight_sum <= 0.0 {
        return CoelutionGradientResult {
            weighted_coelution: 0.0,
            gradient_consistency: 0.0,
            combined: 1.0,
        };
    }

    // Upper bound on active fragments; we fill contiguous rows and track
    // `n_active_coel` / `n_active_grad` ourselves.
    let max_frags = fragments.num_ions();
    scratch.coel_rows.reset_with_value(coel_len, max_frags, 0.0);
    scratch.grad_rows.reset_with_value(diff_len, max_frags, 0.0);
    scratch.weights.clear();
    scratch.diffs.clear();
    scratch.diffs.resize(diff_len, 0.0);

    let mut n_active_coel: usize = 0;
    let mut n_active_grad: usize = 0;

    for ((key, _mz), chrom) in fragments.iter_mzs() {
        let y_hat = linear_get(expected, key).unwrap_or(0.0);
        if y_hat <= 0.0 {
            continue;
        }

        let coel_slice = &chrom[coel_lo..coel_hi];
        let coel_row = scratch
            .coel_rows
            .get_row_mut(n_active_coel)
            .expect("coel_rows pre-sized to max_frags");
        if center_normalize(coel_slice, coel_row) {
            scratch.weights.push(y_hat / weight_sum);
            n_active_coel += 1;
        }

        if diff_len >= 1 && grad_len >= 2 {
            let grad_slice = &chrom[grad_lo..grad_hi];
            for (i, w) in grad_slice.windows(2).enumerate() {
                scratch.diffs[i] = w[1] - w[0];
            }
            let grad_row = scratch
                .grad_rows
                .get_row_mut(n_active_grad)
                .expect("grad_rows pre-sized to max_frags");
            if center_normalize(&scratch.diffs, grad_row) {
                n_active_grad += 1;
            }
        }
    }

    if n_active_coel < 2 {
        return CoelutionGradientResult {
            weighted_coelution: 0.0,
            gradient_consistency: 0.0,
            combined: 1.0,
        };
    }

    // Weighted coelution: pairwise dot products (centered+normalized = correlations).
    let mut wcoel_num = 0.0f32;
    let mut wcoel_den = 0.0f32;
    for i in 0..n_active_coel {
        let row_i = scratch.coel_rows.get_row(i).unwrap();
        for j in (i + 1)..n_active_coel {
            let row_j = scratch.coel_rows.get_row(j).unwrap();
            let corr: f32 = row_i.iter().zip(row_j.iter()).map(|(a, b)| a * b).sum();
            let w = scratch.weights[i] * scratch.weights[j];
            wcoel_num += w * corr;
            wcoel_den += w;
        }
    }
    let weighted_coelution = if wcoel_den > 0.0 {
        wcoel_num / wcoel_den
    } else {
        0.0
    };

    // Gradient consistency: unweighted mean of upper-triangle correlations.
    let gradient_consistency = if n_active_grad >= 2 {
        let mut sum_corr = 0.0f32;
        let mut count = 0u32;
        for i in 0..n_active_grad {
            let row_i = scratch.grad_rows.get_row(i).unwrap();
            for j in (i + 1)..n_active_grad {
                let row_j = scratch.grad_rows.get_row(j).unwrap();
                let corr: f32 = row_i.iter().zip(row_j.iter()).map(|(a, b)| a * b).sum();
                sum_corr += corr;
                count += 1;
            }
        }
        if count > 0 {
            sum_corr / count as f32
        } else {
            0.0
        }
    } else {
        0.0
    };

    // Combined: (1 + 10 * max(wcoel, 0)) * (1 + max(grad, 0))
    let combined =
        (1.0 + 10.0 * weighted_coelution.max(0.0)) * (1.0 + gradient_consistency.max(0.0));

    CoelutionGradientResult {
        weighted_coelution,
        gradient_consistency,
        combined,
    }
}

/// Split product score from independent cosine and scribe apexes (METHODS.md Section 3.1).
///
/// Each profile finds its own argmax, computes area-uniqueness (hw=5) and
/// coelution-gradient (coelution_hw=20, gradient_hw=10) at that apex.
/// base = cos_AU * cos_CG * scr_AU * scr_CG
pub fn compute_split_product<T: KeyLike>(
    cosine_profile: &[f32],
    scribe_profile: &[f32],
    fragments: &MzMajorIntensityArray<T, f32>,
    expected: &[(T, f32)],
    scratch: &mut CoelutionScratch,
) -> SplitProductScore {
    // Independent argmax on each profile
    let cos_apex = argmax(cosine_profile);
    let scr_apex = argmax(scribe_profile);

    // Area-uniqueness with hw=5
    let cos_au = area_uniqueness(cosine_profile, cos_apex, 5);
    let scr_au = area_uniqueness(scribe_profile, scr_apex, 5);

    // Coelution-gradient at each apex — shares the same scratch.
    let cos_cg = coelution_gradient(fragments, expected, cos_apex, 20, 10, scratch);
    let scr_cg = coelution_gradient(fragments, expected, scr_apex, 20, 10, scratch);

    let base_score = cos_au.au_score * cos_cg.combined * scr_au.au_score * scr_cg.combined;

    SplitProductScore {
        cosine_au: cos_au.au_score,
        cosine_cg: cos_cg.combined,
        scribe_au: scr_au.au_score,
        scribe_cg: scr_cg.combined,
        base_score,
        cosine_weighted_coelution: cos_cg.weighted_coelution,
        cosine_gradient_consistency: cos_cg.gradient_consistency,
        scribe_weighted_coelution: scr_cg.weighted_coelution,
        scribe_gradient_consistency: scr_cg.gradient_consistency,
    }
}

/// Find the joint precursor-fragment apex (METHODS.md Section 3.3).
///
/// joint(t) = C(t) * (0.5 + P(t) / max(P))
/// If max(P) == 0, degrades to joint(t) = C(t) * 0.5 (pure fragment apex).
pub fn find_joint_apex(cosine_profile: &[f32], precursor_trace: &[f32]) -> usize {
    let max_p = precursor_trace.iter().copied().fold(0.0f32, f32::max);

    let mut best_val = f32::NEG_INFINITY;
    let mut best_idx = 0usize;

    let n = cosine_profile.len().min(precursor_trace.len());
    for t in 0..n {
        let p_factor = if max_p > 0.0 {
            0.5 + precursor_trace[t] / max_p
        } else {
            0.5
        };
        let joint = cosine_profile[t] * p_factor;
        if joint > best_val {
            best_val = joint;
            best_idx = t;
        }
    }
    best_idx
}

/// Compute all 11 apex-local features at the joint apex (METHODS.md Section 3.4).
///
/// `fragments` and `precursors` are the raw chromatogram data.
/// `expected` contains both fragment and precursor predicted intensities.
/// `cosine_profile` is C(t) = cos(t)^3 * I(t).
/// `precursor_trace` is the summed precursor intensity trace.
/// `joint_apex` is the cycle index from `find_joint_apex`.
/// `n_cycles` is the total number of cycles in the extraction window.
pub fn compute_apex_features<T: KeyLike + Default>(
    fragments: &MzMajorIntensityArray<T, f32>,
    precursors: &MzMajorIntensityArray<i8, f32>,
    expected: &crate::models::ExpectedIntensities<T>,
    cosine_profile: &[f32],
    precursor_trace: &[f32],
    joint_apex: usize,
    n_cycles: usize,
) -> ApexFeatures {
    let apex = joint_apex;

    // ---- Peak Shape (Section 3.4) ----
    let peak_shape = compute_peak_shape(cosine_profile, apex, 10);

    // ---- Ratio CV (Section 3.4) ----
    let ratio_cv = compute_ratio_cv(fragments, expected.fragment_intensities.as_slice(), apex);

    // ---- Centered Apex (Section 3.4) ----
    let centered_apex = if n_cycles > 0 {
        let half = n_cycles as f32 / 2.0;
        (1.0 - (apex as f32 - half).abs() / half).max(0.0)
    } else {
        0.0
    };

    // ---- Precursor Coelution (Section 3.4) ----
    let precursor_coelution =
        compute_precursor_coelution(fragments, precursor_trace, apex, 10, n_cycles);

    // ---- Fragment Coverage (Section 3.4) ----
    let fragment_coverage = compute_fragment_coverage(fragments, apex);

    // ---- Precursor Apex Match (Section 3.4) ----
    let precursor_apex_match = compute_precursor_apex_match(precursor_trace, apex, n_cycles);

    // ---- XIC Quality (Section 3.4) ----
    let xic_quality = compute_xic_quality(fragments, apex, 8);

    // ---- Fragment Apex Agreement (Section 3.4) ----
    let fragment_apex_agreement = compute_fragment_apex_agreement(fragments, apex);

    // ---- Isotope Correlation (Section 3.4) ----
    let isotope_correlation =
        compute_isotope_correlation(precursors, expected.precursor_intensities.as_slice(), apex);

    // ---- Gaussian Correlation (Section 3.4) ----
    let gaussian_correlation = compute_gaussian_correlation(cosine_profile, apex, 15);

    // ---- Per-Fragment Gaussian Correlation (Section 3.4) ----
    let per_frag_gaussian_corr = compute_per_frag_gaussian_corr(fragments, apex, 10);

    ApexFeatures {
        peak_shape,
        ratio_cv,
        centered_apex,
        precursor_coelution,
        fragment_coverage,
        precursor_apex_match,
        xic_quality,
        fragment_apex_agreement,
        isotope_correlation,
        gaussian_correlation,
        per_frag_gaussian_corr,
    }
}

/// Compute the final weighted score (METHODS.md Section 3.5).
///
/// score = base * product(offset_k + scale_k * feature_k)
pub fn compute_weighted_score(base: f32, features: &ApexFeatures) -> f32 {
    let feature_values = [
        features.peak_shape,
        features.ratio_cv,
        features.centered_apex,
        features.precursor_coelution,
        features.fragment_coverage,
        features.precursor_apex_match,
        features.xic_quality,
        features.fragment_apex_agreement,
        features.isotope_correlation,
        features.gaussian_correlation,
        features.per_frag_gaussian_corr,
    ];

    let mut score = base;
    for (i, &fval) in feature_values.iter().enumerate() {
        let (offset, scale) = SCORING_WEIGHTS[i];
        score *= offset + scale * fval;
    }
    score
}

// ---------------------------------------------------------------------------
// Feature implementation helpers (private)
// ---------------------------------------------------------------------------

fn argmax(slice: &[f32]) -> usize {
    if slice.is_empty() {
        return 0;
    }
    let mut best = f32::NEG_INFINITY;
    let mut idx = 0;
    for (i, &v) in slice.iter().enumerate() {
        if v > best {
            best = v;
            idx = i;
        }
    }
    idx
}

/// Peak shape = 0.5 * symmetry + 0.5 * sharpness (Section 3.4).
fn compute_peak_shape(profile: &[f32], apex: usize, hw: usize) -> f32 {
    let n = profile.len();
    if n == 0 {
        return 0.5;
    }

    let lo = apex.saturating_sub(hw);
    let hi = (apex + hw + 1).min(n);

    let left_len = apex - lo;
    let right_len = hi - 1 - apex;
    let flank_len = left_len.min(right_len);

    if flank_len < 2 {
        return 0.5;
    }

    // Symmetry: Pearson correlation between left flank (reversed) and right flank
    let left_rev: Vec<f32> = (0..flank_len).map(|i| profile[apex - 1 - i]).collect();
    let right: Vec<f32> = (0..flank_len).map(|i| profile[apex + 1 + i]).collect();
    let symmetry = pearson_correlation(&left_rev, &right).clamp(0.0, 1.0);

    // Sharpness: 1 - mean(edges) / peak
    let peak_val = profile[apex];
    let sharpness = if peak_val > 0.0 {
        let edge_mean = (profile[lo] + profile[hi - 1]) * 0.5;
        (1.0 - edge_mean / peak_val).max(0.0)
    } else {
        0.0
    };

    0.5 * symmetry + 0.5 * sharpness
}

/// Ratio CV: consistency of observed-to-predicted ratios at apex (Section 3.4).
fn compute_ratio_cv<T: KeyLike>(
    fragments: &MzMajorIntensityArray<T, f32>,
    expected: &[(T, f32)],
    apex: usize,
) -> f32 {
    let n_cycles = fragments.num_cycles();
    if apex >= n_cycles {
        return 0.0;
    }

    let mut ratios = Vec::new();
    for ((key, _mz), chrom) in fragments.iter_mzs() {
        let y_hat = linear_get(expected, key).unwrap_or(0.0);
        if y_hat > 0.0 && apex < chrom.len() && chrom[apex] > 0.0 {
            ratios.push(chrom[apex] / y_hat);
        }
    }

    if ratios.len() < 3 {
        return 0.0;
    }

    let mean: f32 = ratios.iter().sum::<f32>() / ratios.len() as f32;
    if mean <= 0.0 {
        return 0.0;
    }
    let var: f32 = ratios.iter().map(|r| (r - mean).powi(2)).sum::<f32>() / ratios.len() as f32;
    let cv = var.sqrt() / mean;
    1.0 / (1.0 + cv)
}

/// Precursor coelution: Pearson correlation between precursor trace and summed
/// fragment trace in a +/-hw window around apex (Section 3.4).
fn compute_precursor_coelution<T: KeyLike>(
    fragments: &MzMajorIntensityArray<T, f32>,
    precursor_trace: &[f32],
    apex: usize,
    hw: usize,
    n_cycles: usize,
) -> f32 {
    let lo = apex.saturating_sub(hw);
    let hi = (apex + hw + 1).min(n_cycles);
    if hi <= lo || hi > precursor_trace.len() {
        return 0.0;
    }

    // Sum fragment traces in the window
    let win_len = hi - lo;
    let mut frag_sum = vec![0.0f32; win_len];
    for ((_key, _mz), chrom) in fragments.iter_mzs() {
        for i in 0..win_len {
            let idx = lo + i;
            if idx < chrom.len() {
                frag_sum[i] += chrom[idx];
            }
        }
    }

    let prec_win = &precursor_trace[lo..hi];
    pearson_correlation(prec_win, &frag_sum).max(0.0)
}

/// Fragment coverage: fraction of fragments with nonzero intensity at apex (Section 3.4).
fn compute_fragment_coverage<T: KeyLike>(
    fragments: &MzMajorIntensityArray<T, f32>,
    apex: usize,
) -> f32 {
    let n_frags = fragments.num_ions();
    if n_frags == 0 {
        return 0.0;
    }
    let mut count = 0u32;
    for ((_key, _mz), chrom) in fragments.iter_mzs() {
        if apex < chrom.len() && chrom[apex] > 0.0 {
            count += 1;
        }
    }
    count as f32 / n_frags as f32
}

/// Precursor apex match = 0.5 * proximity + 0.5 * fraction (Section 3.4).
fn compute_precursor_apex_match(precursor_trace: &[f32], apex: usize, n_cycles: usize) -> f32 {
    if precursor_trace.is_empty() || n_cycles == 0 {
        return 0.0;
    }

    // Proximity: max(0, 1 - |t*_P - t*| / (T/4))
    let prec_apex = argmax(precursor_trace);
    let quarter = n_cycles as f32 / 4.0;
    let proximity = if quarter > 0.0 {
        (1.0 - (prec_apex as f32 - apex as f32).abs() / quarter).max(0.0)
    } else {
        0.0
    };

    // Fraction: P(t*) / sum(P)
    let total_p: f32 = precursor_trace.iter().sum();
    let fraction = if total_p > 0.0 && apex < precursor_trace.len() {
        precursor_trace[apex] / total_p
    } else {
        0.0
    };

    0.5 * proximity + 0.5 * fraction
}

/// XIC quality: mean per-fragment chromatographic peak quality (Section 3.4).
/// For each fragment in a +/-hw window:
///   alignment = max(0, 1 - d_i / (w/2))
///   sharpness = 1 - mean(edges) / peak
///   xic_i = 0.5 * alignment + 0.5 * sharpness
fn compute_xic_quality<T: KeyLike>(
    fragments: &MzMajorIntensityArray<T, f32>,
    apex: usize,
    hw: usize,
) -> f32 {
    let n_cycles = fragments.num_cycles();
    let n_frags = fragments.num_ions();
    if n_frags == 0 || n_cycles == 0 {
        return 0.0;
    }

    let lo = apex.saturating_sub(hw);
    let hi = (apex + hw + 1).min(n_cycles);
    let win_len = hi - lo;

    if win_len < 2 {
        return 0.0;
    }

    let half_w = win_len as f32 / 2.0;
    let mut sum_quality = 0.0f32;

    for ((_key, _mz), chrom) in fragments.iter_mzs() {
        let window = &chrom[lo..hi];

        // Find local max in window
        let (local_max_val, local_max_idx) =
            window
                .iter()
                .enumerate()
                .fold((0.0f32, 0usize), |(best_v, best_i), (i, &v)| {
                    if v > best_v { (v, i) } else { (best_v, best_i) }
                });

        if local_max_val <= 0.0 {
            // Fragment contributes 0
            continue;
        }

        // Alignment: distance from local max to expected apex position in window
        let apex_in_window = apex - lo;
        let d = (local_max_idx as f32 - apex_in_window as f32).abs();
        let alignment = (1.0 - d / half_w).max(0.0);

        // Sharpness: 1 - mean(edges) / peak
        let edge_mean = (window[0] + window[win_len - 1]) * 0.5;
        let sharpness = (1.0 - edge_mean / local_max_val).max(0.0);

        sum_quality += 0.5 * alignment + 0.5 * sharpness;
    }

    sum_quality / n_frags as f32
}

/// Fragment apex agreement: fraction of fragments whose argmax is within +/-2
/// of the joint apex (Section 3.4).
fn compute_fragment_apex_agreement<T: KeyLike>(
    fragments: &MzMajorIntensityArray<T, f32>,
    apex: usize,
) -> f32 {
    let n_frags = fragments.num_ions();
    if n_frags == 0 {
        return 0.0;
    }

    let mut count = 0u32;
    for ((_key, _mz), chrom) in fragments.iter_mzs() {
        // Single pass: find argmax and check if max > 0
        let (frag_apex, max_val) = chrom.iter().enumerate().fold(
            (0usize, 0.0f32),
            |(bi, bv), (i, &v)| if v > bv { (i, v) } else { (bi, bv) },
        );
        if max_val <= 0.0 {
            continue;
        }
        if (frag_apex as i64 - apex as i64).unsigned_abs() <= 2 {
            count += 1;
        }
    }

    count as f32 / n_frags as f32
}

/// Isotope correlation: cosine between observed and expected precursor isotope
/// envelope at the apex cycle (Section 3.4).
///
/// Uses existing predicted precursor_intensities (keys 0, 1, 2) rather than
/// re-implementing averagine — the codebase already has sequence-specific predictions.
fn compute_isotope_correlation(
    precursors: &MzMajorIntensityArray<i8, f32>,
    expected_precursor: &[(i8, f32)],
    apex: usize,
) -> f32 {
    let n_cycles = precursors.num_cycles();
    if apex >= n_cycles {
        return 0.0;
    }

    // Collect observed and expected for isotope keys 0, 1, 2
    let mut obs = [0.0f32; 3];
    let mut exp = [0.0f32; 3];
    let mut n_valid = 0u32;

    for iso_key in 0i8..=2i8 {
        if let Some(row) = precursors.get_row(&iso_key) {
            if apex < row.len() {
                obs[iso_key as usize] = row[apex];
            }
        }
        if let Some(v) = linear_get(expected_precursor, &iso_key) {
            exp[iso_key as usize] = v;
        }
        if obs[iso_key as usize] > 0.0 && exp[iso_key as usize] > 0.0 {
            n_valid += 1;
        }
    }

    if n_valid < 2 {
        return 0.0;
    }

    // Cosine similarity
    let dot: f32 = obs.iter().zip(exp.iter()).map(|(a, b)| a * b).sum();
    let norm_obs = obs.iter().map(|v| v * v).sum::<f32>().sqrt();
    let norm_exp = exp.iter().map(|v| v * v).sum::<f32>().sqrt();

    if norm_obs < 1e-12 || norm_exp < 1e-12 {
        return 0.0;
    }
    (dot / (norm_obs * norm_exp)).clamp(0.0, 1.0)
}

/// Build a Gaussian reference vector centered at `apex` within window `[lo, lo+win_len)`.
fn build_gaussian_reference(win_len: usize, lo: usize, apex: usize, sigma: f32) -> Vec<f32> {
    (0..win_len)
        .map(|i| {
            let t = (lo + i) as f32 - apex as f32;
            (-t * t / (2.0 * sigma * sigma)).exp()
        })
        .collect()
}

/// Gaussian correlation: Pearson correlation between the combined elution profile
/// and an ideal Gaussian centered at the apex (Section 3.4).
///
/// sigma is estimated from the observed FWHM: sigma = FWHM / 2.355
fn compute_gaussian_correlation(profile: &[f32], apex: usize, hw: usize) -> f32 {
    let n = profile.len();
    if n == 0 {
        return 0.0;
    }

    let lo = apex.saturating_sub(hw);
    let hi = (apex + hw + 1).min(n);
    let win_len = hi - lo;
    if win_len < 3 {
        return 0.0;
    }

    let fwhm = estimate_fwhm(profile, apex);
    let sigma = (fwhm / 2.355).max(0.5);
    let gaussian = build_gaussian_reference(win_len, lo, apex, sigma);

    let window = &profile[lo..hi];
    pearson_correlation(window, &gaussian).max(0.0)
}

/// Per-fragment Gaussian correlation: mean per-fragment correlation with a Gaussian
/// reference (Section 3.4, inspired by Beta-DIA).
///
/// sigma = max(window_size / 6, 1)
fn compute_per_frag_gaussian_corr<T: KeyLike>(
    fragments: &MzMajorIntensityArray<T, f32>,
    apex: usize,
    hw: usize,
) -> f32 {
    let n_cycles = fragments.num_cycles();
    if n_cycles == 0 {
        return 0.0;
    }

    let lo = apex.saturating_sub(hw);
    let hi = (apex + hw + 1).min(n_cycles);
    let win_len = hi - lo;
    if win_len < 3 {
        return 0.0;
    }

    let sigma = (win_len as f32 / 6.0).max(1.0);
    let gaussian = build_gaussian_reference(win_len, lo, apex, sigma);

    let mut sum_corr = 0.0f32;
    let mut n_active = 0u32;

    for ((_key, _mz), chrom) in fragments.iter_mzs() {
        let max_val = chrom.iter().copied().fold(0.0f32, f32::max);
        if max_val <= 0.0 {
            continue;
        }
        let window = &chrom[lo..hi];
        let corr = pearson_correlation(window, &gaussian).max(0.0);
        sum_corr += corr;
        n_active += 1;
    }

    if n_active == 0 {
        return 0.0;
    }
    sum_corr / n_active as f32
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use timsquery::models::Array2D;

    /// Helper: build a MzMajorIntensityArray from a Vec of (key, mz, intensities) tuples.
    fn make_fragments<T: KeyLike>(ions: Vec<(T, f64, Vec<f32>)>) -> MzMajorIntensityArray<T, f32> {
        let n_cycles = ions.first().map(|(_, _, v)| v.len()).unwrap_or(0);
        let n_ions = ions.len();
        let mz_order: Vec<(T, f64)> = ions.iter().map(|(k, mz, _)| (k.clone(), *mz)).collect();
        let flat: Vec<f32> = ions.into_iter().flat_map(|(_, _, v)| v).collect();
        let arr = Array2D::from_flat_vector(flat, n_ions, n_cycles).unwrap();
        MzMajorIntensityArray {
            arr,
            mz_order,
            cycle_offset: 0,
        }
    }

    // --- Test 3: area_uniqueness ---

    #[test]
    fn test_area_uniqueness_clear_peak() {
        // A clear peak at index 10, surrounded by near-zero
        let mut signal = vec![0.0f32; 21];
        signal[8] = 1.0;
        signal[9] = 5.0;
        signal[10] = 10.0;
        signal[11] = 5.0;
        signal[12] = 1.0;

        let result = area_uniqueness(&signal, 10, 5);

        // peak_area = sum of signal[5..=15] = 0+0+0+1+5+10+5+1+0+0+0 = 22
        // total = 22
        // uniqueness = 22/22 = 1.0
        // au_score = 22 * (1 + 200*1.0) = 22 * 201 = 4422
        assert!((result.au_score - 4422.0).abs() < 1.0);
    }

    #[test]
    fn test_area_uniqueness_all_zeros() {
        let signal = vec![0.0f32; 20];
        let result = area_uniqueness(&signal, 10, 5);
        assert_eq!(result.au_score, 0.0);
    }

    #[test]
    fn test_area_uniqueness_boundary_clamping() {
        // Apex near left edge
        let mut signal = vec![0.0f32; 10];
        signal[0] = 5.0;
        signal[1] = 10.0;
        signal[2] = 5.0;

        let result = area_uniqueness(&signal, 1, 5);
        // Window: [0..min(7,10)] = [0..7], peak_area = 5+10+5+0+0+0+0 = 20
        // total = 20
        assert!(result.au_score > 0.0);

        // Apex near right edge
        let mut signal2 = vec![0.0f32; 10];
        signal2[8] = 10.0;
        signal2[9] = 5.0;

        let result2 = area_uniqueness(&signal2, 9, 5);
        assert!(result2.au_score > 0.0);
    }

    // --- Test 4: coelution_gradient ---

    #[test]
    fn test_coelution_gradient_identical_fragments() {
        // Two fragments with identical chromatograms => perfect coelution (wcoel = 1)
        let peak: Vec<f32> = (0..41)
            .map(|i| {
                let x = (i as f32 - 20.0) / 5.0;
                (-x * x / 2.0).exp() * 100.0
            })
            .collect();

        let fragments = make_fragments(vec![
            ("a".to_string(), 100.0, peak.clone()),
            ("b".to_string(), 200.0, peak.clone()),
        ]);

        let expected: Vec<(String, f32)> = vec![("a".to_string(), 1.0), ("b".to_string(), 1.0)];

        let mut scratch = CoelutionScratch::new();
        let result = coelution_gradient(&fragments, &expected, 20, 20, 10, &mut scratch);

        // Identical traces => correlation should be ~1.0
        assert!(
            result.weighted_coelution > 0.95,
            "wcoel = {} should be ~1.0",
            result.weighted_coelution
        );
        assert!(result.combined > 1.0);
    }

    #[test]
    fn test_coelution_gradient_anti_correlated() {
        // Two fragments with anti-correlated chromatograms
        let n = 41;
        let peak_a: Vec<f32> = (0..n)
            .map(|i| {
                let x = (i as f32 - 10.0) / 3.0;
                (-x * x / 2.0).exp() * 100.0
            })
            .collect();
        let peak_b: Vec<f32> = (0..n)
            .map(|i| {
                let x = (i as f32 - 30.0) / 3.0;
                (-x * x / 2.0).exp() * 100.0
            })
            .collect();

        let fragments = make_fragments(vec![
            ("a".to_string(), 100.0, peak_a),
            ("b".to_string(), 200.0, peak_b),
        ]);

        let expected: Vec<(String, f32)> = vec![("a".to_string(), 1.0), ("b".to_string(), 1.0)];

        let mut scratch = CoelutionScratch::new();
        let result = coelution_gradient(&fragments, &expected, 20, 20, 10, &mut scratch);
        // Anti-correlated: wcoel should be negative or near zero
        assert!(
            result.weighted_coelution < 0.5,
            "wcoel = {} should be low for anti-correlated traces",
            result.weighted_coelution
        );
    }

    #[test]
    fn test_coelution_gradient_single_fragment() {
        // Single fragment: fewer than 2 active => combined = 1.0
        let peak: Vec<f32> = (0..41)
            .map(|i| {
                let x = (i as f32 - 20.0) / 5.0;
                (-x * x / 2.0).exp() * 100.0
            })
            .collect();

        let fragments = make_fragments(vec![("a".to_string(), 100.0, peak)]);

        let expected: Vec<(String, f32)> = vec![("a".to_string(), 1.0)];

        let mut scratch = CoelutionScratch::new();
        let result = coelution_gradient(&fragments, &expected, 20, 20, 10, &mut scratch);
        assert!(
            (result.combined - 1.0).abs() < 1e-6,
            "single fragment => combined should be 1.0, got {}",
            result.combined
        );
    }

    // --- Test 5: split product with offset peaks ---

    #[test]
    fn test_split_product_offset_peaks() {
        let n = 50;
        // Cosine profile peaks at 20, scribe at 30
        let cosine_profile: Vec<f32> = (0..n)
            .map(|i| {
                let x = (i as f32 - 20.0) / 3.0;
                (-x * x / 2.0).exp() * 10.0
            })
            .collect();
        let scribe_profile: Vec<f32> = (0..n)
            .map(|i| {
                let x = (i as f32 - 30.0) / 3.0;
                (-x * x / 2.0).exp() * 8.0
            })
            .collect();

        let peak: Vec<f32> = (0..n)
            .map(|i| {
                let x = (i as f32 - 25.0) / 5.0;
                (-x * x / 2.0).exp() * 100.0
            })
            .collect();

        let fragments = make_fragments(vec![
            ("a".to_string(), 100.0, peak.clone()),
            ("b".to_string(), 200.0, peak),
        ]);
        let expected: Vec<(String, f32)> = vec![("a".to_string(), 1.0), ("b".to_string(), 1.0)];

        let mut scratch = CoelutionScratch::new();
        let result = compute_split_product(
            &cosine_profile,
            &scribe_profile,
            &fragments,
            &expected,
            &mut scratch,
        );
        assert!(result.base_score > 0.0);
        assert!(result.cosine_au > 0.0);
        assert!(result.scribe_au > 0.0);
    }

    // --- Test 6: joint apex ---

    #[test]
    fn test_joint_apex_with_precursor() {
        let n = 30;
        // Fragment profile peaks at 15
        let cosine_profile: Vec<f32> = (0..n)
            .map(|i| {
                let x = (i as f32 - 15.0) / 3.0;
                (-x * x / 2.0).exp() * 10.0
            })
            .collect();
        // Precursor trace also peaks at 15
        let precursor_trace: Vec<f32> = (0..n)
            .map(|i| {
                let x = (i as f32 - 15.0) / 4.0;
                (-x * x / 2.0).exp() * 5.0
            })
            .collect();

        let apex = find_joint_apex(&cosine_profile, &precursor_trace);
        assert_eq!(apex, 15);
    }

    #[test]
    fn test_joint_apex_without_precursor() {
        let n = 30;
        let cosine_profile: Vec<f32> = (0..n)
            .map(|i| {
                let x = (i as f32 - 15.0) / 3.0;
                (-x * x / 2.0).exp() * 10.0
            })
            .collect();
        // Zero precursor
        let precursor_trace = vec![0.0f32; n];

        let apex = find_joint_apex(&cosine_profile, &precursor_trace);
        // Should still find the cosine peak at 15
        assert_eq!(apex, 15);
    }

    // --- Test 7: individual feature tests ---

    #[test]
    fn test_fragment_coverage() {
        // 4 fragments, 2 have signal at apex
        let fragments = make_fragments(vec![
            ("a".to_string(), 100.0, vec![0.0, 5.0, 0.0]),
            ("b".to_string(), 200.0, vec![0.0, 0.0, 0.0]),
            ("c".to_string(), 300.0, vec![0.0, 3.0, 0.0]),
            ("d".to_string(), 400.0, vec![0.0, 0.0, 0.0]),
        ]);
        let cov = compute_fragment_coverage(&fragments, 1);
        assert!((cov - 0.5).abs() < 1e-6, "expected 0.5, got {}", cov);
    }

    #[test]
    fn test_ratio_cv_uniform_ratios() {
        // All fragments have the same obs/pred ratio => CV = 0 => ratio_cv = 1.0
        let fragments = make_fragments(vec![
            ("a".to_string(), 100.0, vec![0.0, 2.0, 0.0]),
            ("b".to_string(), 200.0, vec![0.0, 4.0, 0.0]),
            ("c".to_string(), 300.0, vec![0.0, 6.0, 0.0]),
            ("d".to_string(), 400.0, vec![0.0, 8.0, 0.0]),
        ]);
        let expected: Vec<(String, f32)> = vec![
            ("a".to_string(), 1.0),
            ("b".to_string(), 2.0),
            ("c".to_string(), 3.0),
            ("d".to_string(), 4.0),
        ];
        let cv = compute_ratio_cv(&fragments, &expected, 1);
        assert!(
            (cv - 1.0).abs() < 1e-6,
            "uniform ratios => cv should be 1.0, got {}",
            cv
        );
    }

    #[test]
    fn test_isotope_correlation_perfect_match() {
        // Observed matches expected perfectly
        let precursors = make_fragments(vec![
            (0i8, 500.0, vec![0.0, 0.6, 0.0]),
            (1i8, 500.5, vec![0.0, 0.3, 0.0]),
            (2i8, 501.0, vec![0.0, 0.1, 0.0]),
        ]);
        let expected: Vec<(i8, f32)> = vec![(0i8, 0.6), (1, 0.3), (2, 0.1)];
        let iso = compute_isotope_correlation(&precursors, &expected, 1);
        assert!(
            (iso - 1.0).abs() < 1e-4,
            "perfect match => iso should be ~1.0, got {}",
            iso
        );
    }

    // --- Test 8: weighted score ---

    #[test]
    fn test_weighted_score_all_zero_features() {
        let features = ApexFeatures {
            peak_shape: 0.0,
            ratio_cv: 0.0,
            centered_apex: 0.0,
            precursor_coelution: 0.0,
            fragment_coverage: 0.0,
            precursor_apex_match: 0.0,
            xic_quality: 0.0,
            fragment_apex_agreement: 0.0,
            isotope_correlation: 0.0,
            gaussian_correlation: 0.0,
            per_frag_gaussian_corr: 0.0,
        };
        let score = compute_weighted_score(100.0, &features);
        // product of offsets: 1*1*1*1*0.32*1*1*0.43*1*0.27*0.65
        let expected = 100.0 * 1.0 * 1.0 * 1.0 * 1.0 * 0.32 * 1.0 * 1.0 * 0.43 * 1.0 * 0.27 * 0.65;
        assert!(
            (score - expected).abs() < 1e-3,
            "score {} != expected {}",
            score,
            expected
        );
    }

    #[test]
    fn test_weighted_score_all_one_features() {
        let features = ApexFeatures {
            peak_shape: 1.0,
            ratio_cv: 1.0,
            centered_apex: 1.0,
            precursor_coelution: 1.0,
            fragment_coverage: 1.0,
            precursor_apex_match: 1.0,
            xic_quality: 1.0,
            fragment_apex_agreement: 1.0,
            isotope_correlation: 1.0,
            gaussian_correlation: 1.0,
            per_frag_gaussian_corr: 1.0,
        };
        let score = compute_weighted_score(1.0, &features);
        // product of (offset + scale): 4.5 * 4.0 * 5.5 * 3.3 * 1.32 * 9.0 * 8.8 * 1.43 * 3.6 * 1.27 * 1.65
        let expected = 4.5 * 4.0 * 5.5 * 3.3 * 1.32 * 9.0 * 8.8 * 1.43 * 3.6 * 1.27 * 1.65;
        assert!(
            (score - expected).abs() / expected < 1e-3,
            "score {} != expected {}",
            score,
            expected
        );
    }
}
