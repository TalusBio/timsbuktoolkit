//! Apex-features family: the 11 apex-local features.
//!
//! Owns its whole lifecycle in this file: the macro-generated struct +
//! projection ([`crate::score_block!`], all 11 fields `#[raw]`), the compute
//! (`compute_apex_features`, run at the apex stage while the chromatogram
//! buffers are live), and the final weighted score (`compute_weighted_score`).
//! Reusable numeric primitives it leans on live in `crate::scoring::apex_dsp`.

use timsquery::models::MzMajorIntensityArray;
use timsquery::traits::KeyLike;

use crate::models::ExpectedIntensities;
use crate::models::query_item::linear_get;
use crate::score_block;
use crate::scoring::apex_dsp::{
    argmax,
    build_gaussian_reference,
    estimate_fwhm,
    pearson_correlation,
};

score_block! {
    /// The 11 apex-local features.
    /// Stage: apex (computed while chromatogram buffers are live).
    pub struct ApexFeatures {
        #[raw] pub peak_shape: f32,
        #[raw] pub ratio_cv: f32,
        #[raw] pub centered_apex: f32,
        #[raw] pub precursor_coelution: f32,
        #[raw] pub fragment_coverage: f32,
        #[raw] pub precursor_apex_match: f32,
        #[raw] pub xic_quality: f32,
        #[raw] pub fragment_apex_agreement: f32,
        #[raw] pub isotope_correlation: f32,
        #[raw] pub gaussian_correlation: f32,
        #[raw] pub per_frag_gaussian_corr: f32,
    }
}

/// Compute all 11 apex-local features at the apex.
///
/// `fragments` and `precursors` are the raw chromatogram data.
/// `expected` contains both fragment and precursor predicted intensities.
/// `cosine_profile` is C(t) = cos(t)^3 * I(t).
/// `precursor_trace` is the summed precursor intensity trace.
/// `joint_apex` is the apex cycle index (Pass 1's weighted-apex_profile pick).
/// `n_cycles` is the total number of cycles in the extraction window.
pub fn compute_apex_features<T: KeyLike + Default>(
    fragments: &MzMajorIntensityArray<T, f32>,
    precursors: &MzMajorIntensityArray<i8, f32>,
    expected: &ExpectedIntensities<T>,
    cosine_profile: &[f32],
    precursor_trace: &[f32],
    joint_apex: usize,
    n_cycles: usize,
) -> ApexFeatures {
    let apex = joint_apex;

    // ---- Peak Shape ----
    let peak_shape = compute_peak_shape(cosine_profile, apex, 10);

    // ---- Ratio CV ----
    let ratio_cv = compute_ratio_cv(fragments, expected.fragment_intensities.as_slice(), apex);

    // ---- Centered Apex ----
    let centered_apex = if n_cycles > 0 {
        let half = n_cycles as f32 / 2.0;
        (1.0 - (apex as f32 - half).abs() / half).max(0.0)
    } else {
        0.0
    };

    // ---- Precursor Coelution ----
    let precursor_coelution =
        compute_precursor_coelution(fragments, precursor_trace, apex, 10, n_cycles);

    // ---- Fragment Coverage ----
    let fragment_coverage = compute_fragment_coverage(fragments, apex);

    // ---- Precursor Apex Match ----
    let precursor_apex_match = compute_precursor_apex_match(precursor_trace, apex, n_cycles);

    // ---- XIC Quality ----
    let xic_quality = compute_xic_quality(fragments, apex, 8);

    // ---- Fragment Apex Agreement ----
    let fragment_apex_agreement = compute_fragment_apex_agreement(fragments, apex);

    // ---- Isotope Correlation ----
    let isotope_correlation =
        compute_isotope_correlation(precursors, expected.precursor_intensities.as_slice(), apex);

    // ---- Gaussian Correlation ----
    let gaussian_correlation = compute_gaussian_correlation(cosine_profile, apex, 15);

    // ---- Per-Fragment Gaussian Correlation ----
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

/// Compute the final weighted score.
///
/// score = base * product(offset_k + scale_k * feature_k)
pub fn compute_weighted_score(base: f32, features: &ApexFeatures) -> f32 {
    // (feature value, offset, scale) kept as one list so a value can never drift
    // out of alignment with its weight. Final score = base * product(offset +
    // scale * feature_k).
    let terms = [
        (features.peak_shape, 1.0, 3.5),
        (features.ratio_cv, 1.0, 3.0),
        (features.centered_apex, 1.0, 4.5),
        (features.precursor_coelution, 1.0, 2.3),
        (features.fragment_coverage, 0.32, 1.0),
        (features.precursor_apex_match, 1.0, 8.0),
        (features.xic_quality, 1.0, 7.8),
        (features.fragment_apex_agreement, 0.43, 1.0),
        (features.isotope_correlation, 1.0, 2.6),
        (features.gaussian_correlation, 0.27, 1.0),
        (features.per_frag_gaussian_corr, 0.65, 1.0),
    ];

    terms.iter().fold(base, |score, &(fval, offset, scale)| {
        score * (offset + scale * fval)
    })
}

// ---------------------------------------------------------------------------
// Feature implementation helpers (private)
// ---------------------------------------------------------------------------

/// Peak shape = 0.5 * symmetry + 0.5 * sharpness.
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

/// Ratio CV: consistency of observed-to-predicted ratios at apex.
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
/// fragment trace in a +/-hw window around apex.
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
        let avail = chrom.len().saturating_sub(lo).min(win_len);
        let chrom_win = &chrom[lo..lo + avail];
        for (sum, &x) in frag_sum.iter_mut().zip(chrom_win.iter()) {
            *sum += x;
        }
    }

    let prec_win = &precursor_trace[lo..hi];
    pearson_correlation(prec_win, &frag_sum).max(0.0)
}

/// Fragment coverage: fraction of fragments with nonzero intensity at apex.
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

/// Precursor apex match = 0.5 * proximity + 0.5 * fraction.
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

/// XIC quality: mean per-fragment chromatographic peak quality.
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
/// of the joint apex.
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
/// envelope at the apex cycle.
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
        if let Some(row) = precursors.get_row(&iso_key)
            && apex < row.len()
        {
            obs[iso_key as usize] = row[apex];
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

/// Gaussian correlation: Pearson correlation between the combined elution profile
/// and an ideal Gaussian centered at the apex.
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
/// reference (inspired by Beta-DIA).
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
