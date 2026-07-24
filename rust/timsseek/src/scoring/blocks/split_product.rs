//! Split-product family: cosine/scribe apex agreement scores.
//!
//! Owns its whole lifecycle in this file: the derive-generated projection
//! (`#[derive(ScoreBlock)]`), the raw compute (`compute_split_product`, run at
//! the apex stage while chromatogram buffers are live) and its building blocks
//! (`area_uniqueness`, `coelution_gradient`) plus the reusable
//! [`CoelutionScratch`] buffers the scorer owns for allocation reuse. Generic
//! numeric primitives live in `crate::scoring::apex_dsp`.

use array2d::Array2D;
use timsquery::models::MzMajorIntensityArray;
use timsquery::traits::KeyLike;
use timsseek_macros::ScoreBlock;

use crate::models::query_item::linear_get;
use crate::scoring::apex_dsp::{
    argmax,
    center_normalize,
};

/// Stage: apex (built from `SplitProductScore`, computed while chromatogram
/// buffers are live).
#[derive(Debug, Clone, Copy, ::serde::Serialize, ScoreBlock)]
pub struct SplitProduct {
    #[feat(ln1p)]
    pub split_product_score: f32,
    #[feat(ln1p)]
    pub cosine_au: f32,
    #[feat(ln1p)]
    pub scribe_au: f32,
    #[feat(raw)]
    pub cosine_cg: f32,
    #[feat(raw)]
    pub scribe_cg: f32,
    #[feat(raw)]
    pub cosine_weighted_coelution: f32,
    #[feat(raw)]
    pub cosine_gradient_consistency: f32,
    #[feat(raw)]
    pub scribe_weighted_coelution: f32,
    #[feat(raw)]
    pub scribe_gradient_consistency: f32,
}

impl From<&SplitProductScore> for SplitProduct {
    fn from(s: &SplitProductScore) -> Self {
        Self {
            split_product_score: s.base_score,
            cosine_au: s.cosine_au,
            scribe_au: s.scribe_au,
            cosine_cg: s.cosine_cg,
            scribe_cg: s.scribe_cg,
            cosine_weighted_coelution: s.cosine_weighted_coelution,
            cosine_gradient_consistency: s.cosine_gradient_consistency,
            scribe_weighted_coelution: s.scribe_weighted_coelution,
            scribe_gradient_consistency: s.scribe_gradient_consistency,
        }
    }
}

/// Split product score computed from independent cosine and scribe apexes.
/// Raw output of `compute_split_product`, projected
/// into [`SplitProduct`] via [`From`].
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

/// Reusable scratch buffers for `coelution_gradient`. Owned by the scorer so
/// allocations happen once per rayon worker split and are reused across every
/// peptide that worker processes. See [`CoelutionScratch::with_frag_capacity`]
/// for sizing.
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

    /// Reserve scratch capacity for a library whose peptides have at most
    /// `frag_capacity` fragments. Coel window = 2*20+1 = 41 cycles,
    /// gradient diffs = 2*10 = 20. If callers pass the library-wide max
    /// (pre-scanned via `items_to_score.iter().map(|i|
    /// i.expected_intensity.fragment_len()).max()`), no realloc occurs for
    /// any peptide.
    pub fn with_frag_capacity(frag_capacity: usize) -> Self {
        Self {
            coel_rows: Array2D::with_capacity(Self::TYP_COEL_LEN, frag_capacity),
            grad_rows: Array2D::with_capacity(Self::TYP_DIFF_LEN, frag_capacity),
            weights: Vec::with_capacity(frag_capacity),
            diffs: Vec::with_capacity(Self::TYP_DIFF_LEN),
        }
    }
}

/// Result of the area-uniqueness calculation.
#[derive(Debug, Clone, Copy)]
struct AreaUniquenessResult {
    pub au_score: f32,
}

/// Result of the coelution-gradient quality calculation.
#[derive(Debug, Clone, Copy)]
struct CoelutionGradientResult {
    pub weighted_coelution: f32,
    pub gradient_consistency: f32,
    pub combined: f32,
}

/// Split product score from independent cosine and scribe apexes.
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

/// Area-uniqueness score around an apex.
///
/// `signal` is the per-cycle profile (e.g. cosine_profile or scribe_profile).
/// `apex` is the argmax index within `signal`.
/// `hw` is the half-width of the peak window (typically 5).
///
/// AU = peak_area * (1 + 200 * peak_area / total_area)
fn area_uniqueness(signal: &[f32], apex: usize, hw: usize) -> AreaUniquenessResult {
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

/// Coelution-gradient quality at a given apex.
///
/// Measures whether the fragment XICs co-elute and change in the same direction.
/// Weights use raw predicted intensity fractions (NOT sqrt-transformed).
///
/// Returns `combined = 1.0` if fewer than 2 active fragments in the window.
fn coelution_gradient<T: KeyLike>(
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

        let mut scratch = CoelutionScratch::with_frag_capacity(16);
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

        let mut scratch = CoelutionScratch::with_frag_capacity(16);
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

        let mut scratch = CoelutionScratch::with_frag_capacity(16);
        let result = coelution_gradient(&fragments, &expected, 20, 20, 10, &mut scratch);
        assert!(
            (result.combined - 1.0).abs() < 1e-6,
            "single fragment => combined should be 1.0, got {}",
            result.combined
        );
    }

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

        let mut scratch = CoelutionScratch::with_frag_capacity(16);
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
}
