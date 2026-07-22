//! Reusable numeric primitives for apex-local scoring. These carry no scoring
//! semantics — they are the shared math toolkit that the apex-stage score
//! families ([`blocks::apex_features`], [`blocks::split_product`]) build on.
//! Score-specific building blocks (area-uniqueness, coelution-gradient, the
//! weighted-score weights) live in their own family file, not here.

/// Sentinel Scribe score used when no signal is observed (fills the scribe
/// trace so downstream argmax/area math sees a defined floor).
pub const SCRIBE_FLOOR: f32 = -100.0;

/// Index of the maximum element (0 for an empty slice).
pub fn argmax(slice: &[f32]) -> usize {
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

/// Pearson correlation between two equal-length slices.
/// Returns 0.0 if either has zero variance.
pub fn pearson_correlation(a: &[f32], b: &[f32]) -> f32 {
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
pub fn center_normalize(window: &[f32], out: &mut [f32]) -> bool {
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
pub fn estimate_fwhm(profile: &[f32], apex: usize) -> f32 {
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

    (left_dist + right_dist).max(1.0)
}

/// Build a Gaussian reference vector centered at `apex` within window `[lo, lo+win_len)`.
pub fn build_gaussian_reference(win_len: usize, lo: usize, apex: usize, sigma: f32) -> Vec<f32> {
    (0..win_len)
        .map(|i| {
            let t = (lo + i) as f32 - apex as f32;
            (-t * t / (2.0 * sigma * sigma)).exp()
        })
        .collect()
}
