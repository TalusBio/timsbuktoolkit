//! Optional per-cycle weighting passes applied to the composite apex profile
//! before peak-picking.
//!
//! The base profile (`cos^p * I^q * (s_ratio + s_norm)`, optionally blurred)
//! is built in [`crate::scoring::apex_finding`]. The passes here multiply that
//! profile in place: each rewards genuine multi-signal coincidence at a cycle
//! and suppresses random-interferent pileups, without depending on absolute
//! fragment ratios. All knobs and their shipping defaults are bench-validated
//! on the `apex_sim` canonical suite (see [`ApexConfig`]).

use crate::scoring::apex_finding::Extraction;
use timsquery::traits::KeyLike;

mod coelution;
mod vote;

/// Tunable apex-profile knobs. `Default` holds the bench-validated shipping
/// values. Ablation flips individual fields; shipped values live in `Default`.
#[derive(Debug, Clone)]
pub struct ApexConfig {
    /// Cosine exponent in the base profile: `cos^cos_pow`.
    pub cos_pow: f32,
    /// Log-intensity exponent in the base profile: `I^i_exp`.
    pub i_exp: f32,
    /// Additive scribe floor: `apex = C * (s_ratio + s_norm)`.
    pub s_ratio: f32,
    /// Gaussian-blur passes applied to the base profile (0 = none).
    pub blur_passes: usize,
    /// Joint-apex snap window: use joint apex if within this many cycles.
    pub joint_snap: usize,
    /// Ratio-free coelution weight strength (0 = disabled).
    pub coel_k: f32,
    /// Coelution correlation half-window (cycles).
    pub coel_w: f32,
    /// Cross-row coincidence-vote weight strength (0 = disabled).
    pub vote_k: f32,
    /// Vote sigmoid threshold (z-score units).
    pub vote_tau: f32,
    /// Vote sigmoid softness (z-score units).
    pub vote_s: f32,
}

impl Default for ApexConfig {
    fn default() -> Self {
        Self {
            cos_pow: 0.5,
            i_exp: 0.75,
            s_ratio: 1.75,
            blur_passes: 1,
            joint_snap: 1,
            coel_k: 1.0,
            coel_w: 2.0,
            vote_k: 14.0,
            vote_tau: 1.3,
            vote_s: 0.4,
        }
    }
}

/// Reusable per-worker scratch for the weight passes. Owned by the scorer and
/// reused across candidates so the hot loop never allocates.
#[derive(Debug, Default)]
pub struct ApexScratch {
    /// Coelution: centered+normalized windows, row-major `[frag * win + j]`.
    cn: Vec<f32>,
    /// Coelution: indices of fragments with non-degenerate windows this cycle.
    active: Vec<usize>,
    /// Coelution: weighted sum of active unit-window vectors.
    wacc: Vec<f32>,
    /// Vote: per-cycle accumulated support across rows.
    support: Vec<f32>,
    /// Vote: baseline-subtracted, clipped row.
    r: Vec<f32>,
    /// Vote: matched-filtered row.
    m: Vec<f32>,
    /// Vote: selection buffer for median/MAD.
    tmp: Vec<f32>,
}

/// Active fragment rows: `(per-cycle XIC slice, expected-intensity weight)`.
/// Only fragments with positive expected intensity participate.
pub(crate) type Rows<'a> = Vec<(&'a [f32], f32)>;

/// Collect the active fragment rows from an extraction.
pub(crate) fn active_rows<'a, T: KeyLike>(ctx: &'a Extraction<T>) -> Rows<'a> {
    let mut rows: Rows<'a> = Vec::new();
    for ((key, _mz), chrom) in ctx.chromatograms.fragments.iter_mzs() {
        let exp = ctx.expected_intensities.get_fragment(key).unwrap_or(0.0);
        if exp > 0.0 {
            rows.push((chrom, exp));
        }
    }
    rows
}

/// Apply the enabled weight passes to `profile` in place, in order.
pub(crate) fn weight_profile<T: KeyLike>(
    profile: &mut [f32],
    ctx: &Extraction<T>,
    cfg: &ApexConfig,
    scratch: &mut ApexScratch,
) {
    if profile.is_empty() {
        return;
    }
    if cfg.coel_k <= 0.0 && cfg.vote_k <= 0.0 {
        return;
    }
    let rows = active_rows(ctx);
    coelution::weight(profile, &rows, cfg, scratch);
    vote::weight(profile, &rows, cfg, scratch);
}
