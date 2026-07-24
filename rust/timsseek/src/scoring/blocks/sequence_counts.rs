//! Sequence-derived features (features-only, conditionally present).
//!
//! Gated all-or-none on `peptide.aa_counts()`: either every candidate in a run
//! has a parsed sequence (Some) or none do. This makes the feature name-set
//! conditional — the gate adds exactly 22 trailing names (`gate_delta_is_22_dims`
//! locks the delta). Emitted LAST so those names stay at the tail.

use crate::models::AA_COUNT_NAMES;
use crate::models::sequence::Peptide;
use crate::scoring::blocks::{
    FeatSink,
    FrameSink,
    NameSink,
};

/// Push the 22 sequence feature *values* (`peptide_length`, 20 `aa_count_*`,
/// `peptide_n_mods`) iff the peptide has a parsed sequence; otherwise nothing.
/// The gate is speclib-wide, so within a fit either every record emits these or
/// none do — see [`feature_names`] for the set-level names.
pub fn features(peptide: &Peptide, o: &mut FeatSink) {
    if let Some(counts) = peptide.aa_counts() {
        let length = peptide.length().unwrap() as f64;
        let n_mods = peptide.n_mods().unwrap() as f64;
        o.push(length);
        for c in counts.iter() {
            o.push(*c);
        }
        o.push(n_mods);
    }
}

/// The 22 sequence feature *names*, in [`features`] order. Emitted by the
/// set-level name builder only when the run's sequence gate is on; needs no
/// peptide because the name-set is fixed (drives directly off `AA_COUNT_NAMES`).
pub fn feature_names(o: &mut NameSink) {
    o.push("peptide_length");
    for &n in AA_COUNT_NAMES.iter() {
        o.push(n);
    }
    o.push("peptide_n_mods");
}

/// Nonlinear-lane (tree-only) variant of [`features`]: same 22 values, but
/// name-bound into a `FrameSink`. Conditional identically (emits iff parsed).
/// `counts` and `AA_COUNT_NAMES` are both fixed `[_; 20]` arrays, so the zip
/// is compile-time guaranteed to cover all 20 counts.
pub fn nonlinear_features(peptide: &Peptide, o: &mut FrameSink) {
    if let Some(counts) = peptide.aa_counts() {
        let length = peptide.length().unwrap() as f64;
        let n_mods = peptide.n_mods().unwrap() as f64;
        o.push("peptide_length", length);
        for (c, name) in counts.iter().zip(AA_COUNT_NAMES.iter()) {
            o.push(name, *c);
        }
        o.push("peptide_n_mods", n_mods);
    }
}

/// Nonlinear-lane names for [`nonlinear_features`] (same as [`feature_names`]).
pub fn nonlinear_feature_names(o: &mut NameSink) {
    o.push("peptide_length");
    for &n in AA_COUNT_NAMES.iter() {
        o.push(n);
    }
    o.push("peptide_n_mods");
}
