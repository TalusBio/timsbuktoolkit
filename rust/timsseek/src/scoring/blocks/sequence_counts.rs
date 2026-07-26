//! Sequence-derived features (features-only, conditionally present).
//!
//! Gated all-or-none on `peptide.aa_counts()`: either every candidate in a run
//! has a parsed sequence (Some) or none do. This makes the feature name-set
//! conditional — the gate adds `AA_COUNT_NAMES.len() + 2` trailing names
//! (`gate_delta_is_the_sequence_count_block` locks the delta). Emitted LAST so
//! those names stay at the tail.

use crate::models::AA_COUNT_NAMES;
use crate::models::sequence::Peptide;
use crate::scoring::blocks::{
    FrameSink,
    NameSink,
};

/// Push the 22 nonlinear-lane (tree-only) sequence feature *values*
/// (`peptide_length`, 20 `aa_count_*`, `peptide_n_mods`) iff the peptide has a
/// parsed sequence; otherwise nothing. The gate is speclib-wide, so within a
/// fit either every record emits these or none do. `counts` and
/// `AA_COUNT_NAMES` are both fixed `[_; 20]` arrays, so the zip is compile-time
/// guaranteed to cover all 20 counts.
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

/// Nonlinear-lane names for [`nonlinear_features`], same order.
pub fn nonlinear_feature_names(o: &mut NameSink) {
    o.push("peptide_length");
    for &n in AA_COUNT_NAMES.iter() {
        o.push(n);
    }
    o.push("peptide_n_mods");
}
