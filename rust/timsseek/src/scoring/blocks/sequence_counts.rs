//! Sequence-derived features (nonlinear lane, features-only -- no Parquet
//! column).
//!
//! UNCONDITIONAL: always [`LEN`] features wide. A peptide with no parsed
//! sequence emits `f64::NAN` for all of them rather than emitting nothing --
//! NaN is exactly what forust reads as "missing", and a fixed width is what
//! lets this block's contribution to the feature matrix be a compile-time
//! constant like every other block's. Emitted LAST, so these names stay at the
//! tail.

use crate::models::sequence::{
    AA_COUNT_NAMES,
    CANONICAL_AA_LETTERS,
    Peptide,
};
use crate::scoring::blocks::NameSink;

/// `peptide_length`, one count per canonical amino acid, `peptide_n_mods`.
/// Sized off [`CANONICAL_AA_LETTERS`] rather than `AA_COUNT_NAMES` because the
/// latter is a `LazyLock` (its `len()` is not const), and the two are the same
/// 20 residues by construction.
pub const LEN: usize = CANONICAL_AA_LETTERS.len() + 2;

/// The [`LEN`] nonlinear-lane (tree-only) sequence feature *values*, or all
/// `f64::NAN` when the peptide has no parsed sequence. `counts` and
/// `CANONICAL_AA_LETTERS` are both fixed 20-element arrays, so the middle
/// slice is compile-time guaranteed to cover all 20 counts.
pub fn nonlinear_feature_array(peptide: &Peptide) -> [f64; LEN] {
    let mut out = [f64::NAN; LEN];
    if let Some(parsed) = peptide.parse() {
        let counts = parsed.aa_counts();
        out[0] = parsed.residues.len() as f64;
        out[1..1 + counts.len()].copy_from_slice(&counts);
        out[LEN - 1] = parsed.mods.len() as f64;
    }
    out
}

/// Nonlinear-lane names for [`nonlinear_feature_array`], same order.
pub fn nonlinear_feature_names(o: &mut NameSink) {
    o.push("peptide_length");
    for &n in AA_COUNT_NAMES.iter() {
        o.push(n);
    }
    o.push("peptide_n_mods");
}
