//! Sequence-derived features (features-only, conditionally present).
//!
//! Gated all-or-none on `peptide.aa_counts()`: either every candidate in a run
//! has a parsed sequence (Some) or none do. This makes the feature name-set
//! conditional — the golden test asserts two sets: 86 (gate off) and 108
//! (gate on). Emitted LAST so the trailing 22 names are stable.

use crate::models::AA_COUNT_NAMES;
use crate::models::sequence::Peptide;
use crate::scoring::blocks::FeatSink;

/// Push the 22 sequence features (`peptide_length`, 20 `aa_count_*`,
/// `peptide_n_mods`) iff the peptide has a parsed sequence; otherwise nothing.
pub fn features(peptide: &Peptide, o: &mut FeatSink) {
    if let Some(counts) = peptide.aa_counts() {
        let length = peptide.length().unwrap() as f64;
        let n_mods = peptide.n_mods().unwrap() as f64;
        o.push("peptide_length", length);
        for (i, c) in counts.iter().enumerate() {
            o.push(AA_COUNT_NAMES[i], *c);
        }
        o.push("peptide_n_mods", n_mods);
    }
}
