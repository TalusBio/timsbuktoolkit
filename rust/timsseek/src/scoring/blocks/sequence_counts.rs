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
pub fn nonlinear_feature_array(
    peptide: Option<timsquery::chemistry::analyte::PeptideRef<'_>>,
) -> [f64; LEN] {
    let mut out = [f64::NAN; LEN];
    if let Some(peptide) = peptide
        && let Some(modifications) = peptide.modifications.known()
    {
        let mut counts = [0.0; 26];
        for residue in peptide.residues.bytes().filter(u8::is_ascii_uppercase) {
            counts[(residue - b'A') as usize] += 1.0;
        }
        out[0] = peptide.residues.len() as f64;
        for (i, &aa) in CANONICAL_AA_LETTERS.iter().enumerate() {
            out[i + 1] = counts[(aa - b'A') as usize];
        }
        out[LEN - 1] = modifications.len() as f64;
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

#[cfg(test)]
mod tests {
    use super::*;
    use timsquery::chemistry::analyte::Analyte;
    use timsquery::models::capabilities::DecoyPolicy;
    use timsquery::models::{
        Row,
        TargetCapabilities,
        TargetColumnsBuilder,
    };

    #[test]
    fn stored_counts_match_previous_parser_for_supported_sequences() {
        for sequence in [
            "PEPTIDEK",
            "_AC(UniMod:4)M(UniMod:35)K_",
            "[UNIMOD:1]-PEP[+15.99]TIDE-[UNIMOD:2]",
            "M[Oxidation]PEPTIDE",
            "PEPTIDE/2",
        ] {
            let normalized = crate::models::sequence::normalize_to_proforma(sequence);
            let previous = crate::models::sequence::parse_sequence(&normalized).unwrap();
            let analyte = Analyte::from_sequence(sequence);
            let mut builder = TargetColumnsBuilder::<timsquery::ion::IonAnnot>::with_capabilities(
                TargetCapabilities::default_diann(),
            );
            builder.push_row(Row {
                analyte: analyte.as_input(),
                ..Default::default()
            });
            let geom = builder.seal(DecoyPolicy::Never).unwrap();
            let values =
                nonlinear_feature_array(geom.analyte(geom.rows().next().unwrap()).peptide.known());
            assert_eq!(values[0], previous.residues.len() as f64, "{sequence}");
            assert_eq!(&values[1..21], &previous.aa_counts(), "{sequence}");
            assert_eq!(values[21], previous.mods.len() as f64, "{sequence}");
        }
    }
}
