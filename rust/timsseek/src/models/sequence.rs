//! Sequence-feature column names and the isotope helper's parsing fallback.

use timsquery::chemistry::ontologies;

/// Ontology-backed fallback for the existing isotope composition helper.
/// Search sequence features read stored analyte structure instead.
pub fn parse_proforma(
    sequence: &str,
) -> Result<mzcore::sequence::Peptidoform<mzcore::sequence::Linked>, String> {
    mzcore::sequence::Peptidoform::pro_forma(sequence, ontologies())
        .map(|(peptidoform, _warnings)| peptidoform)
        .map_err(|errors| {
            errors
                .iter()
                .map(|e| e.to_string())
                .collect::<Vec<_>>()
                .join("; ")
        })
}

pub const CANONICAL_AA_LETTERS: [u8; 20] = *b"ACDEFGHIKLMNPQRSTVWY";

/// Feature-vector names for the 20-dim AA-count block, derived from
/// [`CANONICAL_AA_LETTERS`] so the order can never drift out of sync.
/// `AA_COUNT_NAMES[i]` is `format!("aa_count_{}", CANONICAL_AA_LETTERS[i] as char)`
/// with a single one-time allocation leaked to `&'static str`.
pub static AA_COUNT_NAMES: std::sync::LazyLock<[&'static str; 20]> =
    std::sync::LazyLock::new(|| {
        let mut out: [&'static str; 20] = [""; 20];
        for (i, &c) in CANONICAL_AA_LETTERS.iter().enumerate() {
            let s = format!("aa_count_{}", c as char);
            out[i] = Box::leak(s.into_boxed_str());
        }
        out
    });

#[cfg(test)]
#[path = "sequence_legacy_tests.rs"]
mod legacy;
#[cfg(test)]
pub(crate) use legacy::parse_sequence;
#[cfg(test)]
pub(crate) use timsquery::chemistry::normalize_to_proforma;
