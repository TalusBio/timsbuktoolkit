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

pub use timsquery::chemistry::CANONICAL_AA_LETTERS;

#[cfg(test)]
mod tests {
    use super::*;

    /// Monoisotopic mass of a peptide with exactly one formula.
    fn neutral_mass(seq: &str) -> f64 {
        use mzcore::prelude::*;
        let pf = parse_proforma(seq).unwrap_or_else(|e| panic!("{seq:?}: {e}"));
        let linear = pf.as_linear().expect("linear").clone();
        let formulas = linear.formulas();
        assert_eq!(formulas.len(), 1, "{seq:?} is not a single formula");
        formulas[0].monoisotopic_mass().value
    }

    /// Compare against atomic masses, not mzcore output. This catches an
    /// accession mapped to the wrong chemistry.
    ///
    /// These four modifications cover the library inputs used here.
    #[test]
    fn unimod_deltas_match_their_compositions() {
        // (bare, modified, delta, composition)
        let cases = [
            ("PEPTCIDEK", "PEPTC[UNIMOD:4]IDEK", 57.021_46, "C2H3NO"),
            ("PEPTMIDEK", "PEPTM[UNIMOD:35]IDEK", 15.994_91, "O"),
            ("PEPTSIDEK", "PEPTS[UNIMOD:21]IDEK", 79.966_33, "HPO3"),
            ("PEPTNIDEK", "PEPTN[UNIMOD:7]IDEK", 0.984_02, "O minus NH"),
        ];
        for (bare, modified, delta, composition) in cases {
            let got = neutral_mass(modified) - neutral_mass(bare);
            assert!(
                (got - delta).abs() < 1e-4,
                "{modified} minus {bare} must be {delta} ({composition}), got {got}"
            );
        }
    }

    /// PEPTIDEK's residue sum is 909.44434 Da. Add H2O to get 927.45491 Da.
    /// The tolerance covers the five-decimal residue masses.
    #[test]
    fn unmodified_peptide_mass_is_the_hand_sum() {
        let got = neutral_mass("PEPTIDEK");
        assert!(
            (got - 927.454_91).abs() < 1e-4,
            "PEPTIDEK must be 927.45491, got {got}"
        );
    }
}
