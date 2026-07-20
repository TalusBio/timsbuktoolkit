use crate::isotopes::peptide_isotopes;
use rustyms::prelude::{
    Element,
    MolecularFormula,
    Peptidoform,
};

/// Super simple 1/k0 prediction.
///
/// Refit on `hela_iccoff_gt20peps` (34k IDs, holdout MAPE 1.36%) after the
/// mobility unit-bug fix. Use `scripts/refit_mobility.py` against a fresh
/// `results.parquet` to refit on a different dataset.
///
/// Example:
/// ```
/// use timsseek::fragment_mass::elution_group_converter::supersimpleprediction;
/// let mass = 1810.917339999999;
/// let charge = 2;
/// let out = supersimpleprediction(mass / charge as f64, charge);
/// assert!((out - 1.144405).abs() < 0.001);
/// ```
pub fn supersimpleprediction(mz: f64, charge: i32) -> f64 {
    let intercept_ = -1.319388e+00;
    let log1p_mz = (mz + 1.).ln();
    let sq_mz_over_charge = mz.powi(2) / charge as f64;
    let log1p_sq_mz_over_charge = (sq_mz_over_charge + 1.).ln();

    intercept_
        + (-2.954677e-01 * log1p_mz)
        + (-9.277763e-05 * mz)
        + (3.219103e-01 * log1p_sq_mz_over_charge)
        + (4.005229e-07 * sq_mz_over_charge)
        + (1.176651e-01 * charge as f64)
}

fn count_carbon_sulphur(form: &MolecularFormula) -> (u16, u16) {
    let mut ncarbon = 0;
    let mut nsulphur = 0;

    // `elements()` yields (element, isotope, count). The middle field is the
    // isotope (nucleon number, `None` for the natural/unspecified isotope), NOT
    // the atom count — the count is the third field. `count` is i32 and can be
    // negative for a loss; clamp before the u16 cast.
    for (elem, _isotope, count) in form.elements() {
        let n = (*count).max(0) as u16;
        match elem {
            Element::C => ncarbon += n,
            Element::S => nsulphur += n,
            _ => {}
        }
    }

    (ncarbon, nsulphur)
}

/// In-chain (C, S) atom counts per standard residue, indexed by `byte - b'A'`.
/// `None` = a non-standard code (B/J/O/U/X/Z) — defer to the rustyms path.
///
/// A residue contributes the same carbon/sulfur as its free amino acid: forming
/// a peptide bond removes one water per bond and the terminal water carries
/// neither C nor S, so a bare-sequence sum equals rustyms' formula exactly.
const RESIDUE_CS: [Option<(u16, u16)>; 26] = {
    // Alphabet offset of an uppercase residue byte (as a fn so `b'A'` maps to 0
    // without a literal `b'A' - b'A'`, which clippy's eq_op denies).
    const fn ri(c: u8) -> usize {
        (c - b'A') as usize
    }
    let mut t = [None; 26];
    t[ri(b'A')] = Some((3, 0));
    t[ri(b'C')] = Some((3, 1));
    t[ri(b'D')] = Some((4, 0));
    t[ri(b'E')] = Some((5, 0));
    t[ri(b'F')] = Some((9, 0));
    t[ri(b'G')] = Some((2, 0));
    t[ri(b'H')] = Some((6, 0));
    t[ri(b'I')] = Some((6, 0));
    t[ri(b'K')] = Some((6, 0));
    t[ri(b'L')] = Some((6, 0));
    t[ri(b'M')] = Some((5, 1));
    t[ri(b'N')] = Some((4, 0));
    t[ri(b'P')] = Some((5, 0));
    t[ri(b'Q')] = Some((5, 0));
    t[ri(b'R')] = Some((6, 0));
    t[ri(b'S')] = Some((3, 0));
    t[ri(b'T')] = Some((4, 0));
    t[ri(b'V')] = Some((5, 0));
    t[ri(b'W')] = Some((11, 0));
    t[ri(b'Y')] = Some((9, 0));
    t
};

/// Fast (C, S) tally over a bare amino-acid sequence via [`RESIDUE_CS`].
/// `None` on an empty string or any non-standard residue, forcing the rustyms
/// fallback so behavior (including the error path) is preserved.
fn count_cs_fast(sequence: &str) -> Option<(u16, u16)> {
    if sequence.is_empty() {
        return None;
    }
    let mut ncarbon = 0u16;
    let mut nsulphur = 0u16;
    for &b in sequence.as_bytes() {
        let idx = b.wrapping_sub(b'A') as usize;
        let (c, s) = RESIDUE_CS.get(idx).copied().flatten()?;
        ncarbon += c;
        nsulphur += s;
    }
    Some((ncarbon, nsulphur))
}

/// (C, S) counts for `sequence` (a bare, mod-stripped peptide on the hot path).
/// Tries the allocation-free table first; defers to the rustyms formula path for
/// empty / non-standard input, which stays the authority.
pub fn count_carbon_sulphur_in_sequence(sequence: &str) -> Result<(u16, u16), String> {
    if let Some(cs) = count_cs_fast(sequence) {
        return Ok(cs);
    }
    count_carbon_sulphur_in_sequence_rustyms(sequence)
}

fn count_carbon_sulphur_in_sequence_rustyms(sequence: &str) -> Result<(u16, u16), String> {
    let peptide = match Peptidoform::pro_forma(sequence, None) {
        Ok(pep) => pep,
        Err(e) => {
            return Err(format!(
                "Error parsing peptide sequence {}: {:?}",
                sequence, e
            ));
        }
    };
    let peptide = match peptide.as_linear() {
        Some(pep) => pep,
        None => return Err("Peptide is not linear.".to_string()),
    }
    .clone();

    let pep_formulas = peptide.formulas();
    if pep_formulas.len() > 1 {
        return Err("Peptide contains more than one formula.".to_string());
    }
    let form = pep_formulas[0].clone();
    Ok(count_carbon_sulphur(&form))
}

pub fn isotope_dist_from_seq(sequence: &str) -> Result<[f32; 3], String> {
    let (ncarbon, nsulphur) = count_carbon_sulphur_in_sequence(sequence)?;
    Ok(peptide_isotopes(ncarbon, nsulphur))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cs_table_matches_rustyms_per_residue() {
        // Every standard residue: the table must equal rustyms' formula count.
        for &aa in b"ACDEFGHIKLMNPQRSTVWY" {
            let seq = String::from_utf8(vec![aa, aa, aa]).unwrap(); // e.g. "AAA"
            let fast = count_cs_fast(&seq).expect("standard residue in table");
            let slow = count_carbon_sulphur_in_sequence_rustyms(&seq)
                .unwrap_or_else(|e| panic!("rustyms failed on {seq}: {e}"));
            assert_eq!(fast, slow, "C/S mismatch for {seq}");
        }
    }

    #[test]
    fn cs_table_matches_rustyms_on_peptides() {
        for seq in ["AAAGAAATHLEVAR", "LEGNSPQGSNQGVK", "MCMCMCK", "PEPTIDEK"] {
            let fast = count_cs_fast(seq).expect("standard peptide");
            let slow = count_carbon_sulphur_in_sequence_rustyms(seq).unwrap();
            assert_eq!(fast, slow, "C/S mismatch for {seq}");
        }
    }

    #[test]
    fn cs_fast_defers_on_nonstandard_and_empty() {
        assert!(count_cs_fast("").is_none(), "empty must defer");
        assert!(count_cs_fast("PEPXK").is_none(), "X must defer");
        assert!(count_cs_fast("pepk").is_none(), "lowercase must defer");
    }
}
