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

    for (elem, count, _extras) in form.elements() {
        match (elem, count) {
            (&Element::C, Some(cnt)) => {
                ncarbon += cnt.get();
            }
            (&Element::S, Some(cnt)) => {
                nsulphur += cnt.get();
            }
            _ => {}
        }
    }

    (ncarbon, nsulphur)
}

pub fn count_carbon_sulphur_in_sequence(sequence: &str) -> Result<(u16, u16), String> {
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
