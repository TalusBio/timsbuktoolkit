use super::elution_group_converter::count_carbon_sulphur_in_sequence;
use crate::isotopes::peptide_isotopes;

/// Which model produced an isotope envelope, for load-time reporting.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IsotopeSource {
    Composition,
    Averagine,
}

// Senko averagine residue (avg amino acid): C4.9384 H7.7583 N1.3577 O1.4773 S0.0417,
// average residue mass ~111.1054 Da. Per-Dalton element counts:
const C_PER_DA: f64 = 4.9384 / 111.1054;
const S_PER_DA: f64 = 0.0417 / 111.1054;

pub fn averagine_cs_from_mass(neutral_mass: f64) -> (u16, u16) {
    let m = neutral_mass.max(0.0);
    let c = (m * C_PER_DA).round().clamp(0.0, u16::MAX as f64) as u16;
    let s = (m * S_PER_DA).round().clamp(0.0, u16::MAX as f64) as u16;
    (c, s)
}

/// Averagine isotope envelope: relative intensity, tallest peak == 1.0.
///
/// Matches `peptide_isotopes`'s own normalization (max peak, not sum), which
/// is what the `Composition` branch of `isotope_dist_or_averagine` returns
/// verbatim. Both isotope sources must share this scale so they're
/// interchangeable at scoring time.
pub fn isotope_dist_from_mass(neutral_mass: f64) -> [f32; 3] {
    let (c, s) = averagine_cs_from_mass(neutral_mass);
    peptide_isotopes(c, s)
}

/// Composition envelope when the sequence is countable, else averagine from mass.
pub fn isotope_dist_or_averagine(seq: &str, neutral_mass: f64) -> (IsotopeSource, [f32; 3]) {
    match count_carbon_sulphur_in_sequence(seq) {
        Ok((c, s)) => (IsotopeSource::Composition, peptide_isotopes(c, s)),
        Err(_) => (
            IsotopeSource::Averagine,
            isotope_dist_from_mass(neutral_mass),
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn averagine_cs_grows_with_mass() {
        let (c1, _s1) = averagine_cs_from_mass(1000.0);
        let (c2, _s2) = averagine_cs_from_mass(2000.0);
        assert!(c2 > c1, "carbon count must grow with mass");
        // ~0.0444 C/Da -> ~44 C at 1000 Da
        assert!((c1 as i32 - 44).abs() <= 2, "got {c1} C at 1000 Da");
    }

    #[test]
    fn averagine_envelope_is_max_normalized() {
        let env = isotope_dist_from_mass(1500.0);
        let max = env.iter().copied().fold(f32::MIN, f32::max);
        assert!((max - 1.0).abs() < 1e-4, "env {env:?} max is {max}");
        assert!(
            env.iter().all(|&v| (0.0..=1.0).contains(&v)),
            "env {env:?} has a value outside [0, 1]"
        );
    }

    #[test]
    fn or_averagine_uses_composition_for_standard_peptide() {
        let (src, _env) = isotope_dist_or_averagine("PEPTIDEK", 900.4);
        assert_eq!(src, IsotopeSource::Composition);
    }

    #[test]
    fn or_averagine_falls_back_on_nonstandard() {
        // `B` (Asx) is genuinely ambiguous between Asp/Asn in rustyms and
        // resolves to more than one formula, which is the real trigger for
        // the rustyms-backed count path to error today. (`X` was tried first
        // but rustyms resolves it to a defined zero-C/S formula rather than
        // erroring, so it does not exercise the fallback.)
        let (src, env) = isotope_dist_or_averagine("PEPBK", 600.0);
        assert_eq!(src, IsotopeSource::Averagine);
        let max = env.iter().copied().fold(f32::MIN, f32::max);
        assert!((max - 1.0).abs() < 1e-4, "env {env:?} max is {max}");
    }
}
