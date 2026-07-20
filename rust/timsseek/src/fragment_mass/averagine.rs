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

/// Averagine isotope envelope, renormalized to sum to 1.0.
///
/// `peptide_isotopes` itself normalizes to its tallest peak (relative
/// intensity, `M0 == 1.0`), matching its other call site
/// (`isotope_dist_from_seq`). The averagine estimate has no composition to
/// anchor a relative-intensity reading against, so it is rescaled here into a
/// genuine probability distribution over the reported isotopologues.
pub fn isotope_dist_from_mass(neutral_mass: f64) -> [f32; 3] {
    let (c, s) = averagine_cs_from_mass(neutral_mass);
    let raw = peptide_isotopes(c, s);
    let sum: f32 = raw.iter().sum();
    if sum > 0.0 {
        [raw[0] / sum, raw[1] / sum, raw[2] / sum]
    } else {
        raw
    }
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
    fn averagine_envelope_sums_to_one() {
        let env = isotope_dist_from_mass(1500.0);
        let sum: f32 = env.iter().sum();
        assert!((sum - 1.0).abs() < 1e-3, "env {env:?} sums to {sum}");
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
        assert!((env.iter().sum::<f32>() - 1.0).abs() < 1e-3);
    }
}
