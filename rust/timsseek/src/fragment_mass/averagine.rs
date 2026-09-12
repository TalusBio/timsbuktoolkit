use crate::isotopes::peptide_isotopes;

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
/// Uses the same C/S calculator and normalization as composition-derived counts.
pub fn isotope_dist_from_mass(neutral_mass: f64) -> [f32; 3] {
    let (c, s) = averagine_cs_from_mass(neutral_mass);
    peptide_isotopes(c, s)
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
}
