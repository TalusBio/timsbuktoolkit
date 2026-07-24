/// Mass difference between the 13C and 12C carbon isotopes (Da), i.e. the
/// spacing between adjacent peaks of a peptide isotope envelope. This is NOT
/// the free-neutron mass (1.00866 Da) — using that overstates the spacing by
/// ~0.0053 Da/step, biasing measured M+1/M+2 m/z errors several ppm negative.
/// Value: 13C = 13.0033548378 u (AME2020 / NIST), 12C ≡ 12; difference below.
pub const C13_C12_MASS_DIFF: f64 = 1.0033548378;
pub const PROTON_MASS: f64 = 1.007276466;
