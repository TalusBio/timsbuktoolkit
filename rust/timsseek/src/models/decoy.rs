use serde::Serialize;

/// Target vs decoy marker on a peptide.
#[derive(Debug, Clone, Copy, Serialize, PartialEq, Eq, std::hash::Hash, PartialOrd, Ord)]
pub enum DecoyMarking {
    Target,
    ReversedDecoy,    // sequence reversal
    MassShiftedDecoy, // same AA sequence, ±mass shift
}

impl DecoyMarking {
    pub fn as_str(&self) -> &'static str {
        match self {
            DecoyMarking::Target => "Target",
            DecoyMarking::ReversedDecoy | DecoyMarking::MassShiftedDecoy => "Decoy",
        }
    }

    pub fn is_decoy(&self) -> bool {
        !matches!(self, DecoyMarking::Target)
    }

    pub fn is_target(&self) -> bool {
        matches!(self, DecoyMarking::Target)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mass_shifted_is_a_decoy_but_not_reversed() {
        assert!(DecoyMarking::MassShiftedDecoy.is_decoy());
        assert!(!DecoyMarking::MassShiftedDecoy.is_target());
        assert!(DecoyMarking::ReversedDecoy.is_decoy());
        assert!(!DecoyMarking::Target.is_decoy());
    }

    #[test]
    fn both_decoy_kinds_print_decoy_for_output_compat() {
        assert_eq!(DecoyMarking::Target.as_str(), "Target");
        assert_eq!(DecoyMarking::ReversedDecoy.as_str(), "Decoy");
        assert_eq!(DecoyMarking::MassShiftedDecoy.as_str(), "Decoy");
    }
}
