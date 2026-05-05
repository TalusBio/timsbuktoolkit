use serde::Serialize;

/// Target vs decoy marker on a peptide.
#[derive(Debug, Clone, Copy, Serialize, PartialEq, Eq, std::hash::Hash, PartialOrd, Ord)]
pub enum DecoyMarking {
    Target,
    ReversedDecoy,
}

impl DecoyMarking {
    pub fn as_str(&self) -> &'static str {
        match self {
            DecoyMarking::Target => "Target",
            DecoyMarking::ReversedDecoy => "Decoy",
        }
    }

    pub fn is_decoy(&self) -> bool {
        matches!(self, DecoyMarking::ReversedDecoy)
    }

    pub fn is_target(&self) -> bool {
        matches!(self, DecoyMarking::Target)
    }
}
