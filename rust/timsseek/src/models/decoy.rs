use serde::Serialize;

/// The different labels that denote if a sequence is a decoy or not.
///
/// NOTE: The main difference between the decoy and reversed decoy is that the reversed decoy
/// has already been reversed, thus converting it to a string can be done as-is.
#[derive(Debug, Clone, Copy, Serialize, PartialEq, Eq, std::hash::Hash, PartialOrd, Ord)]
pub enum DecoyMarking {
    Target,
    NonReversedDecoy,
    ReversedDecoy,
}

impl DecoyMarking {
    pub fn as_str(&self) -> &'static str {
        match self {
            DecoyMarking::Target => "Target",
            DecoyMarking::NonReversedDecoy => "Decoy",
            DecoyMarking::ReversedDecoy => "Decoy",
        }
    }

    pub fn is_decoy(&self) -> bool {
        match self {
            DecoyMarking::Target => false,
            DecoyMarking::NonReversedDecoy => true,
            DecoyMarking::ReversedDecoy => true,
        }
    }

    pub fn is_target(&self) -> bool {
        !self.is_decoy()
    }
}

/// Helper function to convert a sequence into its decoy form
pub(crate) fn as_decoy_string(sequence: &str) -> String {
    let mut sequence = sequence.to_string();
    let inner_rev = sequence[1..(sequence.len() - 1)]
        .chars()
        .rev()
        .collect::<String>();
    sequence.replace_range(1..(sequence.len() - 1), &inner_rev);
    sequence
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decoy() {
        let sequence = "PEPTIDEPINK";
        let decoy = as_decoy_string(sequence);
        assert_eq!(sequence, "PEPTIDEPINK");
        assert_eq!(decoy, "PNIPEDITPEK");
    }
}
