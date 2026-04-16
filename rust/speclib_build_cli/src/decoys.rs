use std::collections::HashMap;
use std::sync::LazyLock;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DecoyMode {
    None,
    Reverse,
    EdgeMutate,
}

impl DecoyMode {
    pub fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "none" => Ok(Self::None),
            "reverse" => Ok(Self::Reverse),
            "edge_mutate" => Ok(Self::EdgeMutate),
            _ => Err(format!("Unknown decoy strategy: {s}")),
        }
    }
}

/// Substitution table ported from Python:
/// "GAVLIFMPWSCTYHKRQEND" → "LLLVVLLLLTSSSSLLNDQE"
pub static MUTATE_TABLE: LazyLock<HashMap<char, char>> = LazyLock::new(|| {
    let from = "GAVLIFMPWSCTYHKRQEND";
    let to = "LLLVVLLLLTSSSSLLNDQE";
    from.chars().zip(to.chars()).collect()
});

/// Keep first and last AA, reverse the middle.
/// "PEPTIDEK" → "PEDITPEK"
/// Sequences of <= 2 chars are returned as-is.
pub fn reverse_decoy(seq: &str) -> String {
    if seq.len() <= 2 {
        return seq.to_string();
    }
    let chars: Vec<char> = seq.chars().collect();
    let first = chars[0];
    let last = *chars.last().unwrap();
    let middle: String = chars[1..chars.len() - 1].iter().rev().collect();
    format!("{first}{middle}{last}")
}

/// Mutate positions 1 and -2 using MUTATE_TABLE.
/// Keep first (0) and last (-1) unchanged.
/// Sequences of <= 3 chars are returned as-is.
pub fn edge_mutate(seq: &str) -> String {
    if seq.len() <= 3 {
        return seq.to_string();
    }
    let mut chars: Vec<char> = seq.chars().collect();
    let n = chars.len();
    chars[1] = *MUTATE_TABLE.get(&chars[1]).unwrap_or(&chars[1]);
    chars[n - 2] = *MUTATE_TABLE.get(&chars[n - 2]).unwrap_or(&chars[n - 2]);
    chars.into_iter().collect()
}

/// Dispatch decoy generation. Panics if called with DecoyMode::None.
pub fn generate_decoy(seq: &str, mode: DecoyMode) -> String {
    match mode {
        DecoyMode::None => panic!("generate_decoy called with DecoyMode::None"),
        DecoyMode::Reverse => reverse_decoy(seq),
        DecoyMode::EdgeMutate => edge_mutate(seq),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_matches_python() {
        assert_eq!(reverse_decoy("PEPTIDEK"), "PEDITPEK");
        assert_eq!(reverse_decoy("PEPTIDEPINK"), "PNIPEDITPEK");
    }

    #[test]
    fn test_edge_mutate() {
        let result = edge_mutate("PEPTIDEK");
        assert_eq!(result.chars().nth(0).unwrap(), 'P'); // first preserved
        assert_ne!(result.chars().nth(1).unwrap(), 'E'); // mutated
        assert_eq!(result.chars().last().unwrap(), 'K'); // last preserved
        // E maps to D in MUTATE_TABLE (from "GAVLIFMPWSCTYHKRQEND" → "LLLVVLLLLTSSSSLLNDQE",
        // position 17: E→D), so position 1 should be 'D'
        assert_eq!(result.chars().nth(1).unwrap(), 'D');
    }

    #[test]
    fn test_mutate_table_completeness() {
        for aa in "GAVLIFMPWSCTYHKRQEND".chars() {
            assert!(MUTATE_TABLE.contains_key(&aa), "Missing mapping for {aa}");
        }
    }

    #[test]
    fn test_generate_decoy_reverse() {
        assert_eq!(generate_decoy("PEPTIDEK", DecoyMode::Reverse), "PEDITPEK");
    }

    #[test]
    fn test_generate_decoy_edge_mutate() {
        let result = generate_decoy("PEPTIDEK", DecoyMode::EdgeMutate);
        assert_ne!(result, "PEPTIDEK");
        assert_eq!(result.len(), "PEPTIDEK".len());
    }

    #[test]
    fn test_reverse_short_sequences() {
        assert_eq!(reverse_decoy(""), "");
        assert_eq!(reverse_decoy("A"), "A");
        assert_eq!(reverse_decoy("AK"), "AK");
        assert_eq!(reverse_decoy("AEK"), "AEK"); // 3 chars: first=A, middle=E reversed=E, last=K
    }

    #[test]
    fn test_edge_mutate_short_sequences() {
        assert_eq!(edge_mutate(""), "");
        assert_eq!(edge_mutate("A"), "A");
        assert_eq!(edge_mutate("AK"), "AK");
        assert_eq!(edge_mutate("AEK"), "AEK");
    }

    #[test]
    fn test_decoy_mode_from_str() {
        assert_eq!(DecoyMode::from_str("none").unwrap(), DecoyMode::None);
        assert_eq!(DecoyMode::from_str("reverse").unwrap(), DecoyMode::Reverse);
        assert_eq!(
            DecoyMode::from_str("edge_mutate").unwrap(),
            DecoyMode::EdgeMutate
        );
        assert!(DecoyMode::from_str("bogus").is_err());
    }

    #[test]
    #[should_panic(expected = "generate_decoy called with DecoyMode::None")]
    fn test_generate_decoy_none_panics() {
        generate_decoy("PEPTIDEK", DecoyMode::None);
    }

    #[test]
    fn test_reverse_three_chars() {
        // "AEK": first=A, middle=[E] reversed=[E], last=K → "AEK"
        assert_eq!(reverse_decoy("AEK"), "AEK");
        // "AEPK": first=A, middle=[E,P] reversed=[P,E], last=K → "APEK"
        assert_eq!(reverse_decoy("AEPK"), "APEK");
    }
}
