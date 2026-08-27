//! Sequence-string helpers shared by the library readers.

/// Strip modification annotations — anything inside `(...)`, `[...]` or
/// `{...}` — leaving the bare residue string.
///
/// One implementation for every reader: DIA-NN, Skyline and Spectronaut all
/// spell modifications differently but nest them the same way, and three
/// copies that differ only in which bracket pairs they recognise is how a
/// library ends up with a "stripped" sequence that still has a `{` in it.
/// Unbalanced closers are ignored rather than driving the depth negative, so a
/// malformed sequence degrades to a partial strip instead of dropping the tail.
pub fn strip_mods(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut depth: u32 = 0;
    for c in s.chars() {
        match c {
            '(' | '[' | '{' => depth += 1,
            ')' | ']' | '}' => depth = depth.saturating_sub(1),
            _ if depth == 0 => out.push(c),
            _ => {}
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strips_every_bracket_style() {
        assert_eq!(strip_mods("PEPTC[UNIMOD:4]IDEK"), "PEPTCIDEK");
        assert_eq!(strip_mods("PEPTC(UniMod:4)IDEK"), "PEPTCIDEK");
        assert_eq!(strip_mods("PEPTC{57.02}IDEK"), "PEPTCIDEK");
        assert_eq!(strip_mods("C[+57.02]AM"), "CAM");
        assert_eq!(strip_mods("PEPTIDEK"), "PEPTIDEK");
    }

    #[test]
    fn nesting_and_unbalanced_closers_do_not_lose_the_tail() {
        assert_eq!(strip_mods("PE[a[b]c]PTIDE"), "PEPTIDE");
        // A stray closer must not drive the depth negative and swallow the
        // rest of the sequence.
        assert_eq!(strip_mods("PEP]TIDEK"), "PEPTIDEK");
    }
}
