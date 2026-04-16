/// Proforma-like modification parsing and application for speclib_build.
///
/// Supports notations like "C[U:4]", "M[U:35]", "S[U:21]", "M[+15.995]".
/// Fixed mods are inserted after every matching residue in the sequence.
/// Variable mods generate all combinations up to `max_mods` sites.

// ── Types ────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Modification {
    /// Uppercase single-letter amino acid code.
    pub target_residue: char,
    /// Bracket notation including brackets, e.g. "[U:4]" or "[+15.995]".
    pub notation: String,
}

impl Modification {
    /// Parse a modification string such as `"C[U:4]"` or `"M[+15.995]"`.
    ///
    /// The first character must be an uppercase ASCII letter (amino acid code).
    /// The remainder must start with `[` and end with `]`.
    pub fn parse(s: &str) -> Result<Self, String> {
        let mut chars = s.chars();
        let first = chars
            .next()
            .ok_or_else(|| "modification string is empty".to_string())?;

        if !first.is_ascii_uppercase() {
            return Err(format!(
                "expected uppercase amino acid letter, got '{first}'"
            ));
        }

        let rest: String = chars.collect();
        if !rest.starts_with('[') || !rest.ends_with(']') {
            return Err(format!(
                "bracket notation must start with '[' and end with ']', got '{rest}'"
            ));
        }
        if rest.len() < 3 {
            // minimum: "[x]"
            return Err(format!("bracket notation too short: '{rest}'"));
        }

        Ok(Self {
            target_residue: first,
            notation: rest,
        })
    }
}

// ── Fixed modifications ───────────────────────────────────────────────────────

/// Insert each modification's notation after every matching unmodified residue.
///
/// "Unmodified" means the residue is not already followed by a `[` bracket.
/// Fixed mods are applied left-to-right; each mod is applied in sequence
/// over the already-modified string (so you can stack multiple fixed mods).
pub fn apply_fixed_mods(sequence: &str, mods: &[Modification]) -> String {
    let mut result = sequence.to_string();
    for m in mods {
        result = apply_one_fixed_mod(&result, m);
    }
    result
}

fn apply_one_fixed_mod(sequence: &str, m: &Modification) -> String {
    let bytes = sequence.as_bytes();
    let target = m.target_residue as u8;
    let mut out = String::with_capacity(sequence.len() + sequence.len() / 4);

    let mut i = 0;
    while i < bytes.len() {
        let b = bytes[i];
        out.push(b as char);
        // Append notation only when: this byte matches the target AND the
        // immediately following character is NOT '[' (already modified).
        if b == target {
            let already_modified = bytes.get(i + 1).copied() == Some(b'[');
            if !already_modified {
                out.push_str(&m.notation);
            }
        }
        i += 1;
    }
    out
}

// ── Variable modifications ────────────────────────────────────────────────────

/// Generate all forms of `sequence` with 0..=`max_mods` variable modifications
/// applied.
///
/// Each element of `mods` specifies a residue and notation. The positions
/// considered are only unmodified sites (not already followed by `[`).
/// Returns at least the unmodified sequence (always the first element).
pub fn expand_variable_mods(
    sequence: &str,
    mods: &[Modification],
    max_mods: usize,
) -> Vec<String> {
    // Collect all modifiable (residue_index_in_sequence, &Modification) pairs.
    // We scan character-by-character, skipping bracket contents.
    let positions = collect_modifiable_positions(sequence, mods);

    if positions.is_empty() || max_mods == 0 {
        return vec![sequence.to_string()];
    }

    let cap = max_mods.min(positions.len());
    let mut results: Vec<String> = Vec::new();

    // k = 0 → unmodified form
    results.push(sequence.to_string());

    // k = 1 ..= cap
    for k in 1..=cap {
        let combos = combinations(&positions, k);
        for combo in combos {
            results.push(build_modified_sequence(sequence, &combo));
        }
    }

    results
}

/// A modifiable site: byte offset in `sequence` of the target residue and the
/// notation to insert after it.
#[derive(Clone)]
struct ModSite<'a> {
    /// Byte offset of the residue character in the original sequence string.
    byte_offset: usize,
    notation: &'a str,
}

/// Walk the sequence, honouring existing bracket groups, and collect all
/// positions where each mod could be applied.
fn collect_modifiable_positions<'a>(
    sequence: &str,
    mods: &'a [Modification],
) -> Vec<ModSite<'a>> {
    let mut sites: Vec<ModSite<'a>> = Vec::new();
    let bytes = sequence.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'[' {
            // Skip over bracket group.
            while i < bytes.len() && bytes[i] != b']' {
                i += 1;
            }
            // skip ']'
            i += 1;
            continue;
        }
        let ch = bytes[i] as char;
        // Check against each mod.
        for m in mods {
            if ch == m.target_residue {
                // Make sure the next char is not '[' (already modified).
                let already = bytes.get(i + 1).copied() == Some(b'[');
                if !already {
                    sites.push(ModSite {
                        byte_offset: i,
                        notation: &m.notation,
                    });
                }
            }
        }
        i += 1;
    }
    sites
}

/// Reconstruct the sequence string with modifications inserted at the chosen
/// sites. Sites must be in ascending `byte_offset` order.
fn build_modified_sequence(sequence: &str, sites: &[&ModSite<'_>]) -> String {
    let bytes = sequence.as_bytes();
    let mut out = String::with_capacity(sequence.len() + sites.len() * 8);
    let mut prev = 0usize;
    for site in sites {
        // Copy up to and including the residue.
        let end = site.byte_offset + 1;
        out.push_str(&sequence[prev..end]);
        out.push_str(site.notation);
        prev = end;
    }
    // Remainder of sequence.
    if prev < bytes.len() {
        out.push_str(&sequence[prev..]);
    }
    out
}

/// Generate all k-element combinations from `items`, preserving order.
fn combinations<'a>(items: &'a [ModSite<'a>], k: usize) -> Vec<Vec<&'a ModSite<'a>>> {
    let mut result = Vec::new();
    let mut combo = Vec::with_capacity(k);
    combinations_inner(items, k, 0, &mut combo, &mut result);
    result
}

fn combinations_inner<'a>(
    items: &'a [ModSite<'a>],
    k: usize,
    start: usize,
    combo: &mut Vec<&'a ModSite<'a>>,
    result: &mut Vec<Vec<&'a ModSite<'a>>>,
) {
    if combo.len() == k {
        result.push(combo.clone());
        return;
    }
    let remaining = k - combo.len();
    for i in start..=(items.len().saturating_sub(remaining)) {
        combo.push(&items[i]);
        combinations_inner(items, k, i + 1, combo, result);
        combo.pop();
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Parsing ──────────────────────────────────────────────────────────────

    #[test]
    fn test_parse_unimod_mod() {
        let m = Modification::parse("C[U:4]").unwrap();
        assert_eq!(m.target_residue, 'C');
        assert_eq!(m.notation, "[U:4]");
    }

    #[test]
    fn test_parse_mass_shift_mod() {
        let m = Modification::parse("M[+15.995]").unwrap();
        assert_eq!(m.target_residue, 'M');
        assert_eq!(m.notation, "[+15.995]");
    }

    #[test]
    fn test_parse_bad_first_char() {
        assert!(Modification::parse("1[U:4]").is_err());
        assert!(Modification::parse("[U:4]").is_err());
    }

    #[test]
    fn test_parse_bad_bracket() {
        assert!(Modification::parse("CU:4]").is_err());
        assert!(Modification::parse("C[U:4").is_err());
    }

    // ── Fixed mods ───────────────────────────────────────────────────────────

    #[test]
    fn test_apply_fixed_mod() {
        let m = Modification::parse("C[U:4]").unwrap();
        let result = apply_fixed_mods("PEPTCIDECK", &[m]);
        assert_eq!(result, "PEPTC[U:4]IDEC[U:4]K");
    }

    #[test]
    fn test_apply_fixed_mod_no_match() {
        let m = Modification::parse("C[U:4]").unwrap();
        let result = apply_fixed_mods("PEPTIDEK", &[m]);
        assert_eq!(result, "PEPTIDEK");
    }

    #[test]
    fn test_apply_fixed_mod_does_not_double_modify() {
        // If a residue is already bracketed it should not be modified again.
        let m = Modification::parse("C[U:4]").unwrap();
        let already = "PEPTC[U:4]IDECK";
        let result = apply_fixed_mods(already, &[m]);
        assert_eq!(result, "PEPTC[U:4]IDEC[U:4]K");
    }

    // ── Variable mods ────────────────────────────────────────────────────────

    #[test]
    fn test_expand_variable_mods_single() {
        let m = Modification::parse("M[U:35]").unwrap();
        let results = expand_variable_mods("PEPTMIDMEK", &[m], 1);
        assert!(results.contains(&"PEPTMIDMEK".to_string()));
        assert!(results.contains(&"PEPTM[U:35]IDMEK".to_string()));
        assert!(results.contains(&"PEPTMIDM[U:35]EK".to_string()));
        // max_mods=1 → doubly-modified form should NOT be present
        assert!(!results.contains(&"PEPTM[U:35]IDM[U:35]EK".to_string()));
        assert_eq!(results.len(), 3);
    }

    #[test]
    fn test_expand_variable_mods_max_two() {
        let m = Modification::parse("M[U:35]").unwrap();
        let results = expand_variable_mods("PEPTMIDMEK", &[m], 2);
        assert!(results.contains(&"PEPTMIDMEK".to_string()));
        assert!(results.contains(&"PEPTM[U:35]IDMEK".to_string()));
        assert!(results.contains(&"PEPTMIDM[U:35]EK".to_string()));
        assert!(results.contains(&"PEPTM[U:35]IDM[U:35]EK".to_string()));
        assert_eq!(results.len(), 4);
    }

    #[test]
    fn test_expand_variable_mods_phospho_sty() {
        // Multiple mod types: phospho on S and T.
        let phospho_s = Modification::parse("S[U:21]").unwrap();
        let phospho_t = Modification::parse("T[U:21]").unwrap();
        // Sequence has one S and one T → positions: S@4, T@6
        let results = expand_variable_mods("PEPSTIEK", &[phospho_s, phospho_t], 1);
        // Unmodified
        assert!(results.contains(&"PEPSTIEK".to_string()));
        // Only S modified
        assert!(results.contains(&"PEPS[U:21]TIEK".to_string()));
        // Only T modified
        assert!(results.contains(&"PEPST[U:21]IEK".to_string()));
        // max_mods=1 → doubly-modified absent
        assert!(!results.contains(&"PEPS[U:21]T[U:21]IEK".to_string()));
        assert_eq!(results.len(), 3);
    }

    #[test]
    fn test_expand_variable_mods_respects_fixed_mods() {
        // Sequence already has a fixed C[U:4] — variable M should still work.
        let var_m = Modification::parse("M[U:35]").unwrap();
        let sequence = "PEPTMC[U:4]IDMEK";
        let results = expand_variable_mods(sequence, &[var_m], 1);
        assert!(results.contains(&"PEPTMC[U:4]IDMEK".to_string()));
        assert!(results.contains(&"PEPTM[U:35]C[U:4]IDMEK".to_string()));
        assert!(results.contains(&"PEPTMC[U:4]IDM[U:35]EK".to_string()));
        assert_eq!(results.len(), 3);
    }

    #[test]
    fn test_expand_variable_mods_no_positions() {
        let m = Modification::parse("M[U:35]").unwrap();
        let results = expand_variable_mods("PEPTIDEK", &[m], 2);
        assert_eq!(results, vec!["PEPTIDEK".to_string()]);
    }

    #[test]
    fn test_expand_variable_mods_max_zero() {
        let m = Modification::parse("M[U:35]").unwrap();
        let results = expand_variable_mods("PEPTMIDMEK", &[m], 0);
        assert_eq!(results, vec!["PEPTMIDMEK".to_string()]);
    }
}
