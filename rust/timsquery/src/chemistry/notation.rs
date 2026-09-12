/// Expand numeric `U:` accessions; preserve named modifications.
fn expand_unimod_accessions(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut rest = s;
    while let Some(at) = rest.to_ascii_lowercase().find("[u:") {
        let body = &rest[at + 3..];
        let numeric = body
            .split(']')
            .next()
            .is_some_and(|tag| !tag.is_empty() && tag.bytes().all(|b| b.is_ascii_digit()));
        out.push_str(&rest[..at]);
        out.push_str(if numeric { "[UNIMOD:" } else { "[U:" });
        rest = body;
    }
    out.push_str(rest);
    out
}

/// Replace every occurrence of `needle` in `haystack`, matching ignoring ASCII
/// case, with `replacement`. `needle` must be ASCII (UNIMOD tags are). ASCII-only
/// lowercasing preserves byte length, so match indices stay aligned with the
/// original (UTF-8-safe) string.
fn replace_ascii_ci(haystack: &str, needle: &str, replacement: &str) -> String {
    debug_assert!(needle.is_ascii());
    let hay_lower = haystack.to_ascii_lowercase();
    let needle_lower = needle.to_ascii_lowercase();
    let mut out = String::with_capacity(haystack.len());
    let mut i = 0;
    while i < haystack.len() {
        if hay_lower[i..].starts_with(&needle_lower) {
            out.push_str(replacement);
            i += needle.len();
        } else {
            let ch = haystack[i..].chars().next().unwrap();
            out.push(ch);
            i += ch.len_utf8();
        }
    }
    out
}

/// Convert DIA-NN parenthesised UNIMOD mods to ProForma brackets: each
/// case-insensitive `(unimod:` / `(u:` opener becomes `[UNIMOD:` and the single
/// `)` that closes it becomes `]`. Any other paren -- notably `)` inside a bracket
/// mod name like `[Carbamidomethyl (C)]` -- is left untouched.
fn convert_paren_unimod(s: &str) -> String {
    let lower = s.to_ascii_lowercase();
    let mut out = String::with_capacity(s.len());
    let mut i = 0;
    while i < s.len() {
        let opener_len = if lower[i..].starts_with("(unimod:") {
            Some("(unimod:".len())
        } else if lower[i..].starts_with("(u:") {
            Some("(u:".len())
        } else {
            None
        };
        match opener_len {
            Some(len) => {
                out.push_str("[UNIMOD:");
                i += len;
                let rest = &s[i..];
                match rest.find(')') {
                    // The id (ASCII digits) up to the matching close paren.
                    Some(rel) => {
                        out.push_str(&rest[..rel]);
                        out.push(']');
                        i += rel + 1; // consume ')'
                    }
                    // Malformed (no closer) -- copy the remainder verbatim.
                    None => {
                        out.push_str(rest);
                        i = s.len();
                    }
                }
            }
            None => {
                let ch = s[i..].chars().next().unwrap();
                out.push(ch);
                i += ch.len_utf8();
            }
        }
    }
    out
}

/// Normalize library sequence notation: strip `_` wrappers, convert DIA-NN
/// parentheses, normalize UNIMOD accessions, and add missing terminal dashes.
///
/// ProForma 2.1 requires full prefixes for accessions: `[U:21]` is invalid,
/// but `M[U:Oxidation]` is valid because `Oxidation` is a name (§6.2.1–2, §7.8).
/// This adapter repairs numeric shorthand; it is not a ProForma validator.
///
/// ```
/// use timsquery::chemistry::normalize_to_proforma;
///
/// assert_eq!(normalize_to_proforma("S[U:21]"), "S[UNIMOD:21]");
/// assert_eq!(normalize_to_proforma("M[U:Oxidation]"), "M[U:Oxidation]");
/// ```
///
/// Other supported library spellings:
/// ```
/// use timsquery::chemistry::normalize_to_proforma;
///
/// for (input, expected) in [
///     ("_PEPTIDEK_", "PEPTIDEK"),
///     ("_LSHPGC[UniMod:4]K_", "LSHPGC[UNIMOD:4]K"),
///     ("_C[Unimod:4]TVPGHK_", "C[UNIMOD:4]TVPGHK"),
///     ("PEPTC[U:4]IDEK", "PEPTC[UNIMOD:4]IDEK"),
///     ("AAC(UniMod:4)DEK", "AAC[UNIMOD:4]DEK"),
///     ("AAC(unimod:4)DEK", "AAC[UNIMOD:4]DEK"),
///     ("AAC(UNIMOD:4)DEK", "AAC[UNIMOD:4]DEK"),
///     ("AAC(UniMod:4)M(UniMod:35)K", "AAC[UNIMOD:4]M[UNIMOD:35]K"),
///     ("(UniMod:1)AACDEK", "[UNIMOD:1]-AACDEK"),
///     ("AAC(UniMod:4)M[Oxidation (M)]K", "AAC[UNIMOD:4]M[Oxidation (M)]K"),
///     ("_C[Carbamidomethyl (C)]PEPK_", "C[Carbamidomethyl (C)]PEPK"),
///     ("[+42]AACDEK", "[+42]-AACDEK"),
///     ("PEPTM[+15.995]IDEK", "PEPTM[+15.995]IDEK"),
///     ("PEPTIDEK", "PEPTIDEK"),
///     ("PEPTC[u:4]IDEK", "PEPTC[UNIMOD:4]IDEK"),
///     ("PEPTC[U:Carbamidomethyl]IDEK", "PEPTC[U:Carbamidomethyl]IDEK"),
/// ] {
///     assert_eq!(normalize_to_proforma(input), expected);
/// }
/// ```
pub fn normalize_to_proforma(raw: &str) -> String {
    let trimmed = raw.trim_matches('_');
    // Fast path: plain-AA sequences (no mod tags) skip the rewrite chain.
    if !trimmed.contains('[') && !trimmed.contains('(') {
        return trimmed.to_owned();
    }

    // Preserve parentheses inside bracketed modification names.
    let s = convert_paren_unimod(trimmed);

    // Normalize any pre-existing bracket casing likewise.
    let mut s = replace_ascii_ci(&s, "[unimod:", "[UNIMOD:");
    s = expand_unimod_accessions(&s);

    // A mod at the very start is N-terminal: ProForma wants `[UNIMOD:n]-SEQ`.
    if s.starts_with('[')
        && let Some(close) = s.find(']')
        && s[close + 1..].chars().next().is_some_and(|c| c != '-')
    {
        s.insert(close + 1, '-');
    }
    s
}
