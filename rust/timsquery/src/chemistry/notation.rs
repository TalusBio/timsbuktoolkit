/// Rewrite `[U:<digits>]` to `[UNIMOD:<digits>]`, leaving `[U:<name>]` alone.
///
/// ProForma 2.1 pairs a name with the one-letter prefix (§6.2.1) and an
/// accession with the long one (§6.2.2), and calls the short-prefix accession
/// `[U:35]` incorrect outright. So this canonicalizes a spelling the spec
/// rejects into the one it wants, which is also the one `classify_mod` reads --
/// and leaves a name alone, since moving it to the long prefix would produce a
/// form the spec does not define.
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

/// Coerce DIA-NN / short-form modified-sequence strings into mzcore-parseable
/// ProForma. Strips `_..._` wrapping used by DIA-NN, converts DIA-NN's
/// parenthesised mods (`C(UniMod:4)`) to ProForma brackets (`C[UNIMOD:4]`), and
/// normalizes UNIMOD tag casing (`[UniMod:`, `[Unimod:`, `[U:` → `[UNIMOD:`).
/// A leading (N-terminal) mod is rewritten to `[UNIMOD:n]-SEQ` as ProForma
/// requires. Pass-through for plain sequences.
///
/// Off the hot path -- allocates on every replacement, which is fine at load.
pub fn normalize_to_proforma(raw: &str) -> String {
    let trimmed = raw.trim_matches('_');
    // Fast path: plain-AA sequences (no mod tags) skip the rewrite chain.
    if !trimmed.contains('[') && !trimmed.contains('(') {
        return trimmed.to_owned();
    }

    // DIA-NN writes mods in parentheses, e.g. `C(UniMod:4)`; ProForma needs
    // brackets, e.g. `C[UNIMOD:4]`. Convert each opener and ONLY its matching
    // `)` (see `convert_paren_unimod`) -- never a blanket `)` replace, which would
    // corrupt parens inside a bracket mod name, e.g. `C[Carbamidomethyl (C)]`.
    let s = convert_paren_unimod(trimmed);

    // Normalize any pre-existing bracket casing likewise.
    let mut s = replace_ascii_ci(&s, "[unimod:", "[UNIMOD:");
    // Expanded only for an accession, which is the form `classify_mod`'s fast
    // path reads. A NAME is left alone, because ProForma 2.1 keeps the two in
    // separate namespaces: §6.2.1 gives names a ONE-LETTER prefix, and §6.2.2
    // requires accessions to use the long one ("full accession numbers MUST be
    // used in all cases"). So `[UNIMOD:Carbamidomethyl]` is in neither, and
    // expanding a name produced a string no parser owes us an answer for.
    //
    // That is what made it a bug rather than a nicety: mzSpecLib names its
    // modifications, the gate is library-wide, and one unparsable row turned
    // sequence features off for a whole library.
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
