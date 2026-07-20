//! Peptide and mod representation for per-entry speclib rows.
//!
//! `Peptide` is used on scoring paths. `ProteinSlice` stays on FASTA digestion.

use crate::models::decoy::DecoyMarking;
use serde::Serialize;
use smallvec::SmallVec;
use std::sync::Arc;

/// Amino acid stored as alphabet offset `c - b'A'` (0..=25). `u8::MAX`
/// means "unrecognized / non-alpha". Unreachable slots in count buffers
/// (B=1, J=9, O=14, U=20, X=23, Z=25) stay zero.
///
/// Chosen over canonical-20 index storage so parse cost is one subtraction
/// per residue (no 20-scan `.find`) and `aa_counts` inner loop is branchless.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AminoAcid(pub u8);

pub const UNKNOWN_AA: AminoAcid = AminoAcid(u8::MAX);

/// Canonical output order for the 20-dim count vector:
/// A C D E F G H I K L M N P Q R S T V W Y.
/// Do not reorder — `FEATURE_NAMES` follows this.
pub const CANONICAL_AA_INDICES: [usize; 20] = [
    0, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 15, 16, 17, 18, 19, 21, 22, 24,
];

pub const CANONICAL_AA_LETTERS: [u8; 20] = *b"ACDEFGHIKLMNPQRSTVWY";

/// Feature-vector names for the 20-dim AA-count block, derived from
/// [`CANONICAL_AA_LETTERS`] so the order can never drift out of sync.
/// `AA_COUNT_NAMES[i]` is `format!("aa_count_{}", CANONICAL_AA_LETTERS[i] as char)`
/// with a single one-time allocation leaked to `&'static str`.
pub static AA_COUNT_NAMES: std::sync::LazyLock<[&'static str; 20]> =
    std::sync::LazyLock::new(|| {
        let mut out: [&'static str; 20] = [""; 20];
        for (i, &c) in CANONICAL_AA_LETTERS.iter().enumerate() {
            let s = format!("aa_count_{}", c as char);
            out[i] = Box::leak(s.into_boxed_str());
        }
        out
    });

impl AminoAcid {
    pub fn from_ascii(c: u8) -> AminoAcid {
        if c.is_ascii_uppercase() {
            AminoAcid(c - b'A')
        } else {
            UNKNOWN_AA
        }
    }

    /// Returns `true` for any ASCII uppercase letter (`A`..=`Z`), including
    /// non-standard codes `B/J/O/U/X/Z`. These non-canonical codes occupy
    /// slots in the 26-wide count buffer but are never projected into the
    /// 20-dim `aa_counts` output. For "is this one of the canonical 20?",
    /// check `CANONICAL_AA_INDICES.contains(&(self.0 as usize))`.
    pub fn is_known(self) -> bool {
        self.0 < 26
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Mod {
    Unimod(u16),
    Mass(f32),
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModEntry {
    /// Residue position. Sentinel values: `POS_N_TERM` (0xFF) / `POS_C_TERM` (0xFE).
    /// Otherwise a 0-based residue index, valid range `0..=253`. Sequences
    /// longer than 254 residues cannot be represented — the parser rejects them
    /// in `parse_sequence` rather than truncating here.
    pub pos: u8,
    pub kind: Mod,
}

pub const POS_N_TERM: u8 = 0xFF;
pub const POS_C_TERM: u8 = 0xFE;

#[derive(Debug, Clone, PartialEq)]
pub struct ParsedSequence {
    pub residues: SmallVec<[AminoAcid; 32]>,
    pub mods: SmallVec<[ModEntry; 2]>,
}

impl ParsedSequence {
    pub fn aa_counts(&self) -> [f64; 20] {
        // Count into 26-wide alphabet buffer (branchless hot loop).
        let mut tmp = [0.0_f64; 26];
        for r in &self.residues {
            if r.is_known() {
                tmp[r.0 as usize] += 1.0;
            }
        }
        // Project to canonical 20-dim output.
        let mut out = [0.0_f64; 20];
        for (i, &idx) in CANONICAL_AA_INDICES.iter().enumerate() {
            out[i] = tmp[idx];
        }
        out
    }
}

#[derive(Debug, Clone)]
pub struct Peptide {
    pub raw: Arc<str>,
    pub parsed: Option<ParsedSequence>,
    pub decoy: DecoyMarking,
    pub decoy_group: u32,
}

impl Peptide {
    pub fn as_str(&self) -> &str {
        &self.raw
    }

    /// Byte length of the raw sequence string (not residue count).
    pub fn len(&self) -> usize {
        self.raw.len()
    }

    pub fn is_empty(&self) -> bool {
        self.raw.is_empty()
    }

    /// `true` when this entry is a decoy (reversed or mass-shifted).
    pub fn is_decoy(&self) -> bool {
        self.decoy.is_decoy()
    }

    pub fn length(&self) -> Option<u8> {
        self.parsed.as_ref().map(|p| p.residues.len() as u8)
    }

    pub fn aa_counts(&self) -> Option<[f64; 20]> {
        self.parsed.as_ref().map(|p| p.aa_counts())
    }

    pub fn n_mods(&self) -> Option<u8> {
        self.parsed.as_ref().map(|p| p.mods.len() as u8)
    }
}

impl Serialize for Peptide {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&self.raw)
    }
}

/// Speclib-level metadata. Lives on `Speclib` struct.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SeqFormat {
    /// Bare AA only, no mods (e.g. `PEPTIDEK`).
    Plain,
    /// Modified ProForma-able form (`[UNIMOD:4]`, `[+15.995]`, `_..._` OK).
    Modified,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpeclibMeta {
    pub parsable_sequences: bool,
    pub sequence_format: SeqFormat,
}

impl Default for SpeclibMeta {
    fn default() -> Self {
        Self {
            parsable_sequences: false,
            sequence_format: SeqFormat::Plain,
        }
    }
}

/// Parse a ProForma-normalized peptide string into our thin representation.
/// Returns `None` on parse error, non-linear peptides, or mods outside the
/// `Unimod`/`Mass` subset.
///
/// A hand-rolled byte walk handles the grammar `normalize_to_proforma` actually
/// emits (bare residues, `[UNIMOD:n]`, `[+/-mass]`, N-/C-terminal forms). It
/// returns `Some` ONLY for inputs it fully recognizes; anything else — named
/// mods, cross-links, unexpected bytes — yields `None` and defers to the rustyms
/// parser, which stays the authority for what is valid. So the fast path can
/// never accept something rustyms would reject, with one deliberate exception:
/// it does not check that a `UNIMOD:n` id exists in the ontology (a
/// syntactically valid id is accepted). Real DIA-NN output only carries real
/// ids, so this never triggers in practice.
pub fn parse_sequence(normalized: &str) -> Option<ParsedSequence> {
    if let Some(parsed) = parse_sequence_fast(normalized) {
        return Some(parsed);
    }
    parse_sequence_rustyms(normalized)
}

/// Classify one bracket body (`UNIMOD:n` or a signed mass like `+15.995`) into a
/// [`Mod`]. `None` for anything else — a named mod, an unsigned number, empty —
/// which forces the rustyms fallback in [`parse_sequence`].
fn classify_mod(body: &str) -> Option<Mod> {
    let body = body.trim();
    if body.len() >= 7 && body[..7].eq_ignore_ascii_case("UNIMOD:") {
        return body[7..].trim().parse::<u16>().ok().map(Mod::Unimod);
    }
    match body.as_bytes().first() {
        Some(b'+') | Some(b'-') => body.parse::<f32>().ok().map(Mod::Mass),
        _ => None,
    }
}

/// Byte-walk parser for the `normalize_to_proforma` output grammar. `None` means
/// "not recognized — defer to rustyms", never "definitively invalid" (that
/// verdict is the fallback's). See [`parse_sequence`] for the contract.
fn parse_sequence_fast(s: &str) -> Option<ParsedSequence> {
    let b = s.as_bytes();
    let mut residues: SmallVec<[AminoAcid; 32]> = SmallVec::new();
    let mut mods: SmallVec<[ModEntry; 2]> = SmallVec::new();
    let mut i = 0usize;

    // Optional leading N-terminal mod: `[..]-SEQ`.
    if b.first() == Some(&b'[') {
        let close = i + 1 + b[i + 1..].iter().position(|&c| c == b']')?;
        // A leading bracket not of the `[..]-` shape is something we do not
        // model; let rustyms decide.
        if b.get(close + 1) != Some(&b'-') {
            return None;
        }
        mods.push(ModEntry {
            pos: POS_N_TERM,
            kind: classify_mod(&s[i + 1..close])?,
        });
        i = close + 2;
    }

    while i < b.len() {
        match b[i] {
            c if c.is_ascii_uppercase() => {
                // A residue index must fit u8 for `ModEntry::pos` (0..=253).
                if residues.len() > 253 {
                    return None;
                }
                residues.push(AminoAcid::from_ascii(c));
                i += 1;
            }
            b'[' => {
                // Residue mod attaches to the residue just pushed.
                if residues.is_empty() {
                    return None;
                }
                let close = i + 1 + b[i + 1..].iter().position(|&c| c == b']')?;
                mods.push(ModEntry {
                    pos: (residues.len() - 1) as u8,
                    kind: classify_mod(&s[i + 1..close])?,
                });
                i = close + 1;
            }
            b'-' if b.get(i + 1) == Some(&b'[') => {
                // Trailing C-terminal mod: `SEQ-[..]`.
                let close = i + 2 + b[i + 2..].iter().position(|&c| c == b']')?;
                mods.push(ModEntry {
                    pos: POS_C_TERM,
                    kind: classify_mod(&s[i + 2..close])?,
                });
                i = close + 1;
            }
            _ => return None, // anything unexpected -> rustyms fallback
        }
    }

    if residues.is_empty() {
        return None;
    }
    Some(ParsedSequence { residues, mods })
}

/// The rustyms-backed parser. Authoritative fallback for [`parse_sequence`]:
/// validates against the ontology, handles named mods, and rejects non-linear
/// peptides. Off the hot path once the fast path covers the common grammar.
fn parse_sequence_rustyms(normalized: &str) -> Option<ParsedSequence> {
    use rustyms::prelude::IsAminoAcid;
    use rustyms::sequence::Peptidoform;

    let pf = Peptidoform::pro_forma(normalized, None).ok()?;
    let linear = pf.into_linear()?;

    let mut residues: SmallVec<[AminoAcid; 32]> = SmallVec::new();
    let mut mods: SmallVec<[ModEntry; 2]> = SmallVec::new();

    for (i, el) in linear.sequence().iter().enumerate() {
        if i > 253 {
            return None;
        }
        let c = el.aminoacid.pro_forma_definition().chars().next()?;
        residues.push(AminoAcid::from_ascii(c as u8));

        for m in &el.modifications {
            let kind = modification_to_mod(m)?;
            mods.push(ModEntry { pos: i as u8, kind });
        }
    }

    for m in linear.get_n_term() {
        let kind = modification_to_mod(m)?;
        mods.push(ModEntry {
            pos: POS_N_TERM,
            kind,
        });
    }
    for m in linear.get_c_term() {
        let kind = modification_to_mod(m)?;
        mods.push(ModEntry {
            pos: POS_C_TERM,
            kind,
        });
    }

    Some(ParsedSequence { residues, mods })
}

fn modification_to_mod(m: &rustyms::sequence::Modification) -> Option<Mod> {
    use rustyms::ontology::Ontology;
    use rustyms::sequence::{
        Modification,
        SimpleModificationInner,
    };
    let simple = match m {
        Modification::Simple(s) => s,
        _ => return None, // Cross-link / ambiguous — out of v1 scope
    };
    match simple.as_ref() {
        SimpleModificationInner::Mass(mass) => Some(Mod::Mass(mass.value as f32)),
        SimpleModificationInner::Database { id, .. } => {
            if id.ontology == Ontology::Unimod {
                Some(Mod::Unimod(id.id? as u16))
            } else {
                None
            }
        }
        _ => None,
    }
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
/// `)` that closes it becomes `]`. Any other paren — notably `)` inside a bracket
/// mod name like `[Carbamidomethyl (C)]` — is left untouched.
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
                    // Malformed (no closer) — copy the remainder verbatim.
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

/// Coerce DIA-NN / short-form modified-sequence strings into rustyms-parseable
/// ProForma. Strips `_..._` wrapping used by DIA-NN, converts DIA-NN's
/// parenthesised mods (`C(UniMod:4)`) to ProForma brackets (`C[UNIMOD:4]`), and
/// normalizes UNIMOD tag casing (`[UniMod:`, `[Unimod:`, `[U:` → `[UNIMOD:`).
/// A leading (N-terminal) mod is rewritten to `[UNIMOD:n]-SEQ` as ProForma
/// requires. Pass-through for plain sequences.
///
/// Off the hot path — allocates on every replacement, which is fine at load.
pub fn normalize_to_proforma(raw: &str) -> String {
    let trimmed = raw.trim_matches('_');
    // Fast path: plain-AA sequences (no mod tags) skip the rewrite chain.
    if !trimmed.contains('[') && !trimmed.contains('(') {
        return trimmed.to_owned();
    }

    // DIA-NN writes mods in parentheses, e.g. `C(UniMod:4)`; ProForma needs
    // brackets, e.g. `C[UNIMOD:4]`. Convert each opener and ONLY its matching
    // `)` (see `convert_paren_unimod`) — never a blanket `)` replace, which would
    // corrupt parens inside a bracket mod name, e.g. `C[Carbamidomethyl (C)]`.
    let s = convert_paren_unimod(trimmed);

    // Normalize any pre-existing bracket casing likewise.
    let mut s = replace_ascii_ci(&s, "[unimod:", "[UNIMOD:");
    s = replace_ascii_ci(&s, "[u:", "[UNIMOD:");

    // A mod at the very start is N-terminal: ProForma wants `[UNIMOD:n]-SEQ`.
    if s.starts_with('[')
        && let Some(close) = s.find(']')
        && s[close + 1..].chars().next().is_some_and(|c| c != '-')
    {
        s.insert(close + 1, '-');
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use smallvec::smallvec;

    fn ci(ch: u8) -> usize {
        CANONICAL_AA_LETTERS.iter().position(|&x| x == ch).unwrap()
    }

    #[test]
    fn aa_counts_peptidek() {
        let seq = ParsedSequence {
            residues: "PEPTIDEK".bytes().map(AminoAcid::from_ascii).collect(),
            mods: smallvec![],
        };
        let c = seq.aa_counts();
        assert_eq!(c[ci(b'P')], 2.0);
        assert_eq!(c[ci(b'E')], 2.0);
        assert_eq!(c[ci(b'T')], 1.0);
        assert_eq!(c[ci(b'I')], 1.0);
        assert_eq!(c[ci(b'D')], 1.0);
        assert_eq!(c[ci(b'K')], 1.0);
        assert_eq!(c[ci(b'A')], 0.0);
        assert_eq!(c.iter().sum::<f64>(), 8.0);
    }

    #[test]
    fn aa_counts_skip_unknown() {
        // 'X' maps to a non-canonical alpha slot (tmp[23]); never projected to out[].
        let seq = ParsedSequence {
            residues: "AXA".bytes().map(AminoAcid::from_ascii).collect(),
            mods: smallvec![],
        };
        let c = seq.aa_counts();
        assert_eq!(c[ci(b'A')], 2.0);
        assert_eq!(c.iter().sum::<f64>(), 2.0);
    }

    #[test]
    fn aa_encoding_alpha_offset() {
        assert_eq!(AminoAcid::from_ascii(b'A').0, 0);
        assert_eq!(AminoAcid::from_ascii(b'Y').0, 24);
        assert_eq!(AminoAcid::from_ascii(b'Z').0, 25);
        assert_eq!(AminoAcid::from_ascii(b'!').0, u8::MAX);
    }

    #[test]
    fn normalize_strips_underscores() {
        assert_eq!(normalize_to_proforma("_PEPTIDEK_"), "PEPTIDEK");
    }

    #[test]
    fn normalize_diann_unimod_case() {
        assert_eq!(
            normalize_to_proforma("_LSHPGC[UniMod:4]K_"),
            "LSHPGC[UNIMOD:4]K"
        );
        assert_eq!(
            normalize_to_proforma("_C[Unimod:4]TVPGHK_"),
            "C[UNIMOD:4]TVPGHK"
        );
    }

    #[test]
    fn normalize_short_u_form() {
        assert_eq!(
            normalize_to_proforma("PEPTC[U:4]IDEK"),
            "PEPTC[UNIMOD:4]IDEK"
        );
    }

    #[test]
    fn normalize_diann_paren_unimod() {
        // DIA-NN parenthesised mods -> ProForma brackets (internal residue mod).
        assert_eq!(
            normalize_to_proforma("AAC(UniMod:4)DEK"),
            "AAC[UNIMOD:4]DEK"
        );
        // Case-insensitive on the tag.
        assert_eq!(
            normalize_to_proforma("AAC(unimod:4)DEK"),
            "AAC[UNIMOD:4]DEK"
        );
        assert_eq!(
            normalize_to_proforma("AAC(UNIMOD:4)DEK"),
            "AAC[UNIMOD:4]DEK"
        );
        // Multiple mods in one peptide.
        assert_eq!(
            normalize_to_proforma("AAC(UniMod:4)M(UniMod:35)K"),
            "AAC[UNIMOD:4]M[UNIMOD:35]K"
        );
    }

    #[test]
    fn normalize_diann_paren_nterm_gets_dash() {
        // A leading (N-terminal) mod must become `[UNIMOD:n]-SEQ`.
        assert_eq!(
            normalize_to_proforma("(UniMod:1)AACDEK"),
            "[UNIMOD:1]-AACDEK"
        );
    }

    #[test]
    fn parse_diann_paren_unimod_roundtrips() {
        // The end-to-end path the load uses: normalize then parse. Parenthesised
        // DIA-NN mods must parse (else the parse gate disables sequence features).
        let norm = normalize_to_proforma("AAC(UniMod:4)DEK");
        let p = parse_sequence(&norm).expect("paren UniMod must parse");
        assert_eq!(p.residues.len(), 6);
        assert_eq!(p.mods.len(), 1);
        assert_eq!(p.mods[0].kind, Mod::Unimod(4));
    }

    #[test]
    fn normalize_mixed_paren_unimod_and_bracket_paren_untouched() {
        // A paren UNIMOD tag AND a bracket mod whose name contains `)`: only the
        // UNIMOD `)` converts; the `(M)` inside the bracket name stays intact.
        assert_eq!(
            normalize_to_proforma("AAC(UniMod:4)M[Oxidation (M)]K"),
            "AAC[UNIMOD:4]M[Oxidation (M)]K"
        );
    }

    #[test]
    fn normalize_spectronaut_parens_in_brackets_untouched() {
        // Spectronaut writes mod names with parens INSIDE brackets. The DIA-NN
        // paren->bracket conversion must not touch these (no `(unimod:` opener).
        assert_eq!(
            normalize_to_proforma("_C[Carbamidomethyl (C)]PEPK_"),
            "C[Carbamidomethyl (C)]PEPK"
        );
    }

    #[test]
    fn normalize_skyline_nterm_mass_gets_dash() {
        // Skyline N-terminal mass mod -> ProForma `[+42]-SEQ`.
        assert_eq!(normalize_to_proforma("[+42]AACDEK"), "[+42]-AACDEK");
    }

    #[test]
    fn normalize_mass_shift_unchanged() {
        assert_eq!(
            normalize_to_proforma("PEPTM[+15.995]IDEK"),
            "PEPTM[+15.995]IDEK"
        );
    }

    #[test]
    fn normalize_plain_unchanged() {
        assert_eq!(normalize_to_proforma("PEPTIDEK"), "PEPTIDEK");
    }

    #[test]
    fn parse_plain() {
        let p = parse_sequence("PEPTIDEK").expect("parse");
        assert_eq!(p.residues.len(), 8);
        assert_eq!(p.mods.len(), 0);
        assert_eq!(p.residues[0], AminoAcid::from_ascii(b'P'));
        assert_eq!(p.residues[7], AminoAcid::from_ascii(b'K'));
    }

    #[test]
    fn parse_with_unimod() {
        // PEPTC[UNIMOD:4]IDEK = 9 residues: P-E-P-T-C-I-D-E-K, C at index 4
        let p = parse_sequence("PEPTC[UNIMOD:4]IDEK").expect("parse");
        assert_eq!(p.residues.len(), 9);
        assert_eq!(p.mods.len(), 1);
        assert_eq!(p.mods[0].pos, 4);
        assert_eq!(p.mods[0].kind, Mod::Unimod(4));
    }

    #[test]
    fn parse_with_mass_shift() {
        // PEPTM[+15.995]IDEK = 9 residues: P-E-P-T-M-I-D-E-K, M at index 4
        let p = parse_sequence("PEPTM[+15.995]IDEK").expect("parse");
        assert_eq!(p.residues.len(), 9);
        assert_eq!(p.mods.len(), 1);
        assert_eq!(p.mods[0].pos, 4);
        match p.mods[0].kind {
            Mod::Mass(m) => assert!((m - 15.995).abs() < 1e-3),
            _ => panic!("expected Mass variant"),
        }
    }

    #[test]
    fn parse_rejects_garbage() {
        assert!(parse_sequence("not a peptide!!!").is_none());
    }

    #[test]
    fn parse_n_term_mod() {
        // Acetyl (UNIMOD:1) on n-term
        let p = parse_sequence("[Acetyl]-PEPTIDEK").expect("parse");
        assert_eq!(p.residues.len(), 8);
        assert_eq!(p.mods.len(), 1);
        assert_eq!(p.mods[0].pos, POS_N_TERM);
        assert_eq!(p.mods[0].kind, Mod::Unimod(1));
    }

    #[test]
    fn fast_path_takes_recognized_grammar() {
        // Bare and numeric-UNIMOD/mass inputs must be served by the fast path
        // (never reach rustyms), else there's no speedup.
        for s in [
            "PEPTIDEK",
            "AAC[UNIMOD:4]DEK",
            "AAC[UNIMOD:4]M[UNIMOD:35]K",
            "[UNIMOD:1]-PEPTIDEK",
            "M[+15.995]PEPTIDEK",
            "[+42]-PEPTIDEK",
        ] {
            assert!(
                parse_sequence_fast(s).is_some(),
                "fast path must accept {s:?}"
            );
        }
    }

    #[test]
    fn fast_path_defers_named_and_garbage() {
        // Named mods and unexpected bytes must defer (fast returns None) so the
        // rustyms authority decides validity + resolves the name.
        for s in [
            "[Acetyl]-PEPTIDEK",
            "C[Carbamidomethyl (C)]PEPK",
            "not a pep!",
        ] {
            assert!(
                parse_sequence_fast(s).is_none(),
                "fast path must defer {s:?}"
            );
        }
    }

    #[test]
    fn fast_matches_rustyms_on_recognized_grammar() {
        // Differential test: wherever the fast path claims an input, it must
        // produce the exact same ParsedSequence rustyms would. Guards against
        // the fast path silently diverging on residue counts or mod mapping.
        let corpus = [
            "PEPTIDEK",
            "AACDEK",
            "AAC[UNIMOD:4]DEK",
            "AAC[UNIMOD:4]M[UNIMOD:35]K",
            "[UNIMOD:1]-AACDEK",
            "M[UNIMOD:35]LEGNSPQGSNQGVK",
            "AAAGAAATHLEVAR",
        ];
        for s in corpus {
            if let Some(fast) = parse_sequence_fast(s) {
                let slow = parse_sequence_rustyms(s)
                    .unwrap_or_else(|| panic!("rustyms must also parse {s:?}"));
                assert_eq!(fast.residues, slow.residues, "residues mismatch for {s:?}");
                assert_eq!(
                    fast.mods.len(),
                    slow.mods.len(),
                    "n_mods mismatch for {s:?}"
                );
                assert_eq!(
                    fast.aa_counts(),
                    slow.aa_counts(),
                    "aa_counts mismatch for {s:?}"
                );
            }
        }
    }
}
