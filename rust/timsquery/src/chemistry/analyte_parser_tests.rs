//! Reader-parser regressions migrated from timsseek's retired sequence parser.
use super::*;

#[test]
fn library_spellings_resolve_to_the_expected_structure() {
    for (input, canonical) in [
        ("PEPTIDEK", "PEPTIDEK"),
        ("AAC(UniMod:4)DEK", "AAC[UNIMOD:4]DEK"),
        ("FAC[U:Carbamidomethyl]HSASLTVR/3", "FAC[UNIMOD:4]HSASLTVR"),
        ("FAC[U:Carbamidomethyl]HSASLTVR", "FAC[UNIMOD:4]HSASLTVR"),
        ("AVC[U:Carbamidomethyl]ASFSLTHR/3", "AVC[UNIMOD:4]ASFSLTHR"),
        ("PEPTC[U:Carbamidomethyl]IDEK", "PEPTC[UNIMOD:4]IDEK"),
        ("PEPTC[Carbamidomethyl]IDEK", "PEPTC[UNIMOD:4]IDEK"),
        ("PEPTC[U:4]IDEK", "PEPTC[UNIMOD:4]IDEK"),
        ("[Acetyl]-PEPTIDEK", "[UNIMOD:1]-PEPTIDEK"),
        ("PEPTM[Oxidation]IDEK", "PEPTM[UNIMOD:35]IDEK"),
    ] {
        let parsed = Analyte::from_sequence(input);
        assert!(
            parsed
                .peptide
                .known()
                .unwrap()
                .modifications
                .known()
                .is_some(),
            "{input}"
        );
        assert_eq!(parsed, Analyte::from_sequence(canonical), "{input}");
    }
}

#[test]
fn invalid_names_and_unreadable_sequences_remain_unresolved() {
    for input in [
        "not a peptide!!!",
        "PEPTC[UNIMOD:Carbamidomethyl]IDEK",
        "PEPTC[U:CAM]IDEK",
        "PEPTC[CAM]IDEK",
        "PEPTC[UNIMOD:CAM]IDEK",
        "PEPTC[UNIMOD:999999]IDEK",
        "PEPTN[GNO:G59626AS]IDEK",
    ] {
        assert!(
            matches!(
                Analyte::from_sequence(input).peptide,
                Property::Unresolved { .. }
            ),
            "{input}"
        );
    }
}

#[test]
fn stored_modifications_retain_sites_and_mass() {
    let parsed = Analyte::from_sequence("[Acetyl]-PEPTM[+15.995]IDEK-[UNIMOD:2]");
    let peptide = parsed.peptide.known().unwrap();
    assert_eq!(peptide.residues, "PEPTMIDEK");
    let mods = peptide.modifications.known().unwrap();
    assert_eq!(mods.len(), 3);
    assert_eq!(mods[0].site, ModificationSite::NTerm);
    assert_eq!(mods[0].definition.kind(), &ModificationKind::Unimod(1));
    assert_eq!(mods[1].site, ModificationSite::Residue(4));
    assert_eq!(mods[1].definition.kind(), &ModificationKind::Mass(15.995));
    assert_eq!(mods[2].site, ModificationSite::CTerm);
}

#[test]
fn fast_parser_matches_ontology_parser_on_explicit_grammar() {
    for sequence in [
        "PEPTIDEK",
        "AACDEK",
        "AAC[UNIMOD:4]DEK",
        "AAC[UNIMOD:4]M[UNIMOD:35]K",
        "[UNIMOD:1]-AACDEK",
        "M[UNIMOD:35]LEGNSPQGSNQGVK",
        "AAAGAAATHLEVAR",
        "M[+15.995]PEPTIDEK",
        "[+42]-PEPTIDEK",
    ] {
        let fast = parse_simple(sequence).expect(sequence);
        let (parsed, _) = mzcore::sequence::Peptidoform::pro_forma(sequence, ontologies()).unwrap();
        let slow = Analyte::from_peptidoform(&parsed);
        assert_eq!(&fast, slow.peptide.known().unwrap(), "{sequence}");
    }
    for sequence in [
        "[Acetyl]-PEPTIDEK",
        "C[Carbamidomethyl (C)]PEPK",
        "not a pep!",
    ] {
        assert!(parse_simple(sequence).is_none(), "{sequence}");
    }
}

#[test]
fn other_ontologies_preserve_facts_without_claiming_unimod_features() {
    for sequence in [
        "PEPTK[MOD:00046]IDEK",
        "PEPTK[XLMOD:02001]IDEK",
        "PEPTK[RESID:AA0038]IDEK",
        "PEPTN[Glycan:HexNAc]IDEK",
    ] {
        let parsed = Analyte::from_sequence(sequence);
        let peptide = parsed.peptide.known().expect(sequence);
        let mods = peptide.modifications.known().expect(sequence);
        assert_eq!(mods.len(), 1);
        assert_eq!(
            mods[0].definition.kind(),
            &ModificationKind::Other,
            "{sequence}"
        );
    }
}
