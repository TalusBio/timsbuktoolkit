use std::path::PathBuf;

fn test_fasta_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("test_data")
        .join("tiny.fasta")
}

#[test]
fn test_pipeline_digestion_and_dedup() {
    let fasta_path = test_fasta_path();
    assert!(fasta_path.exists());

    use timsseek::digest::digestion::*;
    use timsseek::protein::fasta::ProteinSequenceCollection;

    let proteins = ProteinSequenceCollection::from_fasta_file(&fasta_path).unwrap();
    assert_eq!(proteins.sequences.len(), 2);

    let params = DigestionParameters {
        min_length: 6,
        max_length: 25,
        pattern: DigestionPattern::trypsin(),
        digestion_end: DigestionEnd::CTerm,
        max_missed_cleavages: 1,
    };

    let all_digests: Vec<_> = proteins
        .sequences
        .iter()
        .flat_map(|p| params.digest(p.sequence.clone()))
        .collect();
    assert!(all_digests.len() > 0);

    use speclib_build_cli::dedup::PeptideDedup;
    let deduped = PeptideDedup::dedup(all_digests.clone(), 100);
    assert!(deduped.len() <= all_digests.len());
    assert!(deduped.len() > 0);
}

#[test]
fn test_mod_application_chain() {
    use speclib_build_cli::mods::*;

    let fixed = vec![Modification::parse("C[U:4]").unwrap()];
    let var = vec![Modification::parse("M[U:35]").unwrap()];

    let seq = "PEPTMCIDECK";
    let modified = apply_fixed_mods(seq, &fixed);
    assert_eq!(modified, "PEPTMC[U:4]IDEC[U:4]K");

    let expanded = expand_variable_mods(&modified, &var, 1);
    assert!(expanded.len() == 2); // unmodified + M oxidized
    assert!(expanded.contains(&modified));
}

#[test]
fn test_decoy_roundtrip() {
    use speclib_build_cli::decoys::*;

    let original = "PEPTIDEK";
    let reversed = generate_decoy(original, DecoyMode::Reverse);
    assert_ne!(reversed, original);
    assert_eq!(reversed.len(), original.len());
    assert_eq!(reversed.chars().next(), original.chars().next());
    assert_eq!(reversed.chars().last(), original.chars().last());

    let edge_mutated = generate_decoy(original, DecoyMode::EdgeMutate);
    assert_ne!(edge_mutated, original);
    assert_eq!(edge_mutated.len(), original.len());
}
