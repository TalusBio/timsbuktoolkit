use std::path::PathBuf;
use std::process::Command;

#[test]
fn build_library_writes_a_readable_library_and_sidecar() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("tiny.mzspeclib.txt.gz");
    let fasta = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/test_data/tiny.fasta");

    let status = Command::new(env!("CARGO_BIN_EXE_timsseek"))
        .args([
            "build-library",
            "--fasta",
            fasta.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
            "--no-decoys",
            "--no-fixed-mods",
            "--max-fragments",
            "4",
        ])
        .status()
        .unwrap();

    assert!(status.success());
    assert!(out.exists());
    let sidecar = PathBuf::from(format!("{}.config.json", out.display()));
    let provenance: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&sidecar).unwrap()).unwrap();
    assert_eq!(provenance["output"]["path"], out.display().to_string());
    let table = timsquery::serde::read_targets_with(
        &out,
        timsquery::models::capabilities::LoadPolicy::default(),
    )
    .expect("generated library reads back");
    let rows = match table {
        timsquery::serde::TargetTable::Mzpaf { geom, .. } => geom.n_rows(),
        timsquery::serde::TargetTable::Str { geom, .. } => geom.n_rows(),
    };
    assert!(rows > 0);

    let rebuilt = Command::new(env!("CARGO_BIN_EXE_timsseek"))
        .args([
            "build-library",
            "--fasta",
            fasta.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
            "--no-decoys",
            "--no-fixed-mods",
            "--max-fragments",
            "3",
            "--overwrite",
        ])
        .output()
        .unwrap();
    assert!(rebuilt.status.success());
    let stderr = String::from_utf8_lossy(&rebuilt.stderr);
    assert!(stderr.contains("built with different settings"), "{stderr}");
    assert!(stderr.contains("fragments.max_fragments"), "{stderr}");
}

#[test]
fn peptide_tsv_preserves_modified_charges_and_supplied_decoy_groups() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("peptides.tsv");
    std::fs::write(
        &input,
        concat!(
            "proforma\tprotein_ids\tdecoy\tdecoy_group\n",
            "PEC[UNIMOD:4]TIDEK/2\tP1;P2\tfalse\tcam\n",
            "PECTIDEK/3\tP1\tfalse\tplain\n",
            "PEK[UNIMOD:259]TIDEK/2\tPRTC\tfalse\theavy\n",
            "PECTDIEK/3\tP1\ttrue\tplain\n",
        ),
    )
    .unwrap();
    let out = dir.path().join("peptides.mzspeclib.txt");
    let result = Command::new(env!("CARGO_BIN_EXE_timsseek"))
        .args([
            "build-library",
            "--peptides",
            input.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
            "--no-decoys",
            "--max-fragments",
            "4",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let text = std::fs::read_to_string(&out).unwrap();
    for peptide in [
        "PEC[UNIMOD:4]TIDEK/2",
        "PECTIDEK/3",
        "PEK[UNIMOD:259]TIDEK/2",
    ] {
        assert!(text.contains(peptide), "missing {peptide}");
    }
    assert_eq!(text.matches("other attribute value=plain").count(), 2);
    let table = timsquery::serde::read_targets_with(
        &out,
        timsquery::models::capabilities::LoadPolicy::default(),
    )
    .unwrap();
    let timsquery::serde::TargetTable::Mzpaf { geom, .. } = table else {
        panic!("predicted fragments should have mzPAF annotations");
    };
    assert_eq!(geom.n_rows(), 4);
    assert_eq!(geom.rows().filter(|r| geom.is_decoy(*r)).count(), 1);
    assert_eq!(
        geom.rows()
            .map(|r| geom.decoy_group_code(r))
            .collect::<std::collections::HashSet<_>>()
            .len(),
        3
    );
    assert!(geom.rows().any(|r| geom.charge(r) == 3));
}

#[test]
fn peptide_tsv_generates_seeded_decoys_and_assigns_groups() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("peptides.tsv");
    std::fs::write(
        &input,
        "proforma\tprotein_ids\nPEC[UNIMOD:4]TIDEK/2\tP1\nPECTIDEK/2\tP1\n",
    )
    .unwrap();
    let out = dir.path().join("peptides.mzspeclib.txt");
    let result = Command::new(env!("CARGO_BIN_EXE_timsseek"))
        .args([
            "build-library",
            "--peptides",
            input.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
            "--decoy-method",
            "shuffle",
            "--decoy-seed",
            "42",
            "--max-fragments",
            "4",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let table = timsquery::serde::read_targets_with(
        &out,
        timsquery::models::capabilities::LoadPolicy::default(),
    )
    .unwrap();
    let timsquery::serde::TargetTable::Mzpaf { geom, .. } = table else {
        panic!("predicted fragments should have mzPAF annotations");
    };
    assert_eq!(geom.n_rows(), 4);
    assert_eq!(geom.rows().filter(|r| geom.is_decoy(*r)).count(), 2);
    assert_eq!(
        geom.rows()
            .map(|r| geom.decoy_group_code(r))
            .collect::<std::collections::HashSet<_>>()
            .len(),
        2
    );
    let sidecar: serde_json::Value =
        serde_json::from_slice(&std::fs::read(format!("{}.config.json", out.display())).unwrap())
            .unwrap();
    assert_eq!(sidecar["decoys"]["method"], "shuffle");
    assert_eq!(sidecar["decoys"]["seed"], 42);
}

#[test]
fn peptide_tsv_retries_pseudo_reverse_after_a_target_collision() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("peptides.tsv");
    std::fs::write(
        &input,
        "proforma\tprotein_ids\nPEPTIDEK/2\tP1\nPEDITPEK/3\tP2\n",
    )
    .unwrap();
    let out = dir.path().join("peptides.mzspeclib.txt");
    let result = Command::new(env!("CARGO_BIN_EXE_timsseek"))
        .args([
            "build-library",
            "--peptides",
            input.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
            "--decoy-method",
            "pseudo-reverse",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(!String::from_utf8_lossy(&result.stderr).contains("skipped"));
    let table = timsquery::serde::read_targets_with(
        &out,
        timsquery::models::capabilities::LoadPolicy::default(),
    )
    .unwrap();
    let timsquery::serde::TargetTable::Mzpaf { geom, .. } = table else {
        panic!("predicted fragments should have mzPAF annotations");
    };
    assert_eq!(geom.n_rows(), 4);
    assert_eq!(geom.rows().filter(|r| geom.is_decoy(*r)).count(), 2);
    assert!(
        std::fs::read_to_string(&out)
            .unwrap()
            .contains("PEPITDEK/2")
    );
}

#[test]
fn peptide_tsv_keeps_unpaired_targets_when_every_reversal_collides() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("peptides.tsv");
    std::fs::write(&input, "proforma\tprotein_ids\nPEEEEEEK/2\tP1\n").unwrap();
    let out = dir.path().join("peptides.mzspeclib.txt");
    let result = Command::new(env!("CARGO_BIN_EXE_timsseek"))
        .args([
            "build-library",
            "--peptides",
            input.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
            "--decoy-method",
            "pseudo-reverse",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(String::from_utf8_lossy(&result.stderr).contains("skipped 1 generated decoys"));
    let table = timsquery::serde::read_targets_with(
        &out,
        timsquery::models::capabilities::LoadPolicy::default(),
    )
    .unwrap();
    let timsquery::serde::TargetTable::Mzpaf { geom, .. } = table else {
        panic!("predicted fragments should have mzPAF annotations");
    };
    assert_eq!(geom.n_rows(), 1);
    assert_eq!(geom.rows().filter(|r| geom.is_decoy(*r)).count(), 0);
}
