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
        timsquery::serde::TargetTable::Str { geom } => geom.n_rows(),
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
