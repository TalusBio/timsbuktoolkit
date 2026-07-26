use std::path::PathBuf;
use tempfile::TempDir;
use tims_stage::{
    PerRunTempdir,
    StagingConfig,
};
use timscentroid::IndexingCentroidingConfig;
use timsquery::load_index;

fn fixture_dotd() -> PathBuf {
    // Adjust this path to a real fixture if available. The tests are
    // #[ignore]'d so CI does not require the fixture.
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../test_data/sample.d")
}

#[test]
#[ignore = "requires local .d fixture"]
fn load_index_from_local_dotd() {
    let backend = PerRunTempdir::new(StagingConfig::default()).unwrap();
    let uri = fixture_dotd().to_string_lossy().to_string();
    let (_idx, _) = load_index(&uri, &backend, false, IndexingCentroidingConfig::default())
        .expect("load_index should succeed");
}

#[test]
#[ignore = "requires local .d fixture"]
fn sidecar_roundtrip() {
    let backend = PerRunTempdir::new(StagingConfig::default()).unwrap();
    let dotd = fixture_dotd();
    let work = TempDir::new().unwrap();
    let dotd_copy = work.path().join("sample.d");
    copy_dir_recursive(&dotd, &dotd_copy).unwrap();
    let uri = dotd_copy.to_string_lossy().to_string();
    let (_idx1, _) =
        load_index(&uri, &backend, true, IndexingCentroidingConfig::default()).unwrap();
    assert!(
        work.path().join("sample.d.idx/metadata.json").is_file(),
        "sidecar .idx should contain metadata.json after build"
    );
    let (_idx2, _) =
        load_index(&uri, &backend, false, IndexingCentroidingConfig::default()).unwrap();
}

fn copy_dir_recursive(src: &std::path::Path, dst: &std::path::Path) -> std::io::Result<()> {
    std::fs::create_dir_all(dst)?;
    for entry in std::fs::read_dir(src)? {
        let entry = entry?;
        let ft = entry.file_type()?;
        let dst_child = dst.join(entry.file_name());
        if ft.is_dir() {
            copy_dir_recursive(&entry.path(), &dst_child)?;
        } else {
            std::fs::copy(entry.path(), &dst_child)?;
        }
    }
    Ok(())
}
