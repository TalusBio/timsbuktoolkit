use tempfile::TempDir;
use tims_stage::{
    PerRunTempdir,
    SourceSpec,
    StagingBackend,
    StagingConfig,
};
use timscentroid::StorageLocation;

#[test]
fn stage_s3_prefix_on_local_fs() {
    let root = TempDir::new().unwrap();
    let sample_dir = root.path().join("sample.d");
    std::fs::create_dir(&sample_dir).unwrap();
    std::fs::write(sample_dir.join("analysis.tdf"), b"tdf").unwrap();
    std::fs::write(sample_dir.join("analysis.tdf_bin"), b"tdfbin").unwrap();

    let loc = StorageLocation::from_path(root.path());
    let backend = PerRunTempdir::new(StagingConfig::default()).unwrap();
    let spec = SourceSpec::S3Prefix {
        loc,
        prefix: "sample.d".to_string(),
    };
    let staged = backend.stage(&spec).unwrap();

    assert_eq!(
        std::fs::read(staged.as_ref().join("analysis.tdf")).unwrap(),
        b"tdf"
    );
    assert_eq!(
        std::fs::read(staged.as_ref().join("analysis.tdf_bin")).unwrap(),
        b"tdfbin"
    );
}

#[test]
fn stage_s3_prefix_errors_on_missing_files() {
    let root = TempDir::new().unwrap();
    let sample_dir = root.path().join("sample.d");
    std::fs::create_dir(&sample_dir).unwrap();
    std::fs::write(sample_dir.join("analysis.tdf"), b"tdf").unwrap();
    // missing analysis.tdf_bin

    let loc = StorageLocation::from_path(root.path());
    let backend = PerRunTempdir::new(StagingConfig::default()).unwrap();
    let spec = SourceSpec::S3Prefix {
        loc,
        prefix: "sample.d".to_string(),
    };
    let err = match backend.stage(&spec) {
        Ok(_) => panic!("expected error, got Ok"),
        Err(e) => e,
    };
    let msg = format!("{err}");
    assert!(msg.contains("analysis.tdf_bin"), "{msg}");
}
