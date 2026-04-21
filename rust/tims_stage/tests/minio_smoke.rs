mod common;
use common::minio;

#[test]
#[cfg_attr(not(feature = "aws"), ignore = "requires aws feature")]
fn stage_s3_prefix_against_minio() {
    let Some(endpoint) = minio::minio_endpoint() else {
        eprintln!("MINIO_TEST_ENDPOINT not set — skipping");
        return;
    };
    unsafe {
        std::env::set_var("AWS_ENDPOINT_URL", &endpoint);
        std::env::set_var("AWS_S3_FORCE_PATH_STYLE", "true");
    }

    // Requires MinIO bucket pre-seeded with `smoke/sample.d/` prefix
    // containing `analysis.tdf` + `analysis.tdf_bin`. CI seeding is a
    // separate concern — this test does NOT create fixtures in MinIO.
    let bucket = minio::minio_bucket();
    let uri = format!("s3://{bucket}/smoke/sample.d/");

    let backend = tims_stage::PerRunTempdir::new(tims_stage::StagingConfig::default()).unwrap();
    let resolved = tims_stage::resolve(&uri).unwrap();
    let spec = match resolved {
        tims_stage::Resolved::Stageable { spec } => spec,
        other => panic!("expected Stageable, got {other:?}"),
    };
    let staged =
        <tims_stage::PerRunTempdir as tims_stage::StagingBackend>::stage(&backend, &spec).unwrap();
    assert!(staged.as_ref().join("analysis.tdf").exists());
    assert!(staged.as_ref().join("analysis.tdf_bin").exists());
}
