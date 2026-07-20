mod common;
use common::minio;

/// Round-trips a synthetic `.d` prefix through MinIO.
///
/// Marked `#[ignore]` so it never runs in the default suite. CI and local
/// dev must invoke explicitly:
///
///   MINIO_TEST_ENDPOINT=http://localhost:9000 \
///   MINIO_TEST_BUCKET=tims-stage-ci \
///   AWS_ACCESS_KEY_ID=... AWS_SECRET_ACCESS_KEY=... \
///   cargo test -p tims_stage --features aws --test minio_smoke -- --ignored
///
/// Fixtures are seeded by the test itself (two tiny objects) so no external
/// bucket state is assumed. A missing `MINIO_TEST_ENDPOINT` is a hard failure
/// when the test runs — silent passes would mask CI misconfiguration.
#[test]
#[ignore = "requires MinIO endpoint + aws feature; run explicitly with --ignored"]
fn stage_manifest_against_minio() {
    use timscentroid::reader::Manifest;

    let endpoint = minio::minio_endpoint()
        .expect("MINIO_TEST_ENDPOINT must be set when running this ignored test");
    unsafe {
        std::env::set_var("AWS_ENDPOINT_URL", &endpoint);
        std::env::set_var("AWS_S3_FORCE_PATH_STYLE", "true");
    }

    let bucket = minio::minio_bucket();
    // Unique prefix per run so parallel CI jobs / repeated invocations don't collide.
    let run_id = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let entry = format!("s3://{bucket}/smoke/{run_id}/sample.d");

    let tmp = tempfile::TempDir::new().unwrap();
    let tdf = tmp.path().join("analysis.tdf");
    let tdf_bin = tmp.path().join("analysis.tdf_bin");
    std::fs::write(&tdf, b"fake-tdf-payload").unwrap();
    std::fs::write(&tdf_bin, b"fake-tdf-bin-payload").unwrap();

    tims_stage::upload_file(&tdf, &format!("{entry}/analysis.tdf")).unwrap();
    tims_stage::upload_file(&tdf_bin, &format!("{entry}/analysis.tdf_bin")).unwrap();

    // Manifest-driven staging fetches exactly the declared members, by real
    // name, into `<tempdir>/sample.d/`.
    let manifest = Manifest {
        entry: entry.parse().unwrap(),
        required: vec![
            format!("{entry}/analysis.tdf").parse().unwrap(),
            format!("{entry}/analysis.tdf_bin").parse().unwrap(),
        ],
        optional: vec![],
    };

    let backend = tims_stage::PerRunTempdir::new(tims_stage::StagingConfig::default()).unwrap();
    let staged = tims_stage::stage_manifest(&backend, &manifest).unwrap();
    let entry_dir = staged.source().entry_path();

    let got_tdf = std::fs::read(entry_dir.join("analysis.tdf")).unwrap();
    let got_bin = std::fs::read(entry_dir.join("analysis.tdf_bin")).unwrap();
    assert_eq!(got_tdf, b"fake-tdf-payload");
    assert_eq!(got_bin, b"fake-tdf-bin-payload");
}
