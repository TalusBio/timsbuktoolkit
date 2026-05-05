//! MinIO test fixture (requires `MINIO_TEST_ENDPOINT` env var).
//! Returns None when unset; callers should skip the test gracefully.

pub fn minio_endpoint() -> Option<String> {
    std::env::var("MINIO_TEST_ENDPOINT").ok()
}

pub fn minio_bucket() -> String {
    std::env::var("MINIO_TEST_BUCKET").unwrap_or_else(|_| "tims-stage-ci".to_string())
}
