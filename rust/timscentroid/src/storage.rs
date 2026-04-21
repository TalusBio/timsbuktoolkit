//! Storage abstraction for local and cloud object stores
//!
//! This module provides a unified interface for accessing files from local
//! filesystem or cloud object storage (S3, GCS, Azure).
//!
//! # Features
//!
//! - **Unified API**: Same code works for local files and cloud storage
//! - **Lazy Runtime**: Tokio runtime created only when needed
//! - **Feature Flags**: Enable cloud providers via Cargo features
//! - **Sync API**: All async operations hidden behind synchronous interface
//!
//! # Examples
//!
//! ## Local Filesystem
//!
//! ```no_run
//! use timscentroid::{StorageLocation, StorageProvider};
//!
//! // Create from local path
//! let location = StorageLocation::from_path("/path/to/data");
//! let storage = StorageProvider::new(location)?;
//!
//! // Read and write files
//! let data = storage.read_bytes("file.txt")?;
//! storage.write_bytes("output.txt", vec![1, 2, 3])?;
//! # Ok::<(), timscentroid::serialization::SerializationError>(())
//! ```
//!
//! ## Cloud Storage (S3)
//!
//! ```no_run
//! use timscentroid::{StorageLocation, StorageProvider};
//!
//! // Requires "aws" feature flag
//! let location = StorageLocation::from_url("s3://my-bucket/prefix")?;
//! let storage = StorageProvider::new(location)?;
//!
//! // Same API as local filesystem
//! let data = storage.read_bytes("data.parquet")?;
//! # Ok::<(), timscentroid::serialization::SerializationError>(())
//! ```
//!
//! # Authentication
//!
//! Cloud providers use default credential chains:
//! - **AWS**: `~/.aws/credentials`, IAM roles, or environment variables
//! - **GCP**: `GOOGLE_APPLICATION_CREDENTIALS` environment variable
//! - **Azure**: `AZURE_STORAGE_ACCOUNT` and `AZURE_STORAGE_KEY` environment variables

use bytes::Bytes;
use object_store::ObjectStore;
use object_store::local::LocalFileSystem;
use object_store::path::Path as ObjectPath;
use once_cell::sync::Lazy;
use std::path::Path;
use std::sync::Arc;
use tokio::runtime::Runtime;

use crate::instrumentation::{
    InstrumentedStore,
    StorageMetrics,
};
use crate::serialization::SerializationError;
use tracing::{
    debug,
    instrument,
    trace,
};

/// Global tokio runtime for all async operations (created lazily)
///
/// By default, uses 8 worker threads to enable high concurrency for cloud storage I/O.
/// This is critical for performance when querying multiple parquet files concurrently
/// from S3 or other cloud providers.
///
/// Can be configured via the `TIMSCENTROID_WORKER_THREADS` environment variable.
pub(crate) static RUNTIME: Lazy<Runtime> = Lazy::new(|| {
    let worker_threads = std::env::var("TIMSCENTROID_WORKER_THREADS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(8); // Increased from 2 to 8 for better S3 concurrency

    tokio::runtime::Builder::new_multi_thread()
        .worker_threads(worker_threads)
        .enable_all()
        .build()
        .expect("Failed to create Tokio runtime")
});

/// Helper to run async code from sync context, handling nested runtime calls
fn block_on_or_in_place<F: std::future::Future>(future: F) -> F::Output {
    match tokio::runtime::Handle::try_current() {
        Ok(handle) => {
            // We're already in a runtime, use block_in_place to avoid panic
            tokio::task::block_in_place(|| handle.block_on(future))
        }
        Err(_) => {
            // Not in a runtime, use our global runtime
            RUNTIME.block_on(future)
        }
    }
}

/// Storage location - either a local path or a cloud URL
#[derive(Debug, Clone)]
pub enum StorageLocation {
    Local(std::path::PathBuf),
    Url(url::Url),
}

impl StorageLocation {
    /// Create from a local filesystem path
    pub fn from_path(path: impl AsRef<Path>) -> Self {
        Self::Local(path.as_ref().to_path_buf())
    }

    /// Create from a URL (s3://, gs://, az://, or file://)
    pub fn from_url(url: impl AsRef<str>) -> Result<Self, url::ParseError> {
        let parsed = url::Url::parse(url.as_ref())?;

        // Handle file:// URLs as local paths
        if parsed.scheme() == "file"
            && let Ok(path) = parsed.to_file_path()
        {
            return Ok(Self::Local(path));
        }

        Ok(Self::Url(parsed))
    }
}

/// Storage provider wrapping an ObjectStore
#[derive(Clone, Debug)]
pub struct StorageProvider {
    store: Arc<dyn ObjectStore>,
    is_local: bool,
    /// Path prefix to prepend to all file accesses (for cloud URLs with paths)
    prefix: String,
    /// Optional metrics if instrumentation is enabled
    metrics: Option<Arc<StorageMetrics>>,
}

impl StorageProvider {
    /// Create a new storage provider from a location
    pub fn new(location: StorageLocation) -> Result<Self, SerializationError> {
        let (store, is_local, prefix): (Arc<dyn ObjectStore>, bool, String) = match location {
            StorageLocation::Local(path) => {
                // Create directory if it doesn't exist
                std::fs::create_dir_all(&path)?;
                (
                    Arc::new(LocalFileSystem::new_with_prefix(path)?),
                    true,
                    String::new(),
                )
            }
            StorageLocation::Url(url) => {
                // Extract path prefix from URL (e.g., s3://bucket/prefix/path -> "prefix/path")
                let prefix = url.path().trim_start_matches('/').to_string();
                (block_on_or_in_place(parse_url(&url))?, false, prefix)
            }
        };

        Ok(Self {
            store,
            is_local,
            prefix,
            metrics: None,
        })
    }

    /// Non-creating constructor for read-only callers. Does NOT `create_dir_all`
    /// on a local path — use `new` when you want writes (which may implicitly
    /// create the parent directory).
    pub fn open(location: StorageLocation) -> Result<Self, SerializationError> {
        let (store, is_local, prefix): (Arc<dyn ObjectStore>, bool, String) = match location {
            StorageLocation::Local(path) => (
                Arc::new(LocalFileSystem::new_with_prefix(path)?),
                true,
                String::new(),
            ),
            StorageLocation::Url(url) => {
                let prefix = url.path().trim_start_matches('/').to_string();
                (block_on_or_in_place(parse_url(&url))?, false, prefix)
            }
        };
        Ok(Self {
            store,
            is_local,
            prefix,
            metrics: None,
        })
    }

    /// HEAD an object. Returns its metadata.
    pub fn head(&self, key: &str) -> Result<object_store::ObjectMeta, SerializationError> {
        let full = self.qualified_path(key);
        let store = self.store.clone();
        block_on_or_in_place(
            async move { store.head(&full).await.map_err(SerializationError::from) },
        )
    }

    /// Check whether an object exists. Returns `Ok(false)` for NotFound,
    /// `Err(..)` for other transport errors.
    pub fn exists(&self, key: &str) -> Result<bool, SerializationError> {
        match self.head(key) {
            Ok(_) => Ok(true),
            Err(SerializationError::ObjectStore(object_store::Error::NotFound { .. })) => Ok(false),
            Err(e) => Err(e),
        }
    }

    /// Fetch the whole object into memory.
    pub fn get_bytes(&self, key: &str) -> Result<Bytes, SerializationError> {
        let full = self.qualified_path(key);
        let store = self.store.clone();
        block_on_or_in_place(async move {
            let got = store.get(&full).await.map_err(SerializationError::from)?;
            let bytes = got.bytes().await.map_err(SerializationError::from)?;
            Ok(bytes)
        })
    }

    /// Fetch a specific byte range. Errors on short read — S3 returns 416 on
    /// out-of-bounds, but `LocalFileSystem` silently truncates to EOF, so we
    /// post-check the returned length against the requested length. No
    /// pre-HEAD; the tar walker issues many small range GETs and doubling
    /// the request count is expensive.
    pub fn range_get(
        &self,
        key: &str,
        range: std::ops::Range<u64>,
    ) -> Result<Bytes, SerializationError> {
        use object_store::GetRange;
        let full = self.qualified_path(key);
        let store = self.store.clone();
        let expected = range.end.saturating_sub(range.start);
        block_on_or_in_place(async move {
            let opts = object_store::GetOptions {
                range: Some(GetRange::Bounded(range)),
                ..Default::default()
            };
            let got = store
                .get_opts(&full, opts)
                .await
                .map_err(SerializationError::from)?;
            let bytes = got.bytes().await.map_err(SerializationError::from)?;
            if (bytes.len() as u64) < expected {
                return Err(SerializationError::Io(std::io::Error::new(
                    std::io::ErrorKind::UnexpectedEof,
                    format!(
                        "range_get short read: expected {} bytes, got {}",
                        expected,
                        bytes.len()
                    ),
                )));
            }
            Ok(bytes)
        })
    }

    /// Build a fully-qualified `ObjectPath` from a caller-provided key, applying
    /// the provider's URL-path prefix when present.
    fn qualified_path(&self, key: &str) -> ObjectPath {
        if self.prefix.is_empty() {
            ObjectPath::from(key)
        } else {
            ObjectPath::from(format!("{}/{}", self.prefix.trim_end_matches('/'), key))
        }
    }

    /// Enable instrumentation for this storage provider
    ///
    /// This wraps the underlying ObjectStore with an InstrumentedStore that tracks:
    /// - Number of GET/PUT/HEAD operations
    /// - Bytes transferred
    /// - Time spent in each operation
    ///
    /// Returns a new StorageProvider with instrumentation enabled.
    pub fn with_instrumentation(self, label: impl Into<String>) -> Self {
        let metrics = Arc::new(StorageMetrics::new());
        let instrumented = Arc::new(InstrumentedStore::new(
            self.store.clone(),
            metrics.clone(),
            label.into(),
        ));

        Self {
            store: instrumented,
            is_local: self.is_local,
            prefix: self.prefix,
            metrics: Some(metrics),
        }
    }

    /// Add fake latency to simulate network delays
    ///
    /// This method must be called AFTER `with_instrumentation()`.
    ///
    /// Creates a new metrics tracker for the outer layer (with latency), while
    /// the inner layer continues tracking without latency. This allows you to
    /// measure latency overhead: outer_time - inner_time = latency.
    ///
    /// # Panics
    /// Panics if called before `with_instrumentation()`.
    pub fn with_fake_latency(self, latency: std::time::Duration) -> Self {
        // Ensure instrumentation is enabled
        self.metrics
            .as_ref()
            .expect("with_fake_latency() must be called after with_instrumentation()");

        // Create NEW metrics for the outer wrapper (with latency)
        // The inner wrapper keeps its own metrics (without latency)
        // This is a feature: outer_time - inner_time = latency overhead
        let outer_metrics = Arc::new(StorageMetrics::new());

        let new_instrumented = Arc::new(
            InstrumentedStore::new(
                self.store.clone(), // Inner instrumented store (has its own metrics)
                outer_metrics.clone(),
                "with_latency",
            )
            .with_fake_latency(latency),
        );

        Self {
            store: new_instrumented,
            is_local: self.is_local,
            prefix: self.prefix,
            metrics: Some(outer_metrics), // Return outer metrics (includes latency)
        }
    }

    /// Get metrics if instrumentation is enabled
    pub fn metrics(&self) -> Option<&StorageMetrics> {
        self.metrics.as_ref().map(|m| m.as_ref())
    }

    /// Print metrics report if instrumentation is enabled
    pub fn print_metrics(&self, label: &str) {
        if let Some(metrics) = &self.metrics {
            metrics.snapshot().print_report(label);
        }
    }

    /// Check if this storage provider is backed by local filesystem
    ///
    /// Returns true for local paths, false for cloud storage (S3, GCS, Azure)
    pub fn is_local(&self) -> bool {
        self.is_local
    }

    /// Build full path by prepending prefix (for cloud storage)
    pub fn build_path(&self, path: &str) -> String {
        if self.prefix.is_empty() {
            path.to_string()
        } else {
            format!("{}/{}", self.prefix, path)
        }
    }

    /// Read a file as bytes (async version)
    #[instrument(skip(self), fields(path = %path))]
    pub async fn read_bytes_async(&self, path: &str) -> Result<Vec<u8>, SerializationError> {
        let full_path = self.build_path(path);
        trace!("Reading from full path: {}", full_path);
        let object_path = ObjectPath::from(full_path.as_str());
        trace!("Object path: {:?}", object_path);
        let result = match self.store.get(&object_path).await {
            Ok(res) => res,
            Err(e) => {
                // Categorize the error properly based on what it actually is
                let error_str = e.to_string();
                debug!("Error getting object: {:?}", e);

                // Check if it's an authentication/permission error
                let error_kind = if error_str.contains("ExpiredToken")
                    || error_str.contains("The provided token has expired")
                    || error_str.contains("Access Denied")
                    || error_str.contains("InvalidAccessKeyId")
                    || error_str.contains("SignatureDoesNotMatch")
                    || error_str.contains("Forbidden")
                    || error_str.contains("status: 401")
                    || error_str.contains("status: 403")
                {
                    std::io::ErrorKind::PermissionDenied
                } else if error_str.contains("status: 404") || error_str.contains("NotFound") {
                    std::io::ErrorKind::NotFound
                } else {
                    // For other errors, use Other to indicate an unexpected issue
                    std::io::ErrorKind::Other
                };

                return Err(SerializationError::Io(std::io::Error::new(
                    error_kind,
                    format!("Failed to read object at path {}: {}", full_path, e),
                )));
            }
        };
        let bytes = result.bytes().await?;
        Ok(bytes.to_vec())
    }

    /// Read a file as bytes (blocking version)
    pub fn read_bytes(&self, path: &str) -> Result<Vec<u8>, SerializationError> {
        block_on_or_in_place(self.read_bytes_async(path))
    }

    /// Read a file as string (async version)
    pub async fn read_to_string_async(&self, path: &str) -> Result<String, SerializationError> {
        let bytes = self.read_bytes_async(path).await?;
        String::from_utf8(bytes).map_err(|e| {
            SerializationError::Io(std::io::Error::new(std::io::ErrorKind::InvalidData, e))
        })
    }

    /// Read a file as string (blocking version)
    pub fn read_to_string(&self, path: &str) -> Result<String, SerializationError> {
        block_on_or_in_place(self.read_to_string_async(path))
    }

    /// Write bytes to a file
    pub fn write_bytes(&self, path: &str, data: Vec<u8>) -> Result<(), SerializationError> {
        use object_store::buffered::BufWriter;
        use tokio::io::AsyncWriteExt;
        let full = self.qualified_path(path);
        let store = self.store.clone();
        block_on_or_in_place(async move {
            let mut w = BufWriter::new(store, full);
            w.write_all(&data).await.map_err(|e| {
                SerializationError::Io(std::io::Error::other(format!("buffered write failed: {e}")))
            })?;
            w.shutdown().await.map_err(|e| {
                SerializationError::Io(std::io::Error::other(format!(
                    "buffered write shutdown failed: {e}"
                )))
            })?;
            Ok(())
        })
    }

    /// Stream an object into a local file, ticking the progress bar per chunk.
    pub fn get_to_file(
        &self,
        key: &str,
        dst: &Path,
        bar: &indicatif::ProgressBar,
    ) -> Result<(), SerializationError> {
        use futures::StreamExt;
        let full = self.qualified_path(key);
        let store = self.store.clone();
        let dst = dst.to_path_buf();
        block_on_or_in_place(async move {
            if let Some(parent) = dst.parent() {
                std::fs::create_dir_all(parent)?;
            }
            let mut file = std::fs::File::create(&dst)?;
            let got = store.get(&full).await.map_err(SerializationError::from)?;
            let mut stream = got.into_stream();
            while let Some(chunk) = stream.next().await {
                let chunk = chunk.map_err(SerializationError::from)?;
                use std::io::Write;
                file.write_all(&chunk)?;
                bar.inc(chunk.len() as u64);
            }
            Ok(())
        })
    }

    /// Get the underlying ObjectStore
    pub(crate) fn as_object_store(&self) -> Arc<dyn ObjectStore> {
        self.store.clone()
    }

    /// Ensure directory exists (no-op for object stores, but kept for API consistency)
    pub fn ensure_directory(&self, _path: &str) -> Result<(), SerializationError> {
        // Object stores don't have directories, so this is a no-op
        // Local filesystem handles this via put() creating parent "directories" as needed
        Ok(())
    }

    /// List up to `cap` entries under `prefix`. Errors with
    /// `SerializationError::PrefixCapExceeded` as soon as the listing yields
    /// more than `cap` items. This is a fail-fast mechanism to prevent callers
    /// from accidentally listing a whole bucket.
    pub fn list_capped(
        &self,
        prefix: &str,
        cap: usize,
    ) -> Result<Vec<object_store::ObjectMeta>, SerializationError> {
        use futures::StreamExt;
        let full = self.qualified_path(prefix);
        let store = self.store.clone();
        let prefix_display = full.to_string();
        block_on_or_in_place(async move {
            let mut stream = store.list(Some(&full));
            let mut out = Vec::with_capacity(cap);
            while let Some(item) = stream.next().await {
                let meta = item.map_err(SerializationError::from)?;
                out.push(meta);
                if out.len() > cap {
                    return Err(SerializationError::PrefixCapExceeded {
                        prefix: prefix_display.clone(),
                        cap,
                    });
                }
            }
            Ok(out)
        })
    }

    /// Read indexed peaks from a parquet file into SoA columns.
    ///
    /// Uses `ParquetObjectReader` with async streaming for both local and cloud storage.
    /// Validates mobility invariants (non-negative, non-NaN) at the boundary.
    #[instrument(skip(self), fields(path = %path))]
    pub(crate) fn read_parquet_peaks<T: crate::rt_mapping::RTIndex>(
        &self,
        path: &str,
    ) -> Result<crate::indexing::PeakColumns<T>, SerializationError> {
        use futures::stream::StreamExt;
        use parquet::arrow::ParquetRecordBatchStreamBuilder;
        use parquet::arrow::async_reader::ParquetObjectReader;

        let full_path = self.build_path(path);
        let object_path = ObjectPath::from(full_path.as_str());

        block_on_or_in_place(async {
            let reader = ParquetObjectReader::new(self.store.clone(), object_path);
            let builder = ParquetRecordBatchStreamBuilder::new(reader).await?;
            // Pre-size the SoA buffers using the file's row count so the
            // per-batch extends don't re-grow log2(n) times.
            let total_rows: usize = builder
                .metadata()
                .file_metadata()
                .num_rows()
                .try_into()
                .unwrap_or(0);
            let mut columns = crate::indexing::PeakColumns::<T>::with_capacity(total_rows);

            let mut stream = builder.build()?;
            while let Some(batch_result) = stream.next().await {
                let batch = batch_result?;
                crate::serialization::extend_soa_from_batch::<T>(&mut columns, &batch)?;
            }

            Ok(columns)
        })
    }
}

/// Parse a URL into an ObjectStore
///
/// Supports:
/// - s3://bucket/prefix (requires "aws" feature)
/// - gs://bucket/prefix (requires "gcp" feature)
/// - az://container/prefix or azure://container/prefix (requires "azure" feature)
async fn parse_url(url: &url::Url) -> Result<Arc<dyn ObjectStore>, SerializationError> {
    match url.scheme() {
        #[cfg(feature = "aws")]
        "s3" => {
            use aws_config::BehaviorVersion;
            use aws_credential_types::provider::ProvideCredentials;
            use object_store::aws::AmazonS3Builder;

            let bucket = url.host_str().ok_or_else(|| {
                SerializationError::Io(std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "Missing bucket in S3 URL",
                ))
            })?;
            debug!("Creating S3 ObjectStore for bucket: {}", bucket);

            // 1. Load the AWS configuration from the environment (handles Profile, MFA, SSO, etc.)
            let sdk_config = aws_config::load_defaults(BehaviorVersion::latest()).await;

            // 2. Extract the credentials from the resolved config
            //    (This executes the chain: Env Vars -> Profile -> Web Identity -> IMDS)
            let credentials_provider = sdk_config
                .credentials_provider()
                .expect("No credentials provider found");
            let credentials = credentials_provider
                .provide_credentials()
                .await
                .map_err(|e| {
                    SerializationError::Io(std::io::Error::new(std::io::ErrorKind::Other, e))
                })?;

            // 3. Initialize the builder using the resolved credentials
            let mut builder = AmazonS3Builder::new()
                .with_bucket_name(bucket)
                .with_region(
                    sdk_config
                        .region()
                        .map(|r| r.as_ref())
                        .unwrap_or("us-west-2"),
                )
                .with_access_key_id(credentials.access_key_id())
                .with_secret_access_key(credentials.secret_access_key());

            // 4. Important: Attach the session token if it exists (Critical for MFA/SSO)
            if let Some(token) = credentials.session_token() {
                builder = builder.with_token(token);
            }

            Ok(Arc::new(builder.build()?))
        }

        #[cfg(feature = "gcp")]
        "gs" => {
            use object_store::gcp::GoogleCloudStorageBuilder;

            let bucket = url.host_str().ok_or_else(|| {
                SerializationError::Io(std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "Missing bucket in GCS URL",
                ))
            })?;

            Ok(Arc::new(
                GoogleCloudStorageBuilder::from_env()
                    .with_bucket_name(bucket)
                    .build()?,
            ))
        }

        #[cfg(feature = "azure")]
        "az" | "azure" => {
            use object_store::azure::MicrosoftAzureBuilder;

            let container = url.host_str().ok_or_else(|| {
                SerializationError::Io(std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "Missing container in Azure URL",
                ))
            })?;

            Ok(Arc::new(
                MicrosoftAzureBuilder::from_env()
                    .with_container_name(container)
                    .build()?,
            ))
        }

        scheme =>
        {
            #[allow(unreachable_code)]
            Err(SerializationError::Io(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!(
                    "Unsupported URL scheme: '{}'. Available schemes: {}",
                    scheme,
                    {
                        let schemes = ["file"];
                        #[cfg(feature = "aws")]
                        let schemes = [schemes.as_slice(), &["s3"]].concat();
                        #[cfg(feature = "gcp")]
                        let schemes = [schemes.as_slice(), &["gs"]].concat();
                        #[cfg(feature = "azure")]
                        let schemes = [schemes.as_slice(), &["az/azure"]].concat();
                        schemes.join(", ")
                    }
                ),
            )))
        }
    }
}

#[cfg(test)]
mod s3_layer0_tests {
    use super::*;
    use tempfile::TempDir;

    fn local_provider() -> (TempDir, StorageProvider) {
        let dir = TempDir::new().unwrap();
        let loc = StorageLocation::from_path(dir.path());
        let p = StorageProvider::new(loc).unwrap();
        (dir, p)
    }

    #[test]
    fn head_reports_size_for_existing_file() {
        let (_dir, p) = local_provider();
        p.write_bytes("a.bin", vec![1u8; 17]).unwrap();
        let meta = p.head("a.bin").unwrap();
        assert_eq!(meta.size, 17u64);
    }

    #[test]
    fn head_errors_for_missing_file() {
        let (_dir, p) = local_provider();
        assert!(p.head("missing.bin").is_err());
    }

    #[test]
    fn exists_returns_true_for_existing_file() {
        let (_dir, p) = local_provider();
        p.write_bytes("a.bin", vec![1u8; 3]).unwrap();
        assert!(p.exists("a.bin").unwrap());
    }

    #[test]
    fn exists_returns_false_for_missing_file() {
        let (_dir, p) = local_provider();
        assert!(!p.exists("missing.bin").unwrap());
    }

    #[test]
    fn get_bytes_roundtrips() {
        let (_dir, p) = local_provider();
        let payload = b"hello world".to_vec();
        p.write_bytes("a.bin", payload.clone()).unwrap();
        let out = p.get_bytes("a.bin").unwrap();
        assert_eq!(out.as_ref(), payload.as_slice());
    }

    #[test]
    fn range_get_returns_exact_slice() {
        let (_dir, p) = local_provider();
        let payload: Vec<u8> = (0..100u8).collect();
        p.write_bytes("a.bin", payload.clone()).unwrap();
        let slice = p.range_get("a.bin", 10..50).unwrap();
        assert_eq!(slice.as_ref(), &payload[10..50]);
    }

    #[test]
    fn range_get_errors_on_out_of_bounds() {
        let (_dir, p) = local_provider();
        p.write_bytes("a.bin", vec![0u8; 10]).unwrap();
        assert!(p.range_get("a.bin", 0..100).is_err());
    }

    #[test]
    fn list_capped_returns_all_below_cap() {
        let (_dir, p) = local_provider();
        p.write_bytes("sample.d/analysis.tdf", vec![0u8; 3])
            .unwrap();
        p.write_bytes("sample.d/analysis.tdf_bin", vec![0u8; 3])
            .unwrap();
        let list = p.list_capped("sample.d", 10).unwrap();
        assert_eq!(list.len(), 2);
        let names: Vec<String> = list.iter().map(|m| m.location.to_string()).collect();
        assert!(names.iter().any(|n| n.ends_with("/analysis.tdf")));
        assert!(names.iter().any(|n| n.ends_with("/analysis.tdf_bin")));
    }

    #[test]
    fn list_capped_errors_on_cap_exceeded() {
        let (_dir, p) = local_provider();
        for i in 0..5 {
            p.write_bytes(&format!("many/f{i}.bin"), vec![0u8; 1])
                .unwrap();
        }
        let err = p.list_capped("many", 3).unwrap_err();
        assert!(matches!(
            err,
            SerializationError::PrefixCapExceeded { cap: 3, .. }
        ));
    }

    #[test]
    fn get_to_file_copies_object_to_local_path() {
        let (_dir, p) = local_provider();
        p.write_bytes("a.bin", vec![7u8; 100]).unwrap();
        let out_dir = TempDir::new().unwrap();
        let dst = out_dir.path().join("a_copy.bin");
        let bar = indicatif::ProgressBar::hidden();
        p.get_to_file("a.bin", &dst, &bar).unwrap();
        let got = std::fs::read(&dst).unwrap();
        assert_eq!(got, vec![7u8; 100]);
    }

    #[test]
    fn write_bytes_uses_multipart_for_large_buffers() {
        // Local filesystem doesn't care; exercise the path.
        let (_dir, p) = local_provider();
        let big = vec![9u8; 16 * 1024 * 1024];
        p.write_bytes("big.bin", big.clone()).unwrap();
        let got = p.get_bytes("big.bin").unwrap();
        assert_eq!(got.as_ref(), big.as_slice());
    }
}
