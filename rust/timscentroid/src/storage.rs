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
    info,
    instrument,
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
        info!("Reading from full path: {}", full_path);
        let object_path = ObjectPath::from(full_path.as_str());
        info!("Object path: {:?}", object_path);
        let result = match self.store.get(&object_path).await {
            Ok(res) => res,
            Err(e) => {
                // Categorize the error properly based on what it actually is
                let error_str = e.to_string();
                info!("Error getting object: {:?}", e);

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
        let full_path = self.build_path(path);
        let object_path = ObjectPath::from(full_path.as_str());
        block_on_or_in_place(async {
            self.store
                .put(&object_path, Bytes::from(data).into())
                .await?;
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

    /// Read indexed peaks from a parquet file
    ///
    /// Uses `ParquetObjectReader` with async streaming for both local and cloud storage.
    /// This provides efficient reading without intermediate copies or temp files.
    ///
    /// # Arguments
    /// * `path` - Relative path to the parquet file
    ///
    /// # Returns
    /// Vector of indexed peaks loaded from the parquet file
    #[instrument(skip(self), fields(path = %path))]
    pub fn read_parquet_peaks<T: crate::rt_mapping::RTIndex>(
        &self,
        path: &str,
    ) -> Result<Vec<crate::indexing::IndexedPeak<T>>, SerializationError> {
        use futures::stream::StreamExt;
        use parquet::arrow::ParquetRecordBatchStreamBuilder;
        use parquet::arrow::async_reader::ParquetObjectReader;

        let full_path = self.build_path(path);
        let object_path = ObjectPath::from(full_path.as_str());

        block_on_or_in_place(async {
            // Create ParquetObjectReader - works for both local and cloud storage
            let reader = ParquetObjectReader::new(self.store.clone(), object_path);

            // Build stream
            let builder = ParquetRecordBatchStreamBuilder::new(reader).await?;
            let mut stream = builder.build()?;

            let mut peaks = Vec::new();

            // Stream record batches
            while let Some(batch_result) = stream.next().await {
                let batch = batch_result?;
                // Convert batch to peaks (reuse logic from serialization module)
                peaks.extend(crate::serialization::batch_to_peaks::<T>(&batch)?);
            }

            Ok(peaks)
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
            info!("Creating S3 ObjectStore for bucket: {}", bucket);

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
