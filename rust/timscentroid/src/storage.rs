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

use crate::serialization::SerializationError;

/// Global tokio runtime for all async operations (created lazily)
pub(crate) static RUNTIME: Lazy<Runtime> = Lazy::new(|| {
    tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()
        .expect("Failed to create Tokio runtime")
});

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
        if parsed.scheme() == "file" {
            if let Ok(path) = parsed.to_file_path() {
                return Ok(Self::Local(path));
            }
        }

        Ok(Self::Url(parsed))
    }
}

/// Storage provider wrapping an ObjectStore
#[derive(Clone)]
pub struct StorageProvider {
    store: Arc<dyn ObjectStore>,
    is_local: bool,
}

impl StorageProvider {
    /// Create a new storage provider from a location
    pub fn new(location: StorageLocation) -> Result<Self, SerializationError> {
        let (store, is_local): (Arc<dyn ObjectStore>, bool) = match location {
            StorageLocation::Local(path) => {
                // Create directory if it doesn't exist
                std::fs::create_dir_all(&path)?;
                (Arc::new(LocalFileSystem::new_with_prefix(path)?), true)
            }
            StorageLocation::Url(url) => (parse_url(&url)?, false),
        };

        Ok(Self { store, is_local })
    }

    /// Check if this storage provider is backed by local filesystem
    ///
    /// Returns true for local paths, false for cloud storage (S3, GCS, Azure)
    pub fn is_local(&self) -> bool {
        self.is_local
    }

    /// Read a file as bytes (async version)
    pub async fn read_bytes_async(&self, path: &str) -> Result<Vec<u8>, SerializationError> {
        let object_path = ObjectPath::from(path);
        let result = self.store.get(&object_path).await?;
        let bytes = result.bytes().await?;
        Ok(bytes.to_vec())
    }

    /// Read a file as bytes (blocking version)
    pub fn read_bytes(&self, path: &str) -> Result<Vec<u8>, SerializationError> {
        RUNTIME.block_on(self.read_bytes_async(path))
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
        RUNTIME.block_on(self.read_to_string_async(path))
    }

    /// Write bytes to a file
    pub fn write_bytes(&self, path: &str, data: Vec<u8>) -> Result<(), SerializationError> {
        let object_path = ObjectPath::from(path);
        RUNTIME.block_on(async {
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
    pub fn read_parquet_peaks<T: crate::rt_mapping::RTIndex>(
        &self,
        path: &str,
    ) -> Result<Vec<crate::indexing::IndexedPeak<T>>, SerializationError> {
        use futures::stream::StreamExt;
        use parquet::arrow::ParquetRecordBatchStreamBuilder;
        use parquet::arrow::async_reader::ParquetObjectReader;

        let object_path = ObjectPath::from(path);

        RUNTIME.block_on(async {
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
fn parse_url(url: &url::Url) -> Result<Arc<dyn ObjectStore>, SerializationError> {
    match url.scheme() {
        #[cfg(feature = "aws")]
        "s3" => {
            use object_store::aws::AmazonS3Builder;

            let bucket = url.host_str().ok_or_else(|| {
                SerializationError::Io(std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "Missing bucket in S3 URL",
                ))
            })?;

            Ok(Arc::new(
                AmazonS3Builder::from_env()
                    .with_bucket_name(bucket)
                    .build()?,
            ))
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
                        let schemes = vec!["file"];
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
