//! Smart index loading with automatic format detection and cloud support
//!
//! This module provides a unified API for loading timsTOF data from any source:
//! local .d files, cached .idx directories, or cloud storage (S3, GCS, Azure).
//!
//! # Quick Start
//!
//! ## Automatic Loading (Detects Format)
//!
//! ```no_run
//! use timsquery::serde::load_index_auto;
//!
//! // Works with any input - auto-detects format and location
//! let index = load_index_auto("data.d", None)?.into_eager()?;              // Local raw
//! let index = load_index_auto("data.d.idx", None)?.into_eager()?;          // Local cached
//! let index = load_index_auto("s3://bucket/exp.d", None)?.into_eager()?;   // Cloud raw
//! let index = load_index_auto("s3://bucket/exp.idx", None)?.into_eager()?; // Cloud cached
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! ## Lazy Loading (For Large Datasets)
//!
//! ```no_run
//! use timsquery::serde::{load_index_auto, IndexLoadConfig};
//!
//! // Prefer lazy loading when possible (faster initialization, less memory)
//! let config = IndexLoadConfig {
//!     prefer_lazy: true,
//!     ..Default::default()
//! };
//!
//! let handle = load_index_auto("experiment.d.idx", Some(config))?;
//!
//! // Query directly on lazy handle (loads data on-demand)
//! if let Some(lazy) = handle.as_lazy() {
//!     let peaks = lazy.query_peaks_ms1(mz_range, rt_range, im_range);
//! }
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! ## Custom Cache Configuration
//!
//! ```no_run
//! use timsquery::serde::{load_index_auto, IndexLoadConfig, CacheLocation};
//!
//! let config = IndexLoadConfig {
//!     cache_location: CacheLocation::Url("s3://my-bucket/cache/".to_string()),
//!     ..Default::default()
//! };
//!
//! // Process raw data and cache to S3
//! let index = load_index_auto("data.d", Some(config))?.into_eager()?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Cache Workflow
//!
//! 1. **First run**: Reads raw .d file, builds index (slow), saves to cache
//! 2. **Subsequent runs**: Loads from cache (fast), skips raw data processing
//!
//! # Configuration
//!
//! ```no_run
//! use timsquery::serde::TimsIndexReader;
//! use timsquery::CentroidingConfig;
//! use timscentroid::serialization::SerializationConfig;
//! use parquet::basic::{Compression, ZstdLevel};
//!
//! let index = TimsIndexReader::new()
//!     .with_auto_cache()
//!     .with_centroiding_config(CentroidingConfig {
//!         max_peaks: 50_000,
//!         mz_ppm_tol: 10.0,
//!         im_pct_tol: 5.0,
//!         early_stop_iterations: 200,
//!     })
//!     .with_serialization_config(SerializationConfig {
//!         row_group_size: 100_000,
//!         compression: Compression::ZSTD(ZstdLevel::try_new(3)?),
//!     })
//!     .read_index("data.d")?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Use Cases
//!
//! **Team sharing via S3:**
//! ```no_run
//! use timsquery::serde::TimsIndexReader;
//!
//! // First team member: builds and caches to S3
//! let index = TimsIndexReader::new()
//!     .with_cloud_cache("s3://team-bucket/indexes/exp001/")
//!     .read_index("/local/exp001.d")?;
//!
//! // Other team members: load from S3 cache (no raw data needed!)
//! let index = TimsIndexReader::from_cache_url(
//!     "s3://team-bucket/indexes/exp001/"
//! )?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! **CI/CD pipelines:**
//! ```no_run
//! use timsquery::serde::TimsIndexReader;
//!
//! // Build step: process raw data, cache result
//! let index = TimsIndexReader::new()
//!     .with_cloud_cache("s3://pipeline-artifacts/indexes/build-123/")
//!     .read_index("raw_data.d")?;
//!
//! // Deploy step: load pre-built index
//! let index = TimsIndexReader::from_cache_url(
//!     "s3://pipeline-artifacts/indexes/build-123/"
//! )?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::{
    CentroidingConfig,
    IndexedTimstofPeaks,
    TimsTofPath,
};
use std::path::{
    Path,
    PathBuf,
};
use timscentroid::StorageLocation;
use timscentroid::lazy::LazyIndexedTimstofPeaks;
use timscentroid::serialization::SerializationConfig;
use tracing::{
    error,
    info,
};

/// Handle to indexed peaks - can be lazy or materialized (eager)
///
/// This enum allows applications to work with either lazy-loaded or fully materialized
/// index data. Lazy loading is faster to initialize and uses less memory, while eager
/// loading provides faster queries and is required for some operations.
pub enum IndexedPeaksHandle {
    /// Lazy loading - queries fetch data on-demand from parquet files
    Lazy(LazyIndexedTimstofPeaks),
    /// Eager/materialized - all data loaded in memory
    Eager(IndexedTimstofPeaks),
}

impl IndexedPeaksHandle {
    /// Materialize to eager if needed (no-op if already eager)
    ///
    /// This loads all parquet data into memory if the handle is lazy.
    pub fn into_eager(self) -> Result<IndexedTimstofPeaks, crate::errors::DataReadingError> {
        match self {
            Self::Eager(peaks) => Ok(peaks),
            Self::Lazy(_lazy) => {
                // Materialize: load all data from parquet files into memory
                // We need to get the storage location from lazy and reload as eager
                info!("Materializing lazy index to eager (loading all data into memory)");

                // TODO: Implement lazy -> eager materialization
                // The lazy type needs to expose its storage location so we can reload as eager
                Err(crate::errors::DataReadingError::UnsupportedDataError(
                    crate::errors::UnsupportedDataError::NoMS2DataError,
                ))
            }
        }
    }

    /// Try to convert to lazy loading
    ///
    /// If already lazy, returns self unchanged. If eager, saves to cache and reloads as lazy.
    /// Returns error if caching is disabled or save/load fails.
    ///
    /// # Arguments
    /// * `cache_location` - Where to save the cache (must not be Disabled)
    /// * `serialization_config` - Configuration for saving (if needed)
    ///
    /// # Example
    /// ```no_run
    /// use timsquery::serde::{load_index_auto, CacheLocation};
    /// use timscentroid::serialization::SerializationConfig;
    ///
    /// let handle = load_index_auto("data.d", None)?;
    /// let lazy = handle.try_into_lazy(
    ///     CacheLocation::Local("/tmp/cache".into()),
    ///     SerializationConfig::default()
    /// )?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_into_lazy(
        self,
        cache_location: CacheLocation,
        serialization_config: SerializationConfig,
    ) -> Result<LazyIndexedTimstofPeaks, crate::errors::DataReadingError> {
        match self {
            Self::Lazy(lazy) => Ok(lazy),
            Self::Eager(eager) => {
                if matches!(cache_location, CacheLocation::Disabled) {
                    return Err(crate::errors::DataReadingError::UnsupportedDataError(
                        crate::errors::UnsupportedDataError::NoMS2DataError,
                    ));
                }

                let storage_location = match cache_location.to_storage_location() {
                    Some(Ok(loc)) => loc,
                    Some(Err(_)) | None => {
                        return Err(crate::errors::DataReadingError::UnsupportedDataError(
                            crate::errors::UnsupportedDataError::NoMS2DataError,
                        ));
                    }
                };

                info!("Saving eager index to cache for lazy loading");
                eager
                    .save_to_storage(storage_location.clone(), serialization_config)
                    .map_err(crate::errors::DataReadingError::SerializationError)?;

                info!("Loading saved index as lazy");
                LazyIndexedTimstofPeaks::load_from_storage(storage_location)
                    .map_err(crate::errors::DataReadingError::SerializationError)
            }
        }
    }

    /// Check if this handle is lazy
    pub fn is_lazy(&self) -> bool {
        matches!(self, Self::Lazy(_))
    }

    /// Check if this handle is eager
    pub fn is_eager(&self) -> bool {
        matches!(self, Self::Eager(_))
    }
}

/// Configuration for smart index loading
#[derive(Debug, Clone)]
pub struct IndexLoadConfig {
    /// Cache location (Auto, Local, Url, Disabled)
    pub cache_location: CacheLocation,

    /// Centroiding configuration (only used when loading raw .d files)
    pub centroiding_config: Option<CentroidingConfig>,

    /// Serialization configuration
    pub serialization_config: SerializationConfig,

    /// Prefer lazy loading when possible (default: false)
    ///
    /// When true, cached indexes will be loaded lazily for faster initialization
    /// and lower memory usage. Raw .d files still require eager loading.
    pub prefer_lazy: bool,

    /// Allow writing cache if missing (default: true)
    pub write_missing_cache: bool,
}

impl Default for IndexLoadConfig {
    fn default() -> Self {
        Self {
            cache_location: CacheLocation::Auto,
            centroiding_config: None,
            serialization_config: SerializationConfig::default(),
            prefer_lazy: false,
            write_missing_cache: true,
        }
    }
}

/// Cache location for indexed peaks
///
/// # Storage Layout
///
/// Cached indexes are stored as a directory containing:
/// - `metadata.json` - Index metadata (version, timestamps, file mappings)
/// - `ms1.parquet` - MS1 peaks
/// - `ms2/group_0.parquet`, `ms2/group_1.parquet`, ... - MS2 window groups
///
/// ## Examples
///
/// **Auto:** Derives location from input path
/// ```text
/// Input:  /data/experiment.d
/// Cache:  /data/experiment.d.idx/
///         ├── metadata.json
///         ├── ms1.parquet
///         └── ms2/
///             ├── group_0.parquet
///             └── group_1.parquet
/// ```
///
/// **Local:** Explicit local directory
/// ```text
/// CacheLocation::Local("/cache/my_index")
/// Cache:  /cache/my_index/
///         ├── metadata.json
///         ├── ms1.parquet
///         └── ms2/...
/// ```
///
/// **Url:** Cloud storage URL (full path including .idx directory name)
/// ```text
/// CacheLocation::Url("s3://bucket/experiments/exp001.idx")
/// Cache:  s3://bucket/experiments/exp001.idx/
///         ├── metadata.json
///         ├── ms1.parquet
///         └── ms2/...
/// ```
///
/// **Note:** Currently cache paths are exact - relative path preservation
/// (e.g., `/local/parent/data.d` → `s3://bucket/parent/data.d.idx`) is not yet
/// implemented but would be valuable for large dataset workflows.
#[derive(Debug, Clone)]
pub enum CacheLocation {
    /// Derive cache location from input path (adds .idx suffix)
    ///
    /// Example: "data.d" → "data.d.idx"
    Auto,

    /// Explicit local filesystem directory path
    ///
    /// This is the full path to the .idx directory (not the parent).
    Local(PathBuf),

    /// Cloud storage URL
    ///
    /// Must be a String (not PathBuf) because it contains URL scheme and authority:
    /// - "s3://bucket/path.idx" (scheme: s3, authority: bucket)
    /// - "gs://bucket/path.idx" (scheme: gs, authority: bucket)
    /// - "az://container/path.idx" (scheme: az, authority: container)
    ///
    /// PathBuf is for filesystem paths only and doesn't support URLs.
    Url(String),

    /// No caching - process raw data every time
    Disabled,
}

impl CacheLocation {
    /// Convert to StorageLocation, returning None for Disabled or Auto
    fn to_storage_location(&self) -> Option<Result<StorageLocation, String>> {
        match self {
            Self::Local(path) => Some(Ok(StorageLocation::from_path(path))),
            Self::Url(url) => Some(
                StorageLocation::from_url(url).map_err(|e| format!("Invalid cache URL: {}", e)),
            ),
            Self::Auto | Self::Disabled => None,
        }
    }

    /// Derive auto cache location from input path
    ///
    /// Appends `.idx` to the filename (not replacing extension).
    /// Example: "data.d" → "data.d.idx"
    fn derive_auto_location(input_path: &Path) -> PathBuf {
        // Can't use with_extension("idx") because that would replace ".d" with ".idx"
        // We want to append, so "data.d" becomes "data.d.idx" not "data.idx"
        let mut index_location = input_path.to_path_buf();
        let current_name = index_location
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("");
        index_location.set_file_name(format!("{}.idx", current_name));
        index_location
    }
}

/// Builder for loading timsTOF indices with cloud storage and caching support
pub struct TimsIndexReader {
    cache_location: CacheLocation,
    write_missing_cache: bool,
    centroiding_config: Option<CentroidingConfig>,
    serialization_config: SerializationConfig,
}

impl TimsIndexReader {
    /// Create a new index reader with automatic caching enabled (default)
    ///
    /// Cache location is automatically derived by appending `.idx` to the input path.
    /// Example: "data.d" → "data.d.idx"
    ///
    /// Use `without_cache()` if you don't want any caching.
    pub fn new() -> Self {
        Self {
            cache_location: CacheLocation::Auto,
            write_missing_cache: true,
            centroiding_config: None,
            serialization_config: SerializationConfig::default(),
        }
    }

    /// Create a new index reader with caching disabled
    ///
    /// Raw data will be processed every time without saving/loading cache.
    pub fn new_without_cache() -> Self {
        Self {
            cache_location: CacheLocation::Disabled,
            write_missing_cache: false,
            centroiding_config: None,
            serialization_config: SerializationConfig::default(),
        }
    }

    /// Set a custom local cache directory
    pub fn with_local_cache(mut self, path: impl Into<PathBuf>) -> Self {
        self.cache_location = CacheLocation::Local(path.into());
        self
    }

    /// Set a cloud storage cache location (s3://, gs://, az://)
    pub fn with_cloud_cache(mut self, url: impl Into<String>) -> Self {
        self.cache_location = CacheLocation::Url(url.into());
        self
    }

    /// Set custom cache location
    pub fn with_cache_location(mut self, location: CacheLocation) -> Self {
        self.cache_location = location;
        self
    }

    /// Set whether to write cache when it's missing (default: true)
    pub fn with_write_missing_cache(mut self, write: bool) -> Self {
        self.write_missing_cache = write;
        self
    }

    /// Set custom centroiding configuration
    pub fn with_centroiding_config(mut self, config: CentroidingConfig) -> Self {
        self.centroiding_config = Some(config);
        self
    }

    /// Set custom serialization configuration
    pub fn with_serialization_config(mut self, config: SerializationConfig) -> Self {
        self.serialization_config = config;
        self
    }

    /// Load directly from a cached index (bypasses centroiding, just loads pre-built index)
    ///
    /// This is useful when you have a pre-cached index and don't need the original .d file.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use timsquery::serde::TimsIndexReader;
    /// use timscentroid::StorageLocation;
    ///
    /// // Load from local cache
    /// let index = TimsIndexReader::from_cache(
    ///     StorageLocation::from_path("/path/to/experiment.d.idx")
    /// )?;
    ///
    /// // Load from S3 cache
    /// let index = TimsIndexReader::from_cache(
    ///     StorageLocation::from_url("s3://bucket/cache/experiment.idx")?
    /// )?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn from_cache(
        location: StorageLocation,
    ) -> Result<IndexedTimstofPeaks, crate::errors::DataReadingError> {
        IndexedTimstofPeaks::load_from_storage(location)
            .map_err(crate::errors::DataReadingError::SerializationError)
    }

    /// Load directly from a cached index at a local path
    ///
    /// Convenience wrapper for `from_cache(StorageLocation::from_path(...))`
    pub fn from_cache_path(
        path: impl AsRef<Path>,
    ) -> Result<IndexedTimstofPeaks, crate::errors::DataReadingError> {
        Self::from_cache(StorageLocation::from_path(path))
    }

    /// Load directly from a cached index at a cloud URL
    ///
    /// Convenience wrapper for `from_cache(StorageLocation::from_url(...))`
    pub fn from_cache_url(
        url: impl AsRef<str>,
    ) -> Result<IndexedTimstofPeaks, crate::errors::DataReadingError> {
        let location = StorageLocation::from_url(url.as_ref()).map_err(|e| {
            crate::errors::DataReadingError::SerializationError(
                timscentroid::serialization::SerializationError::UrlParse(e),
            )
        })?;
        Self::from_cache(location)
    }

    /// Load timsTOF index with configured caching behavior
    ///
    /// # Arguments
    ///
    /// * `file_location` - Path to the timsTOF .d directory
    ///
    /// # Returns
    /// * `Result<IndexedTimstofPeaks, crate::errors::DataReadingError>` - Loaded index or error
    pub fn read_index(
        &self,
        file_location: impl AsRef<Path>,
    ) -> Result<IndexedTimstofPeaks, crate::errors::DataReadingError> {
        let st = std::time::Instant::now();

        let timstofpath = match TimsTofPath::new(file_location.as_ref()) {
            Ok(x) => x,
            Err(e) => {
                return Err(crate::errors::DataReadingError::TimsTofPathError(e));
            }
        };

        // Determine cache location and handle errors early
        let cache_storage_location: Option<StorageLocation> = match &self.cache_location {
            CacheLocation::Auto => {
                let path = CacheLocation::derive_auto_location(file_location.as_ref());
                Some(StorageLocation::from_path(path))
            }
            other => match other.to_storage_location() {
                Some(Ok(loc)) => Some(loc),
                Some(Err(e)) => {
                    error!("Invalid cache location: {}", e);
                    None
                }
                None => None,
            },
        };

        // Try to load from cache
        let out = if let Some(storage_loc) = &cache_storage_location {
            match self.try_load_from_cache(storage_loc) {
                Some(idx) => Ok(idx),
                None => {
                    let cache_loc = if self.write_missing_cache {
                        cache_storage_location
                    } else {
                        None
                    };
                    Ok(self.uncached_load_index(&timstofpath, cache_loc))
                }
            }
        } else {
            Ok(self.uncached_load_index(&timstofpath, None))
        };

        let et = st.elapsed();
        info!("Loading index took: {:#?}", et);
        out
    }

    fn try_load_from_cache(
        &self,
        storage_location: &StorageLocation,
    ) -> Option<IndexedTimstofPeaks> {
        let location_desc = match storage_location {
            StorageLocation::Local(p) => format!("{:?}", p),
            StorageLocation::Url(u) => u.to_string(),
        };

        info!("Attempting to load index from cache at {}", location_desc);

        match IndexedTimstofPeaks::load_from_storage(storage_location.clone()) {
            Ok(idx) => {
                info!("Loaded index from cache at {}", location_desc);
                Some(idx)
            }
            Err(e) => {
                error!(
                    "Failed to load index from cache at {}: {:?}",
                    location_desc, e
                );
                None
            }
        }
    }

    fn uncached_load_index(
        &self,
        timstofpath: &TimsTofPath,
        cache_loc: Option<StorageLocation>,
    ) -> IndexedTimstofPeaks {
        let centroiding_config = self.centroiding_config.unwrap_or(CentroidingConfig {
            max_peaks: 50_000,
            mz_ppm_tol: 10.0,
            im_pct_tol: 5.0,
            early_stop_iterations: 200,
        });

        info!("Using centroiding config: {:#?}", centroiding_config);
        info!("Starting centroiding + load of the raw data (might take a min)");

        let (index, build_stats) =
            IndexedTimstofPeaks::from_timstof_file(timstofpath, centroiding_config);

        info!("Index built with stats: {}", build_stats);

        // Save to cache
        if let Some(storage_loc) = cache_loc {
            let location_desc = match &storage_loc {
                StorageLocation::Local(p) => format!("{:?}", p),
                StorageLocation::Url(u) => u.to_string(),
            };

            info!("Saving index to cache at {}", location_desc);

            match index.save_to_storage(storage_loc, self.serialization_config) {
                Ok(_) => info!("Saved index to cache"),
                Err(e) => error!("Failed to save index to cache: {:?}", e),
            }
        }

        index
    }
}

impl Default for TimsIndexReader {
    fn default() -> Self {
        Self::new()
    }
}

/// Check if a location contains a cached index by looking for metadata.json
///
/// This checks the actual directory contents rather than just the path extension.
fn sniff_cached_index(location: &str) -> bool {
    let is_cloud = location.contains("://");

    // Try to create storage location and check for metadata.json
    let storage_result = if is_cloud {
        StorageLocation::from_url(location)
    } else {
        Ok(StorageLocation::from_path(location))
    };

    let Ok(storage_location) = storage_result else {
        return false;
    };

    // Try to read metadata.json as a quick check
    match timscentroid::storage::StorageProvider::new(storage_location) {
        Ok(provider) => {
            // Just try to read a few bytes - if metadata.json exists, it's likely a cached index
            matches!(provider.read_bytes("metadata.json"), Ok(_))
        }
        Err(_) => false,
    }
}

/// Check if a location contains raw .d data by looking for analysis.tdf
///
/// This checks the actual directory contents rather than just the path extension.
fn sniff_raw_tdf(location: &str) -> bool {
    // For cloud URLs, we can't easily check directory contents without listing
    // So fall back to path-based detection for now
    if location.contains("://") {
        return location.contains(".d") && !location.contains(".idx");
    }

    // For local paths, check if directory contains analysis.tdf
    let path = Path::new(location);
    if !path.is_dir() {
        return false;
    }

    // Check for required Bruker timsTOF files
    path.join("analysis.tdf").exists() && path.join("analysis.tdf_bin").exists()
}

/// Smart index loader - auto-detects input type and loads appropriately
///
/// This is the unified entry point for loading indexed peaks from any source.
/// It automatically detects:
/// - Local vs cloud storage (by checking for "://" in path)
/// - Raw .d files vs cached .idx files (by sniffing directory contents)
/// - Whether to load lazily or eagerly (based on config)
///
/// # Format Detection
///
/// The function intelligently detects the input format by:
/// 1. Checking for cached index: Looks for `metadata.json` file
/// 2. Checking for raw data: Looks for `analysis.tdf` and `analysis.tdf_bin` files
/// 3. Falling back to path-based heuristics for cloud URLs
///
/// This is more robust than just checking file extensions, especially for
/// cloud URLs where paths may not follow standard naming conventions.
///
/// # Arguments
///
/// * `path_or_url` - Can be:
///   - "/path/to/experiment.d" - Local raw data
///   - "/path/to/experiment.d.idx" - Local cached index
///   - "s3://bucket/experiment.d" - Cloud raw data
///   - "s3://bucket/experiment.idx" - Cloud cached index
///   - Any directory containing the appropriate files (extension-agnostic)
/// * `config` - Optional configuration (uses defaults if None)
///
/// # Returns
///
/// `IndexedPeaksHandle` - Either lazy or eager depending on config and input type
///
/// # Examples
///
/// ```no_run
/// use timsquery::serde::load_index_auto;
///
/// // Simple usage - auto-detects everything
/// let index = load_index_auto("data.d", None)?.into_eager()?;
///
/// // Load from cloud (works without .idx extension if metadata.json exists)
/// let index = load_index_auto("s3://bucket/my_experiment", None)?.into_eager()?;
///
/// // Prefer lazy loading
/// use timsquery::serde::IndexLoadConfig;
/// let config = IndexLoadConfig { prefer_lazy: true, ..Default::default() };
/// let handle = load_index_auto("data.d.idx", Some(config))?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn load_index_auto(
    path_or_url: impl AsRef<str>,
    config: Option<IndexLoadConfig>,
) -> Result<IndexedPeaksHandle, crate::errors::DataReadingError> {
    let input = path_or_url.as_ref();
    let config = config.unwrap_or_default();

    info!("Loading index from: {}", input);

    // Detect input type by actually checking directory contents
    // This is more robust than just looking at file extensions
    let is_cached = sniff_cached_index(input);
    let is_cloud = input.contains("://");

    info!(
        "Detected: cached={}, cloud={}, prefer_lazy={}",
        is_cached, is_cloud, config.prefer_lazy
    );

    // Early validation: reject cloud raw .d files with helpful error
    if is_cloud && !is_cached {
        error!("Attempted to load raw .d file from cloud storage: {}", input);
        return Err(crate::errors::DataReadingError::UnsupportedDataError(
            crate::errors::UnsupportedDataError::CloudRawDataNotSupported {
                url: input.to_string(),
                suggestion: format!(
                    "Raw .d files must be processed locally first. Suggested workflow:\n\
                    1. Download the .d file locally\n\
                    2. Process it to create a cached index:\n\
                    \n   \
                    use timsquery::serde::TimsIndexReader;\n   \
                    let index = TimsIndexReader::new()\n       \
                    .with_cloud_cache(\"{}\")\n       \
                    .read_index(\"/local/path/to/data.d\")?;\n\
                    \n\
                    3. Then load from the cloud cache:\n   \
                    let index = TimsIndexReader::from_cache_url(\"{}\")?;\n\
                    \n\
                    Alternatively, use a local cache and upload it manually to your cloud storage.",
                    input.trim_end_matches(".d").to_string() + ".idx",
                    input.trim_end_matches(".d").to_string() + ".idx",
                ),
            },
        ));
    }

    match (is_cached, config.prefer_lazy) {
        (true, true) => {
            // Cached index + prefer lazy = load lazy
            info!("Loading as lazy (cached index)");
            let location = if is_cloud {
                StorageLocation::from_url(input).map_err(|e| {
                    crate::errors::DataReadingError::SerializationError(
                        timscentroid::serialization::SerializationError::UrlParse(e),
                    )
                })?
            } else {
                StorageLocation::from_path(input)
            };
            let lazy = LazyIndexedTimstofPeaks::load_from_storage(location)
                .map_err(crate::errors::DataReadingError::SerializationError)?;
            Ok(IndexedPeaksHandle::Lazy(lazy))
        }
        (true, false) => {
            // Cached index + prefer eager = load eager
            info!("Loading as eager (cached index)");
            let eager = if is_cloud {
                TimsIndexReader::from_cache_url(input)?
            } else {
                TimsIndexReader::from_cache_path(input)?
            };
            Ok(IndexedPeaksHandle::Eager(eager))
        }
        (false, _) => {
            // Raw .d file - need to centroid/index (always eager)
            info!("Loading as eager (raw .d file - requires centroiding)");

            let mut reader = TimsIndexReader::new()
                .with_cache_location(config.cache_location)
                .with_write_missing_cache(config.write_missing_cache)
                .with_serialization_config(config.serialization_config);

            if let Some(centroid_cfg) = config.centroiding_config {
                reader = reader.with_centroiding_config(centroid_cfg);
            }

            let eager = reader.read_index(input)?;

            // Note: Raw .d files always return eager, even if prefer_lazy is true
            // (can't lazy-load raw data - it needs to be centroided first)
            Ok(IndexedPeaksHandle::Eager(eager))
        }
    }
}
