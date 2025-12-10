use crate::{
    CentroidingConfig,
    IndexedTimstofPeaks,
    TimsTofPath,
};
use std::path::Path;
use timscentroid::serialization::SerializationConfig;
use tracing::{
    error,
    info,
};

fn maybe_cache_load_index(index_cache_loc: impl AsRef<Path>) -> Option<IndexedTimstofPeaks> {
    info!(
        "Attempting to load index from cache at {:?}",
        index_cache_loc.as_ref()
    );
    match IndexedTimstofPeaks::load_from_directory(index_cache_loc.as_ref()) {
        Ok(idx) => {
            info!("Loaded index from cache at {:?}", index_cache_loc.as_ref());
            Some(idx)
        }
        Err(e) => {
            error!(
                "Failed to load index from cache at {:?}: {:?}",
                index_cache_loc.as_ref(),
                e
            );
            None
        }
    }
}

/// Builder for loading timsTOF indices with caching options
pub struct TimsIndexReader {
    write_missing_cache: bool,
    centroiding_config: Option<CentroidingConfig>,
    serialization_config: SerializationConfig,
}

impl TimsIndexReader {
    /// Create a new index reader with default settings
    pub fn new() -> Self {
        Self {
            write_missing_cache: true,
            centroiding_config: None,
            serialization_config: SerializationConfig::default(),
        }
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

        // Create cache directory path by appending .idx to the .d directory
        let mut index_location = file_location.as_ref().to_path_buf();
        let new_name = format!(
            "{}.idx",
            index_location.file_name().unwrap().to_str().unwrap()
        );
        index_location.set_file_name(new_name);

        let out = if let Some(idx) = maybe_cache_load_index(&index_location) {
            Ok(idx)
        } else {
            let cache_loc = if self.write_missing_cache {
                Some(index_location)
            } else {
                None
            };
            Ok(self.uncached_load_index(&timstofpath, &cache_loc))
        };

        let et = st.elapsed();
        info!("Loading index took: {:#?}", et);
        out
    }

    fn uncached_load_index(
        &self,
        timstofpath: &TimsTofPath,
        cache_loc: &Option<std::path::PathBuf>,
    ) -> IndexedTimstofPeaks {
        let centroiding_config = self.centroiding_config.unwrap_or(CentroidingConfig {
            max_peaks: 50_000,
            mz_ppm_tol: 10.0,
            im_pct_tol: 5.0,
            early_stop_iterations: 200,
        });
        info!("Using centroiding config: {:#?}", centroiding_config);
        info!("Starting centroiging + load of the raw data (might take a min)");
        let (index, build_stats) =
            IndexedTimstofPeaks::from_timstof_file(timstofpath, centroiding_config);
        info!("Index built with stats: {}", build_stats);

        // Save to cache
        if let Some(idx_path) = cache_loc {
            info!("Saving index to cache at {:?}", idx_path);
            if let Err(e) = index.save_to_directory_with_config(idx_path, self.serialization_config)
            {
                error!("Failed to save index to cache: {:?}", e);
            } else {
                info!("Saved index to cache");
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

/// Convenience function for loading index with default caching behavior
pub fn load_index_caching(
    file_location: impl AsRef<Path>,
) -> Result<IndexedTimstofPeaks, crate::errors::DataReadingError> {
    TimsIndexReader::new().read_index(file_location)
}
