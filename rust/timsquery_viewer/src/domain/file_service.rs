//! File loading service

use std::path::Path;
use std::sync::Arc;
use timscentroid::IndexedTimstofPeaks;
use timsquery::models::tolerance::Tolerance;
use timsquery::serde::{
    IndexLoadConfig,
    IndexedPeaksHandle,
    load_index_auto,
};
use tracing::info;

use crate::error::ViewerError;
use crate::file_loader::ElutionGroupData;

/// Service for loading files
///
/// JSPP: NGL ... this file is stupid ... should refactor to just have free functions
/// instead of a struct with no state
pub struct FileService;

impl FileService {
    /// Load elution groups from a JSON file
    ///
    /// # Arguments
    /// * `path` - Path to the elution groups file (.json or .txt)
    ///
    /// # Returns
    /// Parsed elution group data
    pub fn load_elution_groups(path: &Path) -> Result<ElutionGroupData, ViewerError> {
        let res = timsquery::serde::read_library_file(path)?;
        info!(
            "Loaded {} elution groups from {}",
            res.len(),
            path.display()
        );
        Ok(ElutionGroupData::new(res))
    }

    /// Load and index raw timsTOF data from a location (path or URL)
    ///
    /// Supports both local paths and cloud URLs (s3://, gs://, az://).
    /// Automatically detects input type and loads appropriately.
    ///
    /// For cloud cached indexes (.idx), uses lazy loading for faster initialization
    /// (loads metadata only, row groups fetched on-demand during queries).
    ///
    /// Note: This operation may take 10-30 seconds for large datasets when
    /// loading from raw .d files. Cached .idx files load much faster.
    ///
    /// # Arguments
    /// * `location` - Path or URL to the .d directory or .idx cache
    ///
    /// # Returns
    /// Indexed peaks loaded into memory
    pub fn load_raw_data_from_location(
        location: &str,
    ) -> Result<Arc<IndexedPeaksHandle>, ViewerError> {
        // Detect if cloud URL
        let is_cloud = location.contains("://") && !location.starts_with("file://");
        // Use lazy loading for cloud cached indexes (fast init, load on query)
        let prefer_lazy = is_cloud;

        info!(
            "Loading raw data from {}: is_cloud={}, prefer_lazy={}",
            location, is_cloud, prefer_lazy
        );

        let config = IndexLoadConfig {
            prefer_lazy,
            ..Default::default()
        };

        let index = load_index_auto(location, Some(config))
            .map_err(|e| ViewerError::General(format!("Failed to load index: {:?}", e)))?;

        Ok(Arc::new(index))
    }

    /// Load and index raw timsTOF data (legacy method for local paths)
    ///
    /// Supports both local paths and cloud URLs (s3://, gs://, az://).
    /// Automatically detects input type and loads appropriately.
    ///
    /// Note: This operation may take 10-30 seconds for large datasets when
    /// loading from raw .d files. Cached .idx files load much faster.
    ///
    /// # Arguments
    /// * `path` - Path to the .d directory, .idx cache, or cloud URL
    ///
    /// # Returns
    /// Indexed peaks loaded into memory
    pub fn load_raw_data(path: &Path) -> Result<Arc<IndexedPeaksHandle>, ViewerError> {
        let path_str = path
            .to_str()
            .ok_or_else(|| ViewerError::General("Invalid path encoding".to_string()))?;

        Self::load_raw_data_from_location(path_str)
    }

    /// Load tolerance settings from a JSON file
    ///
    /// # Arguments
    /// * `path` - Path to the tolerance JSON file
    ///
    /// # Returns
    /// Parsed tolerance settings
    pub fn load_tolerance(path: &Path) -> Result<Tolerance, ViewerError> {
        let file_content = std::fs::read_to_string(path)?;
        let tolerance: Tolerance = serde_json::from_str(&file_content)?;
        Ok(tolerance)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_file_service_structure() {
        let _service = FileService;
    }
}
