//! File loading service

use std::path::Path;
use std::sync::Arc;
use timscentroid::IndexedTimstofPeaks;
use timsquery::models::tolerance::Tolerance;
use timsquery::serde::load_index_auto;
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

    /// Load and index raw timsTOF data
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
    pub fn load_raw_data(path: &Path) -> Result<Arc<IndexedTimstofPeaks>, ViewerError> {
        let path_str = path.to_str().ok_or_else(|| ViewerError::General(
            "Invalid path encoding".to_string()
        ))?;

        let index = load_index_auto(path_str, None)
            .map_err(|e| ViewerError::DataLoading {
                path: path.to_path_buf(),
                source: Box::new(ViewerError::General(format!("{:?}", e))),
            })?
            .into_eager()
            .map_err(|e| ViewerError::DataLoading {
                path: path.to_path_buf(),
                source: Box::new(ViewerError::General(format!("{:?}", e))),
            })?;

        Ok(Arc::new(index))
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
