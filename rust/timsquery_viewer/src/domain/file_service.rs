//! File loading service

use std::path::Path;
use std::sync::Arc;
use timscentroid::{
    IndexedTimstofPeaks,
    TimsTofPath,
};
use timsquery::models::tolerance::Tolerance;
use timsquery::serde::load_index_caching;
use timsrust::MSLevel;
use tracing::info;

use crate::error::ViewerError;
use crate::file_loader::ElutionGroupData;

/// Service for loading files
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
        Ok(ElutionGroupData { inner: res })
    }

    /// Load and index raw timsTOF data
    ///
    /// Note: This operation may take 10-30 seconds for large datasets.
    ///
    /// # Arguments
    /// * `path` - Path to the .d directory
    ///
    /// # Returns
    /// A tuple of (indexed peaks, MS1 retention times in milliseconds)
    pub fn load_raw_data(
        path: &Path,
    ) -> Result<(Arc<IndexedTimstofPeaks>, Arc<[u32]>), ViewerError> {
        let index = load_index_caching(path).map_err(|e| ViewerError::DataLoading {
            path: path.to_path_buf(),
            source: Box::new(ViewerError::General(format!("{:?}", e))),
        })?;

        let rts = Self::get_ms1_rts_as_millis(path)?;

        Ok((Arc::new(index), rts))
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

    /// Retrieves MS1 retention times from a TIMS-TOF file, sorted and deduped
    ///
    /// # Arguments
    /// * `path` - Path to the .d directory
    ///
    /// # Returns
    /// MS1 retention times in milliseconds
    fn get_ms1_rts_as_millis(path: &Path) -> Result<Arc<[u32]>, ViewerError> {
        let ttp = TimsTofPath::new(path).map_err(|e| ViewerError::TimsFileLoad {
            path: path.to_path_buf(),
            source: e,
        })?;
        let reader = ttp.load_frame_reader()?;
        let mut rts: Vec<_> = reader
            .frame_metas
            .iter()
            .filter(|x| x.ms_level == MSLevel::MS1)
            .map(|f| (f.rt_in_seconds * 1000.0).round() as u32)
            .collect();
        rts.sort_unstable();
        rts.dedup();
        Ok(rts.into())
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
