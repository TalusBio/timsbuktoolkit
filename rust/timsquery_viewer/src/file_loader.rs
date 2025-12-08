use std::path::PathBuf;
use std::sync::Arc;
use timscentroid::{
    IndexedTimstofPeaks,
    TimsTofPath,
};
use timsquery::models::tolerance::Tolerance;
use timsquery::serde::{
    ElutionGroupCollection,
    load_index_caching,
};
use timsrust::MSLevel;

use tracing::info;

use crate::error::ViewerError;

/// Handles file dialogs and file loading operations
pub struct FileLoader {
    pub elution_groups_path: Option<PathBuf>,
    pub raw_data_path: Option<PathBuf>,
    pub tolerance_path: Option<PathBuf>,
}

#[derive(Debug)]
pub struct ElutionGroupData {
    pub inner: ElutionGroupCollection,
}

impl ElutionGroupData {
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns indices of all elution groups matching the ID filter.
    ///
    /// If filter is an empty string, returns ALL indices (no filtering applied).
    /// This allows seamless toggling between filtered and unfiltered views.
    pub fn matching_indices_for_id_filter(&self, filter: &str) -> Vec<usize> {
        if filter.is_empty() {
            return (0..self.len()).collect();
        }

        macro_rules! get_ids {
            ($self:expr) => {
                $self
                    .iter()
                    .enumerate()
                    .filter(|(_, eg)| eg.id().to_string().contains(filter))
                    .map(|(idx, _)| idx)
                    .collect()
            };
        }

        match &self.inner {
            ElutionGroupCollection::StringLabels(egs) => get_ids!(egs),
            ElutionGroupCollection::MzpafLabels(egs) => get_ids!(egs),
            ElutionGroupCollection::TinyIntLabels(egs) => get_ids!(egs),
            ElutionGroupCollection::IntLabels(egs) => get_ids!(egs),
        }
    }
}

impl FileLoader {
    pub fn new() -> Self {
        Self {
            elution_groups_path: None,
            raw_data_path: None,
            tolerance_path: None,
        }
    }

    /// Open a file dialog for elution groups JSON file
    pub fn open_elution_groups_dialog(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("Elution Groups File (json/diann txt)", &["json", "txt"])
            .pick_file()
        {
            self.elution_groups_path = Some(path);
        }
    }

    /// Open a file dialog for raw data .d directory
    pub fn open_raw_data_dialog(&mut self) {
        if let Some(path) = rfd::FileDialog::new().pick_folder() {
            self.raw_data_path = Some(path);
        }
    }

    /// Open a file dialog for tolerance settings JSON file
    pub fn open_tolerance_dialog(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("JSON", &["json"])
            .pick_file()
        {
            self.tolerance_path = Some(path);
        }
    }

    /// Load elution groups from a JSON file
    pub fn load_elution_groups(&self, path: &PathBuf) -> Result<ElutionGroupData, ViewerError> {
        let res = timsquery::serde::read_library_file(path)?;
        info!(
            "Loaded {} elution groups from {}",
            res.len(),
            path.display()
        );
        Ok(ElutionGroupData { inner: res })
    }

    /// Load and index raw timsTOF data
    pub fn load_raw_data(
        &self,
        path: &PathBuf,
    ) -> Result<(Arc<IndexedTimstofPeaks>, Arc<[u32]>), ViewerError> {
        let index = load_index_caching(path).map_err(|e| ViewerError::DataLoading {
            path: path.clone(),
            source: Box::new(ViewerError::General(format!("{:?}", e))),
        })?;

        let rts = get_ms1_rts_as_millis(path)?;

        Ok((Arc::new(index), rts))
    }

    /// Load tolerance settings from a JSON file
    pub fn load_tolerance(&self, path: &PathBuf) -> Result<Tolerance, ViewerError> {
        let file_content = std::fs::read_to_string(path)?;
        let tolerance: Tolerance = serde_json::from_str(&file_content)?;
        Ok(tolerance)
    }
}

/// Retrieves MS1 retention times from a TIMS-TOF file, sorted and deduped
fn get_ms1_rts_as_millis(file: &PathBuf) -> Result<Arc<[u32]>, ViewerError> {
    let ttp = TimsTofPath::new(file).map_err(|e| ViewerError::TimsFileLoad {
        path: file.clone(),
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
