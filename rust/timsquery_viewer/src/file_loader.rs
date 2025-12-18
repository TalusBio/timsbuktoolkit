use std::path::PathBuf;
use std::sync::Arc;
use timscentroid::IndexedTimstofPeaks;
use timsquery::models::tolerance::Tolerance;
use timsquery::serde::ElutionGroupCollection;
use std::collections::HashMap;

use crate::domain::FileService;
use crate::error::ViewerError;

/// Handles file dialogs and file loading operations
#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct FileLoader {
    pub elution_groups_path: Option<PathBuf>,
    pub raw_data_path: Option<PathBuf>,
    pub tolerance_path: Option<PathBuf>,
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
            .add_filter(
                "Elution Groups File (json/diann txt/tsv)",
                &["json", "txt", "tsv"],
            )
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
        FileService::load_elution_groups(path)
    }

    /// Load and index raw timsTOF data
    pub fn load_raw_data(
        &self,
        path: &PathBuf,
    ) -> Result<(Arc<IndexedTimstofPeaks>, Arc<[u32]>), ViewerError> {
        FileService::load_raw_data(path)
    }

    /// Load tolerance settings from a JSON file
    pub fn load_tolerance(&self, path: &PathBuf) -> Result<Tolerance, ViewerError> {
        FileService::load_tolerance(path)
    }
}

/// Wrapper around elution group collection with optional library metadata
#[derive(Debug)]
pub struct ElutionGroupData {
    /// The parsed elution groups
    pub inner: ElutionGroupCollection,
    /// Library fragment intensities from library sidecar file.
    /// Maps elution group ID → list of (fragment_label, relative_intensity).
    /// Used for mirror plot visualization comparing observed vs predicted spectra.
    pub extras: Option<HashMap<u64, Vec<(String, f32)>>>,
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

/// Execute a macro with the appropriate elution group collection variant.
///
/// # Example
/// ```ignore
/// macro_rules! process {
///     ($egs:expr) => {{
///         ChromatogramService::generate(&$egs[idx], ...)
///     }};
/// }
/// with_elution_collection!(elution_groups, process)
/// ```
#[macro_export]
macro_rules! with_elution_collection {
    ($data:expr, $macro_name:ident) => {
        match &$data.inner {
            timsquery::serde::ElutionGroupCollection::StringLabels(egs) => $macro_name!(egs),
            timsquery::serde::ElutionGroupCollection::MzpafLabels(egs) => $macro_name!(egs),
            timsquery::serde::ElutionGroupCollection::TinyIntLabels(egs) => $macro_name!(egs),
            timsquery::serde::ElutionGroupCollection::IntLabels(egs) => $macro_name!(egs),
        }
    };
}
