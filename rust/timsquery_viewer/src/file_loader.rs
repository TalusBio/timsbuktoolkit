use std::path::PathBuf;
use std::sync::Arc;
use timscentroid::{
    IndexedTimstofPeaks,
    TimsTofPath,
};
use timsquery::ion::IonAnnot;
use timsquery::models::elution_group::TimsElutionGroup;
use timsquery::models::tolerance::Tolerance;
use timsquery::serde::load_index_caching;
use timsquery::tinyvec::TinyVec;
use timsrust::MSLevel;

use tracing::warn;

use crate::error::ViewerError;

/// Handles file dialogs and file loading operations
pub struct FileLoader {
    pub elution_groups_path: Option<PathBuf>,
    pub raw_data_path: Option<PathBuf>,
    pub tolerance_path: Option<PathBuf>,
}

#[derive(Debug)]
pub enum ElutionGroupData {
    StringLabels(Vec<TimsElutionGroup<String>>),
    MzpafLabels(Vec<TimsElutionGroup<IonAnnot>>),
}

impl ElutionGroupData {
    pub fn len(&self) -> usize {
        match self {
            ElutionGroupData::StringLabels(egs) => egs.len(),
            ElutionGroupData::MzpafLabels(egs) => egs.len(),
        }
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

        match self {
            ElutionGroupData::StringLabels(egs) => egs
                .iter()
                .enumerate()
                .filter(|(_, eg)| eg.id().to_string().contains(filter))
                .map(|(idx, _)| idx)
                .collect(),
            ElutionGroupData::MzpafLabels(egs) => egs
                .iter()
                .enumerate()
                .filter(|(_, eg)| eg.id().to_string().contains(filter))
                .map(|(idx, _)| idx)
                .collect(),
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
            .add_filter("JSON", &["json"])
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
        let file_content = std::fs::read_to_string(path)?;

        if let Ok(eg_inputs) = serde_json::from_str::<Vec<ElutionGroupInput>>(&file_content) {
            let out: Result<Vec<TimsElutionGroup<IonAnnot>>, ViewerError> =
                eg_inputs.into_iter().map(|x| x.try_into()).collect();
            return Ok(ElutionGroupData::MzpafLabels(out?));
        }

        if let Ok(egs) = serde_json::from_str::<Vec<TimsElutionGroup<IonAnnot>>>(&file_content) {
            return Ok(ElutionGroupData::MzpafLabels(egs));
        }

        let egs_string: Vec<TimsElutionGroup<String>> = serde_json::from_str(&file_content)?;
        warn!(
            "Elution groups contained fragment labels as strings that are not interpretable as mzpaf, this can cause performance degradation."
        );
        Ok(ElutionGroupData::StringLabels(egs_string))
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

/// User-friendly format for specifying elution groups in an input file
/// (copied from timsquery_cli for compatibility)
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ElutionGroupInput {
    pub id: u64,
    pub mobility: f32,
    pub rt_seconds: f32,
    pub precursors: Vec<f64>,
    pub fragments: Vec<f64>,
    pub fragment_labels: Option<Vec<IonAnnot>>,
}

impl TryFrom<ElutionGroupInput> for TimsElutionGroup<IonAnnot> {
    type Error = ViewerError;

    fn try_from(val: ElutionGroupInput) -> Result<Self, Self::Error> {
        let builder = TimsElutionGroup::builder()
            .id(val.id)
            .mobility_ook0(val.mobility)
            .rt_seconds(val.rt_seconds)
            .precursor_labels(
                val.precursors
                    .iter()
                    .enumerate()
                    .map(|(i, _)| i as i8)
                    .collect(),
            )
            .precursor_mzs(val.precursors.clone());

        let num_fragments = val.fragments.len();
        let builder = builder.fragment_mzs(val.fragments.clone());
        let builder = match val.fragment_labels {
            Some(ref labels) => {
                if labels.len() != num_fragments {
                    return Err(ViewerError::General(format!(
                        "Fragment labels length {} does not match fragments length {}",
                        labels.len(),
                        num_fragments
                    )));
                }
                builder.fragment_labels(TinyVec::Heap(labels.clone()))
            }
            None => builder.fragment_labels(
                (0..num_fragments)
                    .map(|i| IonAnnot::try_new('?', Some(i as u8), 1, 1).unwrap())
                    .collect(),
            ),
        };

        Ok(builder.try_build().expect("I checked the sizes!"))
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
