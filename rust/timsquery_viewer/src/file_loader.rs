use std::path::PathBuf;
use std::sync::Arc;
use timscentroid::{
    IndexedTimstofPeaks,
    TimsTofPath,
};
use timsquery::models::elution_group::ElutionGroup;
use timsquery::models::tolerance::Tolerance;
use timsquery::serde::load_index_caching;
use timsrust::MSLevel;

use crate::error::ViewerError;

/// Handles file dialogs and file loading operations
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
    pub fn load_elution_groups(
        &self,
        path: &PathBuf,
    ) -> Result<Vec<ElutionGroup<usize>>, ViewerError> {
        let file_content = std::fs::read_to_string(path)?;

        // Try parsing as Vec<ElutionGroupInput> first
        if let Ok(eg_inputs) = serde_json::from_str::<Vec<ElutionGroupInput>>(&file_content) {
            let out: Vec<ElutionGroup<usize>> = eg_inputs.into_iter().map(|x| x.into()).collect();
            return Ok(out);
        }

        // Try parsing as Vec<ElutionGroup<usize>>
        if let Ok(egs) = serde_json::from_str::<Vec<ElutionGroup<usize>>>(&file_content) {
            return Ok(egs);
        }

        // Try parsing as Vec<ElutionGroup<String>> and convert
        let egs_string: Vec<ElutionGroup<String>> = serde_json::from_str(&file_content)?;
        let mut out: Vec<ElutionGroup<usize>> = Vec::with_capacity(egs_string.len());
        for (i, eg) in egs_string.into_iter().enumerate() {
            let eg_usize = ElutionGroup {
                id: i as u64,
                mobility: eg.mobility,
                rt_seconds: eg.rt_seconds,
                precursors: eg.precursors,
                fragments: Arc::from(
                    eg.fragments
                        .iter()
                        .map(|(_lab, mz)| *mz)
                        .enumerate()
                        .collect::<Vec<(usize, f64)>>(),
                ),
            };
            out.push(eg_usize);
        }
        Ok(out)
    }

    /// Load and index raw timsTOF data
    pub fn load_raw_data(
        &self,
        path: &PathBuf,
    ) -> Result<(Arc<IndexedTimstofPeaks>, Arc<[u32]>), ViewerError> {
        // Load and index the data
        let index = load_index_caching(path).map_err(|e| ViewerError::DataLoading {
            path: path.clone(),
            source: Box::new(ViewerError::General(format!("{:?}", e))),
        })?;

        // Extract MS1 retention times
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
}

impl From<ElutionGroupInput> for ElutionGroup<usize> {
    fn from(val: ElutionGroupInput) -> Self {
        let precursors: Arc<[(i8, f64)]> = Arc::from(
            val.precursors
                .into_iter()
                .enumerate()
                .map(|(i, mz)| (i as i8, mz))
                .collect::<Vec<(i8, f64)>>(),
        );
        let fragments: Arc<[(usize, f64)]> = Arc::from(
            val.fragments
                .into_iter()
                .enumerate()
                .collect::<Vec<(usize, f64)>>(),
        );
        ElutionGroup {
            id: val.id,
            mobility: val.mobility,
            rt_seconds: val.rt_seconds,
            precursors,
            fragments,
        }
    }
}

/// Retrieves MS1 retention times from a TIMS-TOF file, sorted and deduped
/// (copied from timsquery_cli)
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
