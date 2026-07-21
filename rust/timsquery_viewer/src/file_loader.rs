use crate::domain::FileService;
use crate::error::ViewerError;
use egui_extras::{
    Table,
    TableBuilder,
};
use std::path::{
    Path,
    PathBuf,
};
use std::sync::Arc;
use timsquery::TimsElutionGroup;
use timsquery::ion::IonAnnot;
use timsquery::models::tolerance::Tolerance;
use timsquery::serde::IndexedPeaksHandle;
use timsquery::traits::QueryGeom;
use timsseek::data_sources::reference_library::ScoredIdentity;
use timsseek::scoring::pipeline::fill_scratch_from;
use timsseek::{
    ExpectedIntensities,
    ExpectedIntensity,
    RefQuery,
    ReferenceLibrary,
};
use tracing::{
    info,
    instrument,
};

/// Handles file dialogs and file loading operations
#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct FileLoader {
    pub elution_groups_path: Option<PathBuf>,
    pub raw_data_path: Option<PathBuf>,
    pub raw_data_url: Option<String>,
    pub tolerance_path: Option<PathBuf>,
}

impl FileLoader {
    pub fn new() -> Self {
        Self {
            elution_groups_path: None,
            raw_data_path: None,
            raw_data_url: None,
            tolerance_path: None,
        }
    }

    pub fn with_initial_paths(
        mut self,
        raw_data_path: &Option<PathBuf>,
        elution_groups_path: &Option<PathBuf>,
    ) -> Self {
        if let Some(raw_data_path) = raw_data_path {
            self.raw_data_path = Some(raw_data_path.clone());
        }
        if let Some(elution_groups_path) = elution_groups_path {
            self.elution_groups_path = Some(elution_groups_path.clone());
        }

        self
    }

    /// Open a file dialog for elution groups JSON file
    pub fn open_elution_groups_dialog(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter(
                "Spectral library (speclib/diann txt/tsv/parquet)",
                &["speclib", "txt", "tsv", "parquet", "json"],
            )
            .pick_file()
        {
            self.elution_groups_path = Some(path);
        }
    }

    /// Open a folder dialog for a raw data `.d` directory (or `.idx` cache).
    pub fn open_raw_data_dialog(&mut self) {
        if let Some(path) = rfd::FileDialog::new().pick_folder() {
            self.set_raw_data_path(path);
        }
    }

    /// Open a file dialog for a single raw data FILE (mzML / Thermo raw).
    /// `.d`/`.idx` are directories and need `open_raw_data_dialog`. The accepted
    /// extensions come from the reader registry — the single source of truth for
    /// which raw file formats are compiled in.
    pub fn open_raw_data_file_dialog(&mut self) {
        let extensions = timscentroid::reader::ReaderRegistry::with_builtins().file_extensions();
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("Raw MS file", &extensions)
            .pick_file()
        {
            self.set_raw_data_path(path);
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

    /// Load elution groups from a spectral library file
    pub fn load_elution_groups(&self, path: &Path) -> Result<ElutionGroupData, ViewerError> {
        FileService::load_elution_groups(path)
    }

    /// Load and index raw timsTOF data from a location (path or URL)
    pub fn load_raw_data_from_location(
        &self,
        location: &str,
    ) -> Result<Arc<IndexedPeaksHandle>, ViewerError> {
        FileService::load_raw_data_from_location(location)
    }

    /// Load tolerance settings from a JSON file
    pub fn load_tolerance(&self, path: &PathBuf) -> Result<Tolerance, ViewerError> {
        FileService::load_tolerance(path)
    }

    /// Set raw data URL for cloud storage
    pub fn set_raw_data_url(&mut self, url: String) {
        self.raw_data_url = Some(url);
        // Clear path when URL is set
        self.raw_data_path = None;
    }

    /// Set raw data path for local storage
    pub fn set_raw_data_path(&mut self, path: PathBuf) {
        self.raw_data_path = Some(path);
        // Clear URL when path is set
        self.raw_data_url = None;
    }

    /// Get the current raw data location (path or URL)
    pub fn get_raw_data_location(&self) -> Option<String> {
        self.raw_data_url
            .clone()
            .or_else(|| self.raw_data_path.as_ref().map(|p| p.display().to_string()))
    }

    /// Clear raw data location (both path and URL)
    pub fn clear_raw_data(&mut self) {
        self.raw_data_path = None;
        self.raw_data_url = None;
    }
}

/// The library the viewer displays and scores against.
///
/// This is the SAME columnar `ReferenceLibrary` arena the timsseek CLI scores
/// (loaded via `Speclib::from_file`). Iteration is via `item_at(flat)`
/// `RefQuery` flyweights: the geometry feeds `build_extraction` + `TraceScorer`
/// and the reference intensities + isotope envelope come from the flyweight's
/// `ExpectedIntensity` impl, which routes the envelope through
/// `isotope_dist_or_averagine` (averagine fallback). There is no private
/// isotope model here — the viewer shows exactly what the CLI scores.
#[derive(Debug)]
pub struct ElutionGroupData {
    inner: ReferenceLibrary,
}

const BASE_LABELS: [&str; 8] = [
    "ID",
    "Sequence",
    "Decoy",
    "RT (s)",
    "Mobility",
    "Precursor m/z",
    "Precursor Charge",
    "Fragments",
];

impl ElutionGroupData {
    pub fn new(inner: ReferenceLibrary) -> Self {
        Self { inner }
    }

    // No public `is_empty`: callers only ever ask for the row count.
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// The `RefQuery` flyweight at flat index `idx` (target/decoy expanded).
    fn item_at(&self, idx: usize) -> RefQuery<'_> {
        self.inner.item_at(idx)
    }

    /// Returns indices of all elution groups matching the ID filter.
    ///
    /// If filter is an empty string, returns ALL indices (no filtering applied).
    /// This allows seamless toggling between filtered and unfiltered views.
    #[instrument(skip(self, buffer))]
    pub fn matching_indices_for_id_filter(&self, filter: &str, buffer: &mut Vec<usize>) {
        buffer.clear();
        if filter.is_empty() {
            buffer.extend(0..self.len());
            return;
        }

        // Case-insensitive substring matching. to_lowercase() allocates per iteration,
        // but eq_ignore_ascii_case only supports full equality, not substring search.
        // Acceptable here since filtering runs only on search input changes, not per-frame.
        let filter_lower = filter.to_lowercase();
        let mut str_buffer = String::new();
        for i in 0..self.len() {
            self.key_onto(i, &mut str_buffer);
            if str_buffer.to_lowercase().contains(&filter_lower) {
                buffer.push(i);
            }
        }
    }

    /// Writes the filterable key (id + sequence) for `idx` into `buffer`.
    fn key_onto(&self, idx: usize, buffer: &mut String) {
        use std::fmt::Write;
        buffer.clear();
        let q = self.item_at(idx);
        let peptide = ScoredIdentity::materialize_peptide(&q);
        let _ = write!(buffer, "{}|{}|{}", q.id(), peptide.raw, !q.is_target());
    }

    pub fn get_elem(
        &self,
        index: usize,
    ) -> Result<(TimsElutionGroup<IonAnnot>, ExpectedIntensities<IonAnnot>), ViewerError> {
        if index >= self.len() {
            return Err(ViewerError::General(format!(
                "Elution group index {index} out of bounds"
            )));
        }

        // Materialize the scratch elution group + expected intensities from the
        // arena flyweight exactly the way the scoring pipeline does
        // (`fill_scratch_from`). The precursor envelope comes from
        // `expected_precursor_envelope()`, which routes through
        // `isotope_dist_or_averagine` — no `[1.0, 0.0, 0.0]` fallback.
        let q = self.item_at(index);
        let mut eg = TimsElutionGroup::empty_like();
        fill_scratch_from(&mut eg, &q);
        let expected = ExpectedIntensities::try_from_pairs(
            q.iter_expected_fragments(),
            q.expected_precursor_envelope().into_iter(),
        )
        .map_err(|e| ViewerError::General(format!("library entry {index}: {e:?}")))?;
        Ok((eg, expected))
    }

    pub fn render_table(
        &self,
        ui: &mut egui::Ui,
        filtered_eg_idxs: &[usize],
        selected_index: &mut Option<usize>,
        scroll_to_selection: bool,
    ) {
        let builder = TableBuilder::new(ui)
            .striped(true)
            .resizable(true)
            .cell_layout(egui::Layout::left_to_right(egui::Align::Center));
        let mut builder = self.add_columns(builder);
        if let Some(row_index) = selected_index.as_ref() {
            // Since the index is the original index, we need to find its position
            // in the filtered list first
            let local_index = match filtered_eg_idxs.binary_search(row_index) {
                Ok(idx) => idx,
                Err(insert_idx) => {
                    info!("Selected index {} not found in filtered indices", row_index);
                    // Set the selection to the closest match
                    let clamped_idx = if insert_idx >= filtered_eg_idxs.len() {
                        filtered_eg_idxs.len().saturating_sub(1) // Prevents underflow
                    } else {
                        insert_idx
                    };
                    if !filtered_eg_idxs.is_empty() {
                        *selected_index = Some(filtered_eg_idxs[clamped_idx]);
                    };
                    clamped_idx
                }
            };
            if scroll_to_selection && !filtered_eg_idxs.is_empty() {
                builder = builder.scroll_to_row(local_index, None);
            }
        }
        let builder = self.add_headers(builder);

        builder.body(|body| {
            let row_height = 18.0;
            body.rows(row_height, filtered_eg_idxs.len(), |mut row| {
                let row_idx = row.index();
                let original_idx = filtered_eg_idxs[row_idx];
                self.add_row_content(original_idx, selected_index, &mut row);
            });
        });
    }

    fn add_columns<'a>(&self, mut table: TableBuilder<'a>) -> TableBuilder<'a> {
        // Max column width ~40 chars at typical font size
        const MAX_COL_WIDTH: f32 = 280.0;
        for _ in BASE_LABELS.iter() {
            table = table.column(
                egui_extras::Column::auto()
                    .at_least(80.0)
                    .at_most(MAX_COL_WIDTH),
            );
        }
        table
    }

    fn add_headers<'a>(&self, builder: TableBuilder<'a>) -> Table<'a> {
        builder.header(20.0, |mut header| {
            for label in BASE_LABELS.iter() {
                header.col(|ui| {
                    ui.strong(*label);
                });
            }
        })
    }

    fn add_row_content(
        &self,
        idx: usize,
        selected_index: &mut Option<usize>,
        table_row: &mut egui_extras::TableRow,
    ) {
        let q = self.item_at(idx);
        let is_selected = Some(idx) == *selected_index;
        let peptide = ScoredIdentity::materialize_peptide(&q);

        let mut clicked = false;
        let mut add_col = |ui: &mut egui::Ui, text: &str| {
            let maybe_highlighted_text = if is_selected {
                egui::RichText::new(text).background_color(ui.visuals().selection.bg_fill)
            } else {
                egui::RichText::new(text)
            };
            let label = egui::Label::new(maybe_highlighted_text)
                .truncate()
                .sense(egui::Sense::click());
            let response = ui.add(label);
            let response = if is_selected {
                response.highlight()
            } else {
                response
            };
            if response.clicked() {
                clicked = true;
            }
        };

        table_row.col(|ui| {
            add_col(ui, &q.id().to_string());
        });
        table_row.col(|ui| {
            add_col(ui, &peptide.raw);
        });
        table_row.col(|ui| {
            add_col(ui, if q.is_target() { "No" } else { "Yes" });
        });
        table_row.col(|ui| {
            add_col(ui, &format!("{:.2}", q.rt_seconds()));
        });
        table_row.col(|ui| {
            add_col(ui, &format!("{:.4}", q.mobility_ook0()));
        });
        table_row.col(|ui| {
            add_col(ui, &format!("{:.4}", q.mono_precursor_mz()));
        });
        table_row.col(|ui| {
            add_col(ui, &format!("{}", q.precursor_charge()));
        });
        table_row.col(|ui| {
            add_col(ui, &format!("{}", q.fragment_count()));
        });

        if clicked {
            *selected_index = Some(idx);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use timsquery::models::QueryCollection;
    use timsquery::models::capabilities::LibCapabilities;
    use timsquery::utils::constants::PROTON_MASS;
    use timsseek::fragment_mass::isotope_dist_or_averagine;

    /// Build a one-entry `ReferenceLibrary` whose STRIPPED sequence has an
    /// uncountable composition (`B` is not a real residue), forcing the
    /// averagine isotope path.
    fn uncountable_lib() -> ReferenceLibrary {
        let mut geom = QueryCollection::with_capabilities(LibCapabilities::default_diann());
        geom.push_target(
            600.0,
            1,
            10.0,
            1.0,
            &[
                (IonAnnot::try_from("y3").unwrap(), 300.0),
                (IonAnnot::try_from("y5").unwrap(), 500.0),
            ],
            "PEPBK",
            "PEPBK",
            &[],
        );
        geom.seal();
        ReferenceLibrary {
            geom,
            frag_intens: vec![1.0, 0.5],
        }
    }

    /// The envelope the viewer DISPLAYS (obtained via `get_elem`, i.e. the
    /// flyweight's `expected_precursor_envelope`) must equal what the scoring
    /// path computes via `isotope_dist_or_averagine` — NOT the deleted
    /// `[1.0, 0.0, 0.0]` fallback.
    #[test]
    fn displayed_envelope_matches_isotope_dist_or_averagine() {
        let data = ElutionGroupData::new(uncountable_lib());
        let (_eg, expected) = data.get_elem(0).unwrap();

        let charge = 1.0_f64;
        let neutral = 600.0 * charge - charge * PROTON_MASS;
        let (_src, env) = isotope_dist_or_averagine("PEPBK", neutral);

        // Every displayed precursor isotope equals the shared model output.
        for (iso_idx, ref_intensity) in env.iter().enumerate() {
            let displayed = expected
                .get_precursor(iso_idx as i8)
                .unwrap_or_else(|| panic!("missing displayed precursor isotope {iso_idx}"));
            assert!(
                (displayed - ref_intensity).abs() < 1e-6,
                "isotope {iso_idx}: displayed {displayed} != shared model {ref_intensity}",
            );
        }

        // Prove the [1,0,0] fallback is gone: an averagine peptide of this mass
        // has non-zero +1/+2 isotope peaks.
        assert!(
            expected.get_precursor(1).unwrap() > 0.0,
            "the +1 isotope must be populated by the averagine model (not the [1,0,0] fallback)",
        );
    }
}
