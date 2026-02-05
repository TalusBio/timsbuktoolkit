use crate::domain::FileService;
use crate::error::ViewerError;
use egui_extras::{
    Table,
    TableBuilder,
};
use std::collections::HashMap;
use std::path::{
    Path,
    PathBuf,
};
use std::sync::Arc;
use timsquery::ion::IonAnnot;
use timsquery::models::tolerance::Tolerance;
use timsquery::serde::{
    ElutionGroupCollection,
    FileReadingExtras,
    IndexedPeaksHandle,
};
use timsquery::{
    KeyLike,
    TimsElutionGroup,
};
use timsseek::ExpectedIntensities;
use timsseek::fragment_mass::elution_group_converter::isotope_dist_from_seq;
use tracing::{
    info,
    instrument,
    warn,
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
                "Elution Groups File (json/diann txt/tsv)",
                &["json", "txt", "tsv", "parquet"],
            )
            .pick_file()
        {
            self.elution_groups_path = Some(path);
        }
    }

    /// Open a file dialog for raw data .d directory
    pub fn open_raw_data_dialog(&mut self) {
        if let Some(path) = rfd::FileDialog::new().pick_folder() {
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

    /// Load elution groups from a JSON file
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

#[derive(Debug)]
pub struct ElutionGroupData {
    inner: ElutionGroupCollection,
}
const BASE_LABELS: [&str; 6] = [
    "ID",
    "RT (s)",
    "Mobility",
    "Precursor m/z",
    "Precursor Charge",
    "Fragments",
];

/// Labels for library-specific extra columns (DIA-NN and Spectronaut share the same structure)
const LIBRARY_EXTRA_LABELS: [&str; 3] = ["Modified Peptide", "Protein ID(s)", "Is Decoy"];

/// Common view structure for library extras that both DIA-NN and Spectronaut can map to
struct LibraryExtrasView {
    modified_peptide: String,
    protein_id: String,
    is_decoy: bool,
}

impl ElutionGroupData {
    pub fn new(inner: ElutionGroupCollection) -> Self {
        Self { inner }
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }

    #[must_use]
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
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
            if self.key_onto(i, &mut str_buffer).is_ok()
                && str_buffer.to_lowercase().contains(&filter_lower)
            {
                buffer.push(i);
            }
        }
    }

    /// Adds the key contents to the string, the idea here is to avoid allocations
    fn key_onto(&self, idx: usize, buffer: &mut String) -> Result<(), ()> {
        use std::fmt::Write;
        buffer.clear();
        match &self.inner {
            ElutionGroupCollection::StringLabels(egs, _) => {
                write!(buffer, "{}", egs[idx].id()).map_err(|_| ())?;
            }
            ElutionGroupCollection::MzpafLabels(egs, _) => {
                write!(buffer, "{}", egs[idx].id()).map_err(|_| ())?;
            }
            ElutionGroupCollection::TinyIntLabels(egs, _) => {
                write!(buffer, "{}", egs[idx].id()).map_err(|_| ())?;
            }
            ElutionGroupCollection::IntLabels(egs, _) => {
                write!(buffer, "{}", egs[idx].id()).map_err(|_| ())?;
            }
        }
        let extras_view = self.get_library_extras_view(idx);
        if let Some(extra) = extras_view {
            write!(
                buffer,
                "|{}|{}|{}",
                extra.modified_peptide, extra.protein_id, extra.is_decoy
            )
            .map_err(|_| ())?;
        }
        Ok(())
    }

    /// Extracts a common view of library extras from either DIA-NN or Spectronaut format
    fn get_library_extras_view(&self, idx: usize) -> Option<LibraryExtrasView> {
        let extras = match &self.inner {
            ElutionGroupCollection::StringLabels(_, extras)
            | ElutionGroupCollection::MzpafLabels(_, extras)
            | ElutionGroupCollection::TinyIntLabels(_, extras)
            | ElutionGroupCollection::IntLabels(_, extras) => extras.as_ref()?,
        };
        match extras {
            FileReadingExtras::Diann(diann_extras) => {
                let de = diann_extras.get(idx)?;
                Some(LibraryExtrasView {
                    modified_peptide: de.modified_peptide.clone(),
                    protein_id: de.protein_id.clone(),
                    is_decoy: de.is_decoy,
                })
            }
            FileReadingExtras::Spectronaut(spectronaut_extras) => {
                let se = spectronaut_extras.get(idx)?;
                Some(LibraryExtrasView {
                    modified_peptide: se.modified_peptide.clone(),
                    protein_id: se.protein_id.clone(),
                    is_decoy: se.is_decoy,
                })
            }
        }
    }

    /// Checks if the collection has library extras (DIA-NN or Spectronaut)
    fn has_library_extras(&self) -> bool {
        match &self.inner {
            ElutionGroupCollection::StringLabels(_, extras)
            | ElutionGroupCollection::MzpafLabels(_, extras)
            | ElutionGroupCollection::TinyIntLabels(_, extras)
            | ElutionGroupCollection::IntLabels(_, extras) => extras.is_some(),
        }
    }

    /// Builds ExpectedIntensities from library extras (shared by DIA-NN and Spectronaut)
    fn build_expected_intensities(
        stripped_peptide: &str,
        relative_intensities: &[(IonAnnot, f32)],
        eg: &mut TimsElutionGroup<String>,
    ) -> ExpectedIntensities<String> {
        let fragment_intensities = HashMap::from_iter(
            relative_intensities
                .iter()
                .cloned()
                .map(|(k, v)| (k.to_string(), v)),
        );

        let isotopes = match isotope_dist_from_seq(stripped_peptide) {
            Ok(isotopes) => isotopes,
            Err(e) => {
                warn!(
                    "Failed to calculate isotope distribution for sequence {}: {}",
                    stripped_peptide, e
                );
                [1.0, 0.0, 0.0]
            }
        };

        eg.set_precursor_labels([0, 1, 2].iter().cloned());
        let precursor_intensities: HashMap<i8, f32> = isotopes
            .iter()
            .cloned()
            .enumerate()
            .map(|(i, intensity)| (i as i8, intensity))
            .collect();

        ExpectedIntensities {
            precursor_intensities,
            fragment_intensities,
        }
    }

    pub fn get_elem(
        &self,
        index: usize,
    ) -> Result<(TimsElutionGroup<String>, ExpectedIntensities<String>), ViewerError> {
        let (eg, extras) = match &self.inner {
            ElutionGroupCollection::StringLabels(egs, ext) => (egs.get(index).cloned(), ext),
            ElutionGroupCollection::MzpafLabels(egs, ext) => {
                (egs.get(index).map(|eg| eg.cast(|x| x.to_string())), ext)
            }
            ElutionGroupCollection::TinyIntLabels(egs, ext) => {
                (egs.get(index).map(|eg| eg.cast(|x| x.to_string())), ext)
            }
            ElutionGroupCollection::IntLabels(egs, ext) => {
                (egs.get(index).map(|eg| eg.cast(|x| x.to_string())), ext)
            }
        };
        let mut eg = eg.ok_or(ViewerError::General(format!(
            "Elution group index {} out of bounds",
            index
        )))?;

        let extra = match extras {
            Some(FileReadingExtras::Diann(diann_extras)) => {
                let de = diann_extras.get(index).ok_or(ViewerError::General(format!(
                    "Diann extras index {} out of bounds",
                    index
                )))?;
                Self::build_expected_intensities(
                    &de.stripped_peptide,
                    &de.relative_intensities,
                    &mut eg,
                )
            }
            Some(FileReadingExtras::Spectronaut(spectronaut_extras)) => {
                let se = spectronaut_extras
                    .get(index)
                    .ok_or(ViewerError::General(format!(
                        "Spectronaut extras index {} out of bounds",
                        index
                    )))?;
                Self::build_expected_intensities(
                    &se.stripped_peptide,
                    &se.relative_intensities,
                    &mut eg,
                )
            }
            None => ExpectedIntensities {
                precursor_intensities: eg.iter_precursors().map(|(idx, _mz)| (idx, 1.0)).collect(),
                fragment_intensities: eg
                    .iter_fragments()
                    .map(|(label, _mz)| (label.to_string(), 1.0))
                    .collect(),
            },
        };
        Ok((eg, extra))
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
        if self.has_library_extras() {
            for _ in LIBRARY_EXTRA_LABELS.iter() {
                table = table.column(
                    egui_extras::Column::auto()
                        .at_least(100.0)
                        .at_most(MAX_COL_WIDTH),
                );
            }
        }
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
            if self.has_library_extras() {
                for label in LIBRARY_EXTRA_LABELS.iter() {
                    header.col(|ui| {
                        ui.strong(*label);
                    });
                }
            }
            for label in BASE_LABELS.iter() {
                header.col(|ui| {
                    ui.strong(*label);
                });
            }
        })
    }

    /// Helper function to add row content
    /// `is_selected` indicates if the row is currently selected
    /// This function adds the appropriate columns based on available extras
    /// Returns true if any of the content was clicked (for selection handling)
    fn add_row_content_inner<T: KeyLike>(
        eg: &TimsElutionGroup<T>,
        extras: Option<LibraryExtrasView>,
        table_row: &mut egui_extras::TableRow,
        is_selected: bool,
    ) -> bool {
        let mut clicked = false;
        let mut add_col = |ui: &mut egui::Ui, text: &str| {
            // Highlight if selected
            let maybe_highlighted_text = if is_selected {
                egui::RichText::new(text).background_color(ui.visuals().selection.bg_fill)
            } else {
                egui::RichText::new(text)
            };
            // Use Label with truncate to show "..." when text overflows column width
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
        if let Some(extra) = extras {
            table_row.col(|ui| {
                add_col(ui, &extra.modified_peptide);
            });
            table_row.col(|ui| {
                add_col(ui, &extra.protein_id);
            });
            table_row.col(|ui| {
                add_col(ui, if extra.is_decoy { "Yes" } else { "No" });
            });
        }
        table_row.col(|ui| {
            add_col(ui, &eg.id().to_string());
        });
        table_row.col(|ui| {
            let text = format!("{:.2}", eg.rt_seconds());
            add_col(ui, &text);
        });
        table_row.col(|ui| {
            let text = format!("{:.4}", eg.mobility_ook0());
            add_col(ui, &text);
        });
        table_row.col(|ui| {
            let display_text = format!("{:.4}", eg.mono_precursor_mz());
            add_col(ui, &display_text);
        });
        table_row.col(|ui| {
            let text = format!("{}", eg.precursor_charge());
            add_col(ui, &text);
        });
        table_row.col(|ui| {
            let text = format!("{}", eg.fragment_count());
            add_col(ui, &text);
        });
        clicked
    }

    fn add_row_content(
        &self,
        idx: usize,
        selected_index: &mut Option<usize>,
        table_row: &mut egui_extras::TableRow,
    ) {
        let library_extras = self.get_library_extras_view(idx);
        let is_selected = Some(idx) == *selected_index;
        let clicked = match &self.inner {
            ElutionGroupCollection::StringLabels(egs, _) => {
                Self::add_row_content_inner(&egs[idx], library_extras, table_row, is_selected)
            }
            ElutionGroupCollection::MzpafLabels(egs, _) => {
                Self::add_row_content_inner(&egs[idx], library_extras, table_row, is_selected)
            }
            ElutionGroupCollection::TinyIntLabels(egs, _) => {
                Self::add_row_content_inner(&egs[idx], library_extras, table_row, is_selected)
            }
            ElutionGroupCollection::IntLabels(egs, _) => {
                Self::add_row_content_inner(&egs[idx], library_extras, table_row, is_selected)
            }
        };
        if clicked {
            *selected_index = Some(idx);
        }
    }
}
