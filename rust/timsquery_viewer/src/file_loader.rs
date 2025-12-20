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
use timscentroid::IndexedTimstofPeaks;
use timsquery::models::tolerance::Tolerance;
use timsquery::serde::{
    DiannPrecursorExtras,
    ElutionGroupCollection,
    FileReadingExtras,
};
use timsquery::{
    KeyLike,
    TimsElutionGroup,
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
    pub fn load_elution_groups(&self, path: &Path) -> Result<ElutionGroupData, ViewerError> {
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

const DIANN_EXTRA_LABELS: [&str; 3] = ["Modified Peptide", "Protein ID(s)", "Is Decoy"];

impl ElutionGroupData {
    pub fn new(inner: ElutionGroupCollection) -> Self {
        Self { inner }
    }

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
    #[instrument(skip(self, buffer))]
    pub fn matching_indices_for_id_filter(&self, filter: &str, buffer: &mut Vec<usize>) {
        buffer.clear();
        if filter.is_empty() {
            buffer.extend(0..self.len());
            return;
        }

        let mut str_buffer = String::new();
        for i in 0..self.len() {
            if self.key_onto(i, &mut str_buffer).is_ok() && str_buffer.contains(filter) {
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
        let extras = match &self.inner {
            ElutionGroupCollection::StringLabels(_, extras)
            | ElutionGroupCollection::MzpafLabels(_, extras)
            | ElutionGroupCollection::TinyIntLabels(_, extras)
            | ElutionGroupCollection::IntLabels(_, extras) => match extras {
                Some(FileReadingExtras::Diann(diann_extras)) => Some(&diann_extras[idx]),
                _ => None,
            },
        };
        if let Some(diann_extra) = extras {
            write!(
                buffer,
                "|{}|{}|{}",
                diann_extra.modified_peptide, diann_extra.protein_id, diann_extra.is_decoy
            )
            .map_err(|_| ())?;
        }
        Ok(())
    }

    pub fn get_elem(
        &self,
        index: usize,
    ) -> (
        Option<TimsElutionGroup<String>>,
        Option<HashMap<String, f32>>,
    ) {
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

        let extra = match extras {
            Some(FileReadingExtras::Diann(diann_extras)) => diann_extras.get(index).map(|de| {
                HashMap::from_iter(
                    de.relative_intensities
                        .iter()
                        .cloned()
                        .map(|(k, v)| (k.clone().to_string(), v)),
                )
            }),
            _ => None,
        };
        (eg, extra)
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
            let local_index = filtered_eg_idxs.binary_search(&row_index);
            if let Ok(local_idx) = local_index {
                // we want to scroll to this row only once
                if scroll_to_selection {
                    builder = builder.scroll_to_row(local_idx, None);
                }
            } else {
                info!("Selected index {} not found in filtered indices", row_index);
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
        let has_diann_extras = match &self.inner {
            ElutionGroupCollection::StringLabels(_, extras)
            | ElutionGroupCollection::MzpafLabels(_, extras)
            | ElutionGroupCollection::TinyIntLabels(_, extras)
            | ElutionGroupCollection::IntLabels(_, extras) => match extras {
                Some(FileReadingExtras::Diann(_)) => true,
                _ => false,
            },
        };
        if has_diann_extras {
            for _ in DIANN_EXTRA_LABELS.iter() {
                table = table.column(egui_extras::Column::auto().at_least(100.0));
            }
        }
        for _ in BASE_LABELS.iter() {
            table = table.column(egui_extras::Column::auto().at_least(80.0));
        }
        table
    }

    fn add_headers<'a>(&self, builder: TableBuilder<'a>) -> Table<'a> {
        let has_diann_extras = match &self.inner {
            ElutionGroupCollection::StringLabels(_, extras)
            | ElutionGroupCollection::MzpafLabels(_, extras)
            | ElutionGroupCollection::TinyIntLabels(_, extras)
            | ElutionGroupCollection::IntLabels(_, extras) => match extras {
                Some(FileReadingExtras::Diann(_)) => true,
                _ => false,
            },
        };

        builder.header(20.0, |mut header| {
            if has_diann_extras {
                for label in DIANN_EXTRA_LABELS.iter() {
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

    fn get_column_strings(&self, index: usize) -> Vec<String> {
        match &self.inner {
            ElutionGroupCollection::StringLabels(egs, _) => {
                egs.iter().map(|eg| eg.id().to_string()).collect()
            }
            ElutionGroupCollection::MzpafLabels(egs, _) => {
                egs.iter().map(|eg| eg.id().to_string()).collect()
            }
            ElutionGroupCollection::TinyIntLabels(egs, _) => {
                egs.iter().map(|eg| eg.id().to_string()).collect()
            }
            ElutionGroupCollection::IntLabels(egs, _) => {
                egs.iter().map(|eg| eg.id().to_string()).collect()
            }
        }
    }

    /// Helper function to add row content
    /// `is_selected` indicates if the row is currently selected
    /// This function adds the appropriate columns based on available extras
    /// Returns true if any of the content was clicked (for selection handling)
    fn add_row_content_inner<T: KeyLike>(
        eg: &TimsElutionGroup<T>,
        extras: Option<DiannPrecursorExtras>,
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
            let label = ui.selectable_label(is_selected, maybe_highlighted_text);
            let label = if is_selected {
                label.highlight()
            } else {
                label
            };

            if label.clicked() {
                clicked = true;
            }
        };
        match extras {
            Some(diann_extra) => {
                table_row.col(|ui| {
                    add_col(ui, &diann_extra.modified_peptide);
                });
                table_row.col(|ui| {
                    add_col(ui, &diann_extra.protein_id);
                });
                table_row.col(|ui| {
                    add_col(ui, if diann_extra.is_decoy { "Yes" } else { "No" });
                });
            }
            None => { /* No extra columns */ }
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
        let diann_extra = match &self.inner {
            ElutionGroupCollection::StringLabels(_, extras)
            | ElutionGroupCollection::MzpafLabels(_, extras)
            | ElutionGroupCollection::TinyIntLabels(_, extras)
            | ElutionGroupCollection::IntLabels(_, extras) => match extras {
                Some(FileReadingExtras::Diann(diann_extras)) => Some(diann_extras[idx].clone()),
                _ => None,
            },
        };
        let is_selected = Some(idx) == *selected_index;
        let clicked = match &self.inner {
            ElutionGroupCollection::StringLabels(egs, _) => {
                Self::add_row_content_inner(&egs[idx], diann_extra, table_row, is_selected)
            }
            ElutionGroupCollection::MzpafLabels(egs, _) => {
                Self::add_row_content_inner(&egs[idx], diann_extra, table_row, is_selected)
            }
            ElutionGroupCollection::TinyIntLabels(egs, _) => {
                Self::add_row_content_inner(&egs[idx], diann_extra, table_row, is_selected)
            }
            ElutionGroupCollection::IntLabels(egs, _) => {
                Self::add_row_content_inner(&egs[idx], diann_extra, table_row, is_selected)
            }
        };
        if clicked {
            *selected_index = Some(idx);
        }
    }
}
