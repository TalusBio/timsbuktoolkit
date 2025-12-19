use crate::domain::FileService;
use crate::error::ViewerError;
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;
use timscentroid::IndexedTimstofPeaks;
use timsquery::models::tolerance::Tolerance;
use timsquery::serde::{
    ElutionGroupCollection,
    FileReadingExtras,
};
use timsquery::{
    KeyLike,
    TimsElutionGroup,
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

#[derive(Debug)]
pub struct ElutionGroupData {
    inner: ElutionGroupCollection,
}

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
            ElutionGroupCollection::StringLabels(egs, _) => get_ids!(egs),
            ElutionGroupCollection::MzpafLabels(egs, _) => get_ids!(egs),
            ElutionGroupCollection::TinyIntLabels(egs, _) => get_ids!(egs),
            ElutionGroupCollection::IntLabels(egs, _) => get_ids!(egs),
        }
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
    ) {
        match &self.inner {
            ElutionGroupCollection::StringLabels(egs, _) => {
                render_precursor_table_filtered(ui, filtered_eg_idxs, egs, selected_index)
            }
            ElutionGroupCollection::MzpafLabels(egs, _) => {
                render_precursor_table_filtered(ui, filtered_eg_idxs, egs, selected_index)
            }
            ElutionGroupCollection::TinyIntLabels(egs, _) => {
                render_precursor_table_filtered(ui, filtered_eg_idxs, egs, selected_index)
            }
            ElutionGroupCollection::IntLabels(egs, _) => {
                render_precursor_table_filtered(ui, filtered_eg_idxs, egs, selected_index)
            }
        }
    }
}

fn render_precursor_table_filtered<T: KeyLike>(
    ui: &mut egui::Ui,
    filtered_eg_idxs: &[usize],
    reference_eg_slice: &[TimsElutionGroup<T>],
    selected_index: &mut Option<usize>,
) {
    use egui_extras::{
        Column,
        TableBuilder,
    };

    TableBuilder::new(ui)
        .striped(true)
        .resizable(true)
        .cell_layout(egui::Layout::left_to_right(egui::Align::Center))
        .column(Column::auto().at_least(60.0)) // ID
        .column(Column::auto().at_least(80.0)) // RT
        .column(Column::auto().at_least(80.0)) // Mobility
        .column(Column::auto().at_least(120.0)) // Precursor m/z
        .column(Column::auto().at_least(100.0)) // Fragment count
        .header(20.0, |mut header| {
            header.col(|ui| {
                ui.strong("ID");
            });
            header.col(|ui| {
                ui.strong("RT (s)");
            });
            header.col(|ui| {
                ui.strong("Mobility");
            });
            header.col(|ui| {
                ui.strong("Precursor m/z");
            });
            header.col(|ui| {
                ui.strong("Fragments");
            });
        })
        .body(|body| {
            let row_height = 18.0;
            body.rows(row_height, filtered_eg_idxs.len(), |mut row| {
                let row_idx = row.index();
                let original_idx = filtered_eg_idxs[row_idx];
                let eg = &reference_eg_slice[original_idx];
                let is_selected = Some(original_idx) == *selected_index;

                row.col(|ui| {
                    if ui
                        .selectable_label(is_selected, format!("{}", eg.id()))
                        .clicked()
                    {
                        *selected_index = Some(original_idx);
                    }
                });

                row.col(|ui| {
                    let text = format!("{:.2}", eg.rt_seconds());
                    ui.label(text);
                });

                row.col(|ui| {
                    let text = format!("{:.4}", eg.mobility_ook0());
                    ui.label(text);
                });

                row.col(|ui| {
                    let lims = eg.get_precursor_mz_limits();
                    let display_text = format!("{:.4} - {:.4}", lims.0, lims.1);

                    ui.label(display_text);
                });

                row.col(|ui| {
                    let text = format!("{}", eg.fragment_count());
                    ui.label(text);
                });
            });
        });
}
