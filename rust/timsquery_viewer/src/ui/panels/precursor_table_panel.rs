use eframe::egui;

use crate::app::{
    AppCommand,
    UiState,
};
use crate::file_loader::ElutionGroupData;
use crate::ui::components::precursor_table;

/// Panel for displaying and filtering the precursor table
pub struct TablePanel;

impl TablePanel {
    pub fn new() -> Self {
        Self
    }

    /// Renders the precursor table panel with filtering and selection.
    ///
    /// Reads: `elution_groups` (all elution group data)
    ///        `ui_state.table_filter` (current filter string)
    ///        `ui_state.search_mode` (whether search is active)
    ///        `ui_state.selected_index` (currently selected row)
    /// Writes: `ui_state.selected_index` (if user clicks a row)
    /// Returns: Command to select elution group if selection changed
    pub fn render(
        &self,
        ui: &mut egui::Ui,
        elution_groups: &ElutionGroupData,
        ui_state: &mut UiState,
    ) -> Option<AppCommand> {
        ui.heading("Precursor Table");
        ui.separator();

        if ui_state.search_mode {
            self.render_search_ui(ui, ui_state);
        } else {
            self.render_filter_ui(ui, ui_state);
        }

        let filtered_indices =
            elution_groups.matching_indices_for_id_filter(&ui_state.table_filter);

        ui.label(format!(
            "Showing {} of {} precursors",
            filtered_indices.len(),
            elution_groups.len()
        ));

        self.render_table(ui, &filtered_indices, elution_groups, ui_state)
    }

    fn render_search_ui(&self, ui: &mut egui::Ui, ui_state: &mut UiState) {
        ui.horizontal(|ui| {
            ui.label("Search:");
            let response = ui.text_edit_singleline(&mut ui_state.search_input);
            response.request_focus();
            ui.label("(Enter to apply, Esc to cancel)");
        });
        ui.separator();
    }

    fn render_filter_ui(&self, ui: &mut egui::Ui, ui_state: &mut UiState) {
        ui.horizontal(|ui| {
            ui.label("Filter by ID:");
            ui.text_edit_singleline(&mut ui_state.table_filter);
            if ui.button("Clear").clicked() {
                ui_state.table_filter.clear();
            }
        });
        ui.add_space(5.0);
        ui.label(
            egui::RichText::new("Vim keys: j/k=navigate, /=search, g/G=first/last")
                .small()
                .italics(),
        );
        ui.separator();
    }

    fn render_table(
        &self,
        ui: &mut egui::Ui,
        filtered_indices: &[usize],
        elution_groups: &ElutionGroupData,
        ui_state: &mut UiState,
    ) -> Option<AppCommand> {
        let old_selection = ui_state.selected_index;

        egui::ScrollArea::vertical()
            .auto_shrink([false; 2])
            .show(ui, |ui| match elution_groups {
                ElutionGroupData::StringLabels(egs) => {
                    precursor_table::render_precursor_table_filtered(
                        ui,
                        filtered_indices,
                        egs,
                        &mut ui_state.selected_index,
                    );
                }
                ElutionGroupData::MzpafLabels(egs) => {
                    precursor_table::render_precursor_table_filtered(
                        ui,
                        filtered_indices,
                        egs,
                        &mut ui_state.selected_index,
                    );
                }
            });

        if old_selection != ui_state.selected_index {
            ui_state.selected_index.map(AppCommand::SelectElutionGroup)
        } else {
            None
        }
    }
}

impl Default for TablePanel {
    fn default() -> Self {
        Self::new()
    }
}
