use eframe::egui;

use crate::file_loader::ElutionGroupData;
/// Panel for displaying and filtering the precursor table
pub struct TablePanel;

pub enum TablePanelMessage {
    SelectElutionGroup(usize),
    None,
}

impl TablePanel {
    pub fn new() -> Self {
        Self
    }

    fn render_search_ui(&self, ui: &mut egui::Ui, search_line: &mut String) {
        ui.horizontal(|ui| {
            ui.label("Search:");
            let response = ui.text_edit_singleline(search_line);
            response.request_focus();
            ui.label("(Enter to apply, Esc to cancel)");
        });
        ui.separator();
    }

    fn render_filter_ui(&self, ui: &mut egui::Ui, table_filter: &mut String) {
        ui.horizontal(|ui| {
            ui.label("Filter by ID:");
            ui.text_edit_singleline(table_filter);
            if ui.button("Clear").clicked() {
                table_filter.clear();
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
        selected_index: &mut Option<usize>,
    ) -> TablePanelMessage {
        let old_selection = *selected_index;

        egui::ScrollArea::vertical()
            .auto_shrink([false; 2])
            .show(ui, |ui| {
                elution_groups.render_table(ui, filtered_indices, selected_index);
            });

        if old_selection != *selected_index
            && let Some(new_idx) = *selected_index
        {
            TablePanelMessage::SelectElutionGroup(new_idx)
        } else {
            TablePanelMessage::None
        }
    }

    pub fn render(
        &mut self,
        ui: &mut egui::Ui,
        // ctx: &mut PanelContext,
        elution_groups: &Option<ElutionGroupData>,
        search_mode: bool,
        search_line: &mut String,
        selected_index: &mut Option<usize>,
    ) {
        ui.heading("Precursor Table");
        ui.separator();

        // Check if we have data first
        if elution_groups.is_none() {
            ui.label("Load elution groups to see the table");
            return;
        }

        // Render search/filter UI (needs mutable ctx, doesn't need elution_groups)
        if search_mode {
            self.render_search_ui(ui, search_line);
        } else {
            self.render_filter_ui(ui, search_line);
        }

        // Now borrow elution_groups for the rest
        let elution_groups = elution_groups.as_ref().unwrap();
        let filtered_indices = elution_groups.matching_indices_for_id_filter(search_line);

        ui.label(format!(
            "Showing {} of {} precursors",
            filtered_indices.len(),
            elution_groups.len()
        ));

        self.render_table(ui, &filtered_indices, elution_groups, selected_index);
    }

    pub fn title(&self) -> &str {
        "Table"
    }
}

impl Default for TablePanel {
    fn default() -> Self {
        Self::new()
    }
}
