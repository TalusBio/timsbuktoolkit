use eframe::egui;

use crate::file_loader::ElutionGroupData;
/// Panel for displaying and filtering the precursor table
pub struct TablePanel {
    filtered_indices: Option<Vec<usize>>,
    last_search: Option<String>,
}

impl TablePanel {
    pub fn new() -> Self {
        Self {
            filtered_indices: None,
            last_search: None,
        }
    }

    pub fn filtered_indices(&self) -> &[usize] {
        match &self.filtered_indices {
            Some(indices) => indices,
            None => &[],
        }
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
    ) {
        egui::ScrollArea::vertical()
            .auto_shrink([false; 2])
            .show(ui, |ui| {
                elution_groups.render_table(ui, filtered_indices, selected_index);
            });
    }

    pub fn render(
        &mut self,
        ui: &mut egui::Ui,
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

        if search_mode {
            self.render_search_ui(ui, search_line);
        } else {
            self.render_filter_ui(ui, search_line);
        }

        let elution_groups = elution_groups.as_ref().unwrap();
        // Invalidate cache if search line changed
        // TODO: I think I can optimize search using the fact that when typing
        // letters are added, so I can filter from previous results.
        if self.last_search.as_ref() != Some(search_line) {
            // Update filtered indices
            elution_groups.matching_indices_for_id_filter(
                search_line,
                self.filtered_indices.get_or_insert_with(Vec::new),
            );
            Some(());
            self.last_search = Some(search_line.clone());
        }
        let filtered_indices = self.filtered_indices.as_ref().unwrap();

        ui.label(format!(
            "Showing {} of {} precursors",
            filtered_indices.len(),
            elution_groups.len()
        ));

        self.render_table(ui, filtered_indices, elution_groups, selected_index);
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
