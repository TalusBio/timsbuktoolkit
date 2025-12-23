use eframe::egui;

use crate::file_loader::ElutionGroupData;

enum Modes {
    Insert,
    Normal,
}

/// Panel for displaying and filtering the precursor table
pub struct TablePanel {
    filtered_indices: Option<Vec<usize>>,
    last_search: Option<String>,
    last_displayed_mode: Modes,
    last_selected_index: Option<usize>,
}

impl TablePanel {
    pub fn new() -> Self {
        Self {
            filtered_indices: None,
            last_search: None,
            last_displayed_mode: Modes::Normal,
            last_selected_index: None,
        }
    }

    pub fn filtered_indices(&self) -> &[usize] {
        match &self.filtered_indices {
            Some(indices) => indices,
            None => &[],
        }
    }

    fn render_search_ui(&mut self, ui: &mut egui::Ui, search_line: &mut String, search_mode: bool) {
        ui.horizontal(|ui| {
            ui.label("Search:");
            let response = ui.text_edit_singleline(search_line);
            // TODO make some color change to indicate mode

            let last_mode = &mut self.last_displayed_mode;
            if search_mode {
                response.request_focus();
                *last_mode = Modes::Insert;
            } else if let Modes::Insert = last_mode {
                // Just exited insert mode
                *last_mode = Modes::Normal;
            }
        });
        ui.separator();
    }

    fn render_keybinding_ui(&self, ui: &mut egui::Ui) {
        // TODO: Add ...
        ui.separator();
    }

    fn render_table(
        &self,
        ui: &mut egui::Ui,
        filtered_indices: &[usize],
        elution_groups: &ElutionGroupData,
        selected_index: &mut Option<usize>,
        scroll_to_selection: bool,
    ) {
        egui::ScrollArea::vertical()
            .auto_shrink([false; 2])
            .show(ui, |ui| {
                elution_groups.render_table(
                    ui,
                    filtered_indices,
                    selected_index,
                    scroll_to_selection,
                );
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

        self.render_search_ui(ui, search_line, search_mode);
        self.render_keybinding_ui(ui);

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

        // We scroll to selection only if the selection changed
        let scroll_to_selection = match (self.last_selected_index, *selected_index) {
            (Some(last), Some(current)) => last != current,
            (None, Some(_)) => true,
            _ => false,
        };
        if scroll_to_selection {
            self.last_selected_index = *selected_index;
        }

        self.render_table(
            ui,
            filtered_indices,
            elution_groups,
            selected_index,
            scroll_to_selection,
        );
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
