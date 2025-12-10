use eframe::egui;

use crate::app::AppCommand;
use crate::file_loader::ElutionGroupData;
use crate::ui::components::precursor_table;
use crate::ui::{
    Panel,
    PanelContext,
};
/// Panel for displaying and filtering the precursor table
pub struct TablePanel;

impl TablePanel {
    pub fn new() -> Self {
        Self
    }

    fn render_search_ui(&self, ui: &mut egui::Ui, ctx: &mut PanelContext) {
        ui.horizontal(|ui| {
            ui.label("Search:");
            let response = ui.text_edit_singleline(&mut ctx.ui.search_input);
            response.request_focus();
            ui.label("(Enter to apply, Esc to cancel)");
        });
        ui.separator();
    }

    fn render_filter_ui(&self, ui: &mut egui::Ui, ctx: &mut PanelContext) {
        ui.horizontal(|ui| {
            ui.label("Filter by ID:");
            ui.text_edit_singleline(&mut ctx.ui.table_filter);
            if ui.button("Clear").clicked() {
                ctx.ui.table_filter.clear();
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
        commands: &mut crate::ui::CommandSink,
    ) {
        let old_selection = *selected_index;

        egui::ScrollArea::vertical()
            .auto_shrink([false; 2])
            .show(ui, |ui| {
                macro_rules! render_table {
                    ($egs:expr) => {
                        precursor_table::render_precursor_table_filtered(
                            ui,
                            filtered_indices,
                            $egs,
                            selected_index,
                        )
                    };
                }
                crate::with_elution_collection!(elution_groups, render_table);
            });

        if old_selection != *selected_index
            && let Some(new_idx) = *selected_index
        {
            commands.push(AppCommand::SelectElutionGroup(new_idx));
        }
    }
}

impl Panel for TablePanel {
    fn render(&mut self, ui: &mut egui::Ui, ctx: &mut PanelContext) {
        ui.heading("Precursor Table");
        ui.separator();

        // Check if we have data first
        if ctx.data.elution_groups.is_none() {
            ui.label("Load elution groups to see the table");
            return;
        }

        // Render search/filter UI (needs mutable ctx, doesn't need elution_groups)
        if ctx.ui.search_mode {
            self.render_search_ui(ui, ctx);
        } else {
            self.render_filter_ui(ui, ctx);
        }

        // Now borrow elution_groups for the rest
        let elution_groups = ctx.data.elution_groups.as_ref().unwrap();
        let filtered_indices = elution_groups.matching_indices_for_id_filter(&ctx.ui.table_filter);

        ui.label(format!(
            "Showing {} of {} precursors",
            filtered_indices.len(),
            elution_groups.len()
        ));

        self.render_table(
            ui,
            &filtered_indices,
            elution_groups,
            &mut ctx.ui.selected_index,
            &mut ctx.commands,
        );
    }

    fn title(&self) -> &str {
        "Table"
    }
}

impl Default for TablePanel {
    fn default() -> Self {
        Self::new()
    }
}
