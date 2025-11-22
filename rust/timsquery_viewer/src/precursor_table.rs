use eframe::egui;
use timsquery::models::elution_group::ElutionGroup;

/// Renders the precursor table with selectable rows (filtered version)
pub fn render_precursor_table_filtered(
    ui: &mut egui::Ui,
    filtered_groups: &[(usize, &ElutionGroup<usize>)],
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
        .body(|mut body| {
            for &(idx, eg) in filtered_groups {
                let row_height = 18.0;
                let is_selected = Some(idx) == *selected_index;

                body.row(row_height, |mut row| {
                    // ID column
                    row.col(|ui| {
                        if ui
                            .selectable_label(is_selected, format!("{}", eg.id))
                            .clicked()
                        {
                            *selected_index = Some(idx);
                        }
                    });

                    // RT column
                    row.col(|ui| {
                        if ui
                            .selectable_label(is_selected, format!("{:.2}", eg.rt_seconds))
                            .clicked()
                        {
                            *selected_index = Some(idx);
                        }
                    });

                    // Mobility column
                    row.col(|ui| {
                        if ui
                            .selectable_label(is_selected, format!("{:.4}", eg.mobility))
                            .clicked()
                        {
                            *selected_index = Some(idx);
                        }
                    });

                    // Precursor m/z column (may have multiple)
                    row.col(|ui| {
                        let precursor_mzs: Vec<String> = eg
                            .precursors
                            .iter()
                            .map(|(_, mz)| format!("{:.4}", mz))
                            .collect();
                        let display_text = if precursor_mzs.len() > 2 {
                            format!(
                                "{}, {} (+{})",
                                precursor_mzs[0],
                                precursor_mzs[1],
                                precursor_mzs.len() - 2
                            )
                        } else {
                            precursor_mzs.join(", ")
                        };

                        if ui.selectable_label(is_selected, display_text).clicked() {
                            *selected_index = Some(idx);
                        }
                    });

                    // Fragment count column
                    row.col(|ui| {
                        if ui
                            .selectable_label(is_selected, format!("{}", eg.fragments.len()))
                            .clicked()
                        {
                            *selected_index = Some(idx);
                        }
                    });
                });
            }
        });
}
