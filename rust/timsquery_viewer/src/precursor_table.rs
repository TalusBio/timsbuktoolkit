use eframe::egui;
use timsquery::models::elution_group::ElutionGroup;

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
        .body(| body| {
            // for &(idx, eg) in filtered_groups {
            let row_height = 18.0;
            body.rows(
                row_height,
                filtered_groups.len(),
                |mut row| {
                    let idx = row.index();
                    let eg = filtered_groups[idx].1;
                    let is_selected = Some(idx) == *selected_index;
                    // ID column - clickable
                    row.col(|ui| {
                        if ui
                            .selectable_label(is_selected, format!("{}", eg.id))
                            .clicked()
                        {
                            *selected_index = Some(idx);
                        }
                    });

                    // RT column - plain text with optional highlighting
                    row.col(|ui| {
                        let text = format!("{:.2}", eg.rt_seconds);
                        ui.label(text);
                    });

                    // Mobility column - plain text with optional highlighting
                    row.col(|ui| {
                        let text = format!("{:.4}", eg.mobility);
                        ui.label(text);
                    });

                    // Precursor m/z column - plain text with optional highlighting
                    row.col(|ui| {
                        let display_text = format!( "{:.4} ({})",
                                eg.precursors.first().unwrap().1,
                                eg.precursors.len(),
                            );

                        ui.label(display_text);
                    });

                    // Fragment count column - plain text with optional highlighting
                    row.col(|ui| {
                        let text = format!("{}", eg.fragments.len());
                        ui.label(text);
                    });
                });
        });
}

