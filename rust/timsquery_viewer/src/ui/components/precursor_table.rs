use eframe::egui;
use timsquery::KeyLike;
use timsquery::models::elution_group::TimsElutionGroup;

pub fn render_precursor_table_filtered<T: KeyLike>(
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
