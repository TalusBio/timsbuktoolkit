use eframe::egui;

use crate::app::{
    AppCommand,
    DataState,
};
use crate::chromatogram_processor::SmoothingMethod;
use crate::tolerance_editor;

/// Panel for data loading and settings on the left side
pub struct LeftPanel;

impl LeftPanel {
    pub fn new() -> Self {
        Self
    }

    /// Render tolerance editor section
    /// Returns commands if settings changed
    pub fn render_tolerance_editor(
        &self,
        ui: &mut egui::Ui,
        data: &mut DataState,
    ) -> Vec<AppCommand> {
        let mut commands = Vec::new();

        ui.heading("Tolerance Settings");

        let changed = tolerance_editor::render_tolerance_editor(ui, &mut data.tolerance);
        if changed {
            commands.push(AppCommand::UpdateTolerance);
        }

        commands
    }

    /// Render smoothing section
    /// Returns commands if settings changed
    pub fn render_smoothing(
        &self,
        ui: &mut egui::Ui,
        data: &mut DataState,
    ) -> Vec<AppCommand> {
        let mut commands = Vec::new();
        ui.heading("Smoothing");

        let mut changed = false;

        ui.label("Method:");
        let current_method = match data.smoothing {
            SmoothingMethod::None => 0,
            SmoothingMethod::SavitzkyGolay { .. } => 1,
            SmoothingMethod::Gaussian { .. } => 2,
        };

        let mut selected = current_method;
        egui::ComboBox::from_id_salt("smoothing_method")
            .selected_text(match current_method {
                0 => "None",
                1 => "Savitzky-Golay",
                2 => "Gaussian",
                _ => "Unknown",
            })
            .show_ui(ui, |ui| {
                if ui.selectable_value(&mut selected, 0, "None").clicked() {
                    data.smoothing = SmoothingMethod::None;
                    changed = true;
                }
                if ui
                    .selectable_value(&mut selected, 1, "Savitzky-Golay")
                    .clicked()
                {
                    data.smoothing = SmoothingMethod::SavitzkyGolay {
                        window: 5,
                        polynomial: 2,
                    };
                    changed = true;
                }
                if ui.selectable_value(&mut selected, 2, "Gaussian").clicked() {
                    data.smoothing = SmoothingMethod::Gaussian { sigma: 2.0 };
                    changed = true;
                }
            });

        ui.add_space(10.0);

        match &mut data.smoothing {
            SmoothingMethod::None => {
                ui.label("No smoothing applied");
            }
            SmoothingMethod::SavitzkyGolay { window, polynomial } => {
                ui.label("Parameters:");
                ui.horizontal(|ui| {
                    ui.label("Window size:");
                    let mut window_val = *window as i32;
                    if ui
                        .add(egui::Slider::new(&mut window_val, 3..=21).step_by(2.0))
                        .changed()
                    {
                        *window = window_val as usize;
                        changed = true;
                    }
                });
                ui.horizontal(|ui| {
                    ui.label("Polynomial:");
                    let mut poly_val = *polynomial as i32;
                    if ui.add(egui::Slider::new(&mut poly_val, 0..=5)).changed() {
                        *polynomial = poly_val as usize;
                        changed = true;
                    }
                });
            }
            SmoothingMethod::Gaussian { sigma } => {
                ui.label("Parameters:");
                ui.horizontal(|ui| {
                    ui.label("Sigma:");
                    if ui.add(egui::Slider::new(sigma, 0.5..=10.0)).changed() {
                        changed = true;
                    }
                });
            }
        }

        if changed {
            commands.push(AppCommand::UpdateSmoothing);
        }

        commands
    }
}

impl Default for LeftPanel {
    fn default() -> Self {
        Self::new()
    }
}
