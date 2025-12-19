use eframe::egui;

use crate::app::AppCommand;
use crate::chromatogram_processor::SmoothingMethod;
use crate::ui::{
    Panel,
    PanelContext,
    tolerance_editor,
};

/// Panel for data loading and settings on the left side
pub struct ConfigPanel;

impl ConfigPanel {
    pub fn new() -> Self {
        Self
    }

    /// Renders tolerance editor section
    fn render_tolerance_editor(&self, ui: &mut egui::Ui, ctx: &mut PanelContext) {
        ui.heading("Tolerance Settings");

        let changed = tolerance_editor::render_tolerance_editor(ui, &mut ctx.data.tolerance);
        if changed {
            ctx.commands.push(AppCommand::UpdateTolerance);
        }
    }

    /// Renders smoothing configuration section
    fn render_smoothing(&self, ui: &mut egui::Ui, ctx: &mut PanelContext) {
        ui.heading("Smoothing");

        let mut changed = false;

        ui.label("Method:");
        let current_method = match ctx.data.smoothing {
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
                    ctx.data.smoothing = SmoothingMethod::None;
                    changed = true;
                }
                if ui
                    .selectable_value(&mut selected, 1, "Savitzky-Golay")
                    .clicked()
                {
                    ctx.data.smoothing = SmoothingMethod::default_savitzky_golay();
                    changed = true;
                }
                if ui.selectable_value(&mut selected, 2, "Gaussian").clicked() {
                    ctx.data.smoothing = SmoothingMethod::default_gaussian();
                    changed = true;
                }
            });

        ui.add_space(10.0);

        match &mut ctx.data.smoothing {
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
            ctx.commands.push(AppCommand::UpdateSmoothing);
        }
    }
}

impl Panel for ConfigPanel {
    fn render(&mut self, ui: &mut egui::Ui, ctx: &mut PanelContext) {
        self.render_tolerance_editor(ui, ctx);
        ui.add_space(20.0);
        ui.separator();

        self.render_smoothing(ui, ctx);
        ui.add_space(20.0);
        ui.separator();

        // Plot visibility controls
        ui.heading("Plot Settings");
        ui.checkbox(&mut ctx.ui.show_ms2_spectrum, "Show MS2 Spectrum");
    }

    fn title(&self) -> &str {
        "Settings"
    }
}

impl Default for ConfigPanel {
    fn default() -> Self {
        Self::new()
    }
}
