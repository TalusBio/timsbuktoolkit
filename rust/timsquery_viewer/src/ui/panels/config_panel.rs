use eframe::egui;

use crate::chromatogram_processor::SmoothingMethod;
use crate::plot_renderer::AutoZoomMode;
use crate::ui::tolerance_editor;
use timsquery::Tolerance;

// UI spacing constants
const SECTION_MARGIN: i8 = 10;
const SECTION_SPACING: f32 = 12.0;
const INTERNAL_SPACING: f32 = 8.0;

/// Panel for configuration settings
pub struct ConfigPanel;

impl ConfigPanel {
    pub fn new() -> Self {
        Self
    }

    /// Renders tolerance editor section
    fn render_tolerance_editor(&self, ui: &mut egui::Ui, tolerance: &mut Tolerance) {
        egui::Frame::group(ui.style())
            .inner_margin(egui::Margin::same(SECTION_MARGIN))
            .show(ui, |ui| {
                ui.heading("Tolerance Settings");
                ui.add_space(INTERNAL_SPACING);

                tolerance_editor::render_tolerance_editor(ui, tolerance);
            });
    }

    /// Renders smoothing configuration section
    fn render_smoothing(&self, ui: &mut egui::Ui, smoothing_method: &mut SmoothingMethod) {
        egui::Frame::group(ui.style())
            .inner_margin(egui::Margin::same(SECTION_MARGIN))
            .show(ui, |ui| {
                ui.heading("Smoothing");
                ui.add_space(INTERNAL_SPACING);

                ui.label("Method:");
                let current_method = match smoothing_method {
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
                            *smoothing_method = SmoothingMethod::None;
                        }
                        if ui
                            .selectable_value(&mut selected, 1, "Savitzky-Golay")
                            .clicked()
                        {
                            *smoothing_method = SmoothingMethod::default_savitzky_golay();
                        }
                        if ui.selectable_value(&mut selected, 2, "Gaussian").clicked() {
                            *smoothing_method = SmoothingMethod::default_gaussian();
                        }
                    });

                ui.add_space(INTERNAL_SPACING);

                match smoothing_method {
                    SmoothingMethod::None => {
                        ui.label(egui::RichText::new("No smoothing applied").weak().italics());
                    }
                    SmoothingMethod::SavitzkyGolay { window, polynomial } => {
                        ui.label(egui::RichText::new("Parameters:").strong().small());
                        ui.horizontal(|ui| {
                            ui.label("Window size:");
                            let mut window_val = *window as i32;
                            ui.add(egui::Slider::new(&mut window_val, 3..=21).step_by(2.0));
                            *window = window_val as usize;
                        });
                        ui.horizontal(|ui| {
                            ui.label("Polynomial:");
                            let mut poly_val = *polynomial as i32;
                            ui.add(egui::Slider::new(&mut poly_val, 0..=5));
                            *polynomial = poly_val as usize;
                        });
                    }
                    SmoothingMethod::Gaussian { sigma } => {
                        ui.label(egui::RichText::new("Parameters:").strong().small());
                        ui.horizontal(|ui| {
                            ui.label("Sigma:");
                            ui.add(egui::Slider::new(sigma, 0.5..=10.0));
                        });
                    }
                }
            });
    }

    pub fn render(
        &mut self,
        ui: &mut egui::Ui,
        tolerance: &mut Tolerance,
        smoothing_method: &mut SmoothingMethod,
        auto_zoom_mode: &mut AutoZoomMode,
    ) {
        // Configuration Section Header
        ui.label(egui::RichText::new("CONFIGURATION").strong().size(13.0));
        ui.add_space(INTERNAL_SPACING);

        self.render_tolerance_editor(ui, tolerance);
        ui.add_space(SECTION_SPACING);

        self.render_smoothing(ui, smoothing_method);
        ui.add_space(SECTION_SPACING);

        // Auto Zoom section
        egui::Frame::group(ui.style())
            .inner_margin(egui::Margin::same(SECTION_MARGIN))
            .show(ui, |ui| {
                ui.heading("Auto Zoom");
                ui.add_space(INTERNAL_SPACING);

                ui.horizontal(|ui| {
                    ui.label("Mode:");
                    egui::ComboBox::from_id_salt("auto_zoom_mode")
                        .selected_text(match auto_zoom_mode {
                            AutoZoomMode::Disabled => "None",
                            AutoZoomMode::PeakApex => "Peak Score Apex",
                            AutoZoomMode::QueryRange => "Query Range",
                        })
                        .show_ui(ui, |ui| {
                            if ui
                                .selectable_value(auto_zoom_mode, AutoZoomMode::Disabled, "None")
                                .clicked()
                            {
                                *auto_zoom_mode = AutoZoomMode::Disabled;
                            }
                            if ui
                                .selectable_value(
                                    auto_zoom_mode,
                                    AutoZoomMode::PeakApex,
                                    "Peak Score Apex",
                                )
                                .clicked()
                            {
                                *auto_zoom_mode = AutoZoomMode::PeakApex;
                            }
                            if ui
                                .selectable_value(
                                    auto_zoom_mode,
                                    AutoZoomMode::QueryRange,
                                    "Query Range",
                                )
                                .clicked()
                            {
                                *auto_zoom_mode = AutoZoomMode::QueryRange;
                            }
                        });
                });
            });
    }

    pub fn title(&self) -> &str {
        "Settings"
    }
}

impl Default for ConfigPanel {
    fn default() -> Self {
        Self::new()
    }
}
