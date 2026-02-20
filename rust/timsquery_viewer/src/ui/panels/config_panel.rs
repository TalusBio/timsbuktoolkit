use eframe::egui;
use std::sync::Arc;
use std::time::Instant;

use crate::chromatogram_processor::SmoothingMethod;
use crate::plot_renderer::AutoZoomMode;
use crate::ui::tolerance_editor;
use timsquery::Tolerance;

// UI spacing constants
const SECTION_MARGIN: i8 = 10;
const SECTION_SPACING: f32 = 12.0;
const INTERNAL_SPACING: f32 = 8.0;

/// Screenshot capture lifecycle
pub enum ScreenshotState {
    /// Nothing happening
    Idle,
    /// Timer running, show remaining seconds overlay
    Countdown { deadline: Instant },
    /// Deadline reached, scale applied, screenshot command sent this frame
    Capturing,
    /// Screenshot received, saving to file
    Saving(Arc<egui::ColorImage>),
}

impl Default for ScreenshotState {
    fn default() -> Self {
        Self::Idle
    }
}

/// Actions the export UI can request
pub enum ScreenshotAction {
    None,
    Start,
    Cancel,
}

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

    /// Renders the Export/Screenshot section at the bottom of the config panel
    pub fn render_export_section(
        ui: &mut egui::Ui,
        screenshot_delay_secs: &mut f32,
        screenshot_state: &ScreenshotState,
    ) -> ScreenshotAction {
        let mut action = ScreenshotAction::None;

        ui.label(egui::RichText::new("EXPORT").strong().size(13.0));
        ui.add_space(INTERNAL_SPACING);

        egui::Frame::group(ui.style())
            .inner_margin(egui::Margin::same(SECTION_MARGIN))
            .show(ui, |ui| {
                ui.heading("Screenshot");
                ui.add_space(INTERNAL_SPACING);

                let is_active = !matches!(screenshot_state, ScreenshotState::Idle);

                // Delay selector
                ui.horizontal(|ui| {
                    ui.label("Delay:");
                    ui.add_enabled_ui(!is_active, |ui| {
                        egui::ComboBox::from_id_salt("screenshot_delay")
                            .selected_text(if *screenshot_delay_secs == 0.0 {
                                "None".to_string()
                            } else {
                                format!("{}s", *screenshot_delay_secs as u32)
                            })
                            .show_ui(ui, |ui| {
                                ui.selectable_value(screenshot_delay_secs, 0.0, "None");
                                ui.selectable_value(screenshot_delay_secs, 3.0, "3s");
                                ui.selectable_value(screenshot_delay_secs, 5.0, "5s");
                                ui.selectable_value(screenshot_delay_secs, 10.0, "10s");
                            });
                    });
                });

                ui.add_space(INTERNAL_SPACING);

                // Action buttons
                match screenshot_state {
                    ScreenshotState::Countdown { deadline } => {
                        let remaining = deadline
                            .saturating_duration_since(Instant::now())
                            .as_secs_f32()
                            .ceil() as u32;
                        ui.horizontal(|ui| {
                            ui.add_enabled(false, egui::Button::new(
                                format!("Capturing in {}s...", remaining),
                            ));
                            if ui.button("Cancel").clicked() {
                                action = ScreenshotAction::Cancel;
                            }
                        });
                    }
                    ScreenshotState::Capturing => {
                        ui.add_enabled(false, egui::Button::new("Capturing..."));
                    }
                    ScreenshotState::Saving(_) => {
                        ui.add_enabled(false, egui::Button::new("Saving..."));
                    }
                    ScreenshotState::Idle => {
                        if ui.button("Save Screenshot").clicked() {
                            action = ScreenshotAction::Start;
                        }
                    }
                }
            });

        action
    }
}

impl Default for ConfigPanel {
    fn default() -> Self {
        Self::new()
    }
}
