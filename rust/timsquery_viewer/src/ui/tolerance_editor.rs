use eframe::egui;
use timsquery::models::tolerance::{
    MobilityTolerance,
    MzTolerance,
    QuadTolerance,
    RtTolerance,
    Tolerance,
};

/// Renders the tolerance editor UI
pub fn render_tolerance_editor(ui: &mut egui::Ui, tolerance: &mut Tolerance) {
    ui.collapsing("M/Z Tolerance", |ui| {
        render_mz_tolerance(ui, &mut tolerance.ms);
    });

    ui.collapsing("RT Tolerance", |ui| {
        render_rt_tolerance(ui, &mut tolerance.rt);
    });

    ui.collapsing("Mobility Tolerance", |ui| {
        render_mobility_tolerance(ui, &mut tolerance.mobility);
    });

    ui.collapsing("Quadrupole Tolerance", |ui| {
        render_quad_tolerance(ui, &mut tolerance.quad);
    });
}

fn render_mz_tolerance(ui: &mut egui::Ui, tol: &mut MzTolerance) {
    let mut changed = false;

    let mut is_ppm = matches!(tol, MzTolerance::Ppm(_));
    let old_is_ppm = is_ppm;

    ui.horizontal(|ui| {
        ui.label("Type:");
        ui.radio_value(&mut is_ppm, false, "Absolute (Da)");
        ui.radio_value(&mut is_ppm, true, "PPM");
    });

    if is_ppm != old_is_ppm {
        *tol = if is_ppm {
            MzTolerance::Ppm((15.0, 15.0))
        } else {
            MzTolerance::Absolute((0.01, 0.01))
        };
        changed = true;
    }

    match tol {
        MzTolerance::Ppm((lower, upper)) => {
            ui.horizontal(|ui| {
                ui.label("Lower (ppm):");
                changed |= ui.add(egui::DragValue::new(lower).speed(0.1)).changed();
            });
            ui.horizontal(|ui| {
                ui.label("Upper (ppm):");
                changed |= ui.add(egui::DragValue::new(upper).speed(0.1)).changed();
            });
        }
        MzTolerance::Absolute((lower, upper)) => {
            ui.horizontal(|ui| {
                ui.label("Lower (Da):");
                changed |= ui.add(egui::DragValue::new(lower).speed(0.001)).changed();
            });
            ui.horizontal(|ui| {
                ui.label("Upper (Da):");
                changed |= ui.add(egui::DragValue::new(upper).speed(0.001)).changed();
            });
        }
    }
}

fn render_rt_tolerance(ui: &mut egui::Ui, tol: &mut RtTolerance) {
    let mut selected = match tol {
        RtTolerance::Minutes(_) => 0,
        RtTolerance::Pct(_) => 1,
        RtTolerance::Unrestricted => 2,
    };
    let old_selected = selected;

    ui.horizontal(|ui| {
        ui.label("Type:");
        ui.radio_value(&mut selected, 0, "Minutes");
        ui.radio_value(&mut selected, 1, "Percent");
        ui.radio_value(&mut selected, 2, "Unrestricted");
    });

    if selected != old_selected {
        *tol = match selected {
            0 => RtTolerance::Minutes((0.5, 0.5)),
            1 => RtTolerance::Pct((10.0, 10.0)),
            2 => RtTolerance::Unrestricted,
            _ => RtTolerance::Unrestricted,
        };
    }

    match tol {
        RtTolerance::Minutes((lower, upper)) => {
            ui.horizontal(|ui| {
                ui.label("Lower (min):");
                ui.add(egui::DragValue::new(lower).speed(0.1));
            });
            ui.horizontal(|ui| {
                ui.label("Upper (min):");
                ui.add(egui::DragValue::new(upper).speed(0.1));
            });
        }
        RtTolerance::Pct((lower, upper)) => {
            ui.horizontal(|ui| {
                ui.label("Lower (%):");
                ui.add(egui::DragValue::new(lower).speed(0.5));
            });
            ui.horizontal(|ui| {
                ui.label("Upper (%):");
                ui.add(egui::DragValue::new(upper).speed(0.5));
            });
        }
        RtTolerance::Unrestricted => {
            ui.label("No RT restriction");
        }
    }
}

fn render_mobility_tolerance(ui: &mut egui::Ui, tol: &mut MobilityTolerance) -> bool {
    let mut changed = false;

    let mut selected = match tol {
        MobilityTolerance::Absolute(_) => 0,
        MobilityTolerance::Pct(_) => 1,
        MobilityTolerance::Unrestricted => 2,
    };
    let old_selected = selected;

    ui.horizontal(|ui| {
        ui.label("Type:");
        ui.radio_value(&mut selected, 0, "Absolute");
        ui.radio_value(&mut selected, 1, "Percent");
        ui.radio_value(&mut selected, 2, "Unrestricted");
    });

    if selected != old_selected {
        *tol = match selected {
            0 => MobilityTolerance::Absolute((0.01, 0.01)),
            1 => MobilityTolerance::Pct((3.0, 3.0)),
            2 => MobilityTolerance::Unrestricted,
            _ => MobilityTolerance::Unrestricted,
        };
        changed = true;
    }

    match tol {
        MobilityTolerance::Absolute((lower, upper)) => {
            ui.horizontal(|ui| {
                ui.label("Lower:");
                changed |= ui.add(egui::DragValue::new(lower).speed(0.001)).changed();
            });
            ui.horizontal(|ui| {
                ui.label("Upper:");
                changed |= ui.add(egui::DragValue::new(upper).speed(0.001)).changed();
            });
        }
        MobilityTolerance::Pct((lower, upper)) => {
            ui.horizontal(|ui| {
                ui.label("Lower (%):");
                changed |= ui.add(egui::DragValue::new(lower).speed(0.5)).changed();
            });
            ui.horizontal(|ui| {
                ui.label("Upper (%):");
                changed |= ui.add(egui::DragValue::new(upper).speed(0.5)).changed();
            });
        }
        MobilityTolerance::Unrestricted => {
            ui.label("No mobility restriction");
        }
    }

    changed
}

fn render_quad_tolerance(ui: &mut egui::Ui, tol: &mut QuadTolerance) -> bool {
    let mut changed = false;

    match tol {
        QuadTolerance::Absolute((lower, upper)) => {
            ui.horizontal(|ui| {
                ui.label("Lower (abs):");
                changed |= ui.add(egui::DragValue::new(lower).speed(0.01)).changed();
            });
            ui.horizontal(|ui| {
                ui.label("Upper (abs):");
                changed |= ui.add(egui::DragValue::new(upper).speed(0.01)).changed();
            });
        }
    }

    changed
}
