//! eframe application: control panel + stacked visualizations.

use eframe::egui::{
    self,
    TextureHandle,
};

use crate::bench::{
    self,
    ScorePopulations,
};
use crate::plots;
use crate::scorer::{
    self,
    ScoreResult,
};
use crate::sim::{
    self,
    FragmentSpec,
    SimData,
    SimParams,
};

pub struct SimApp {
    params: SimParams,
    data: Option<SimData>,
    score: Option<ScoreResult>,
    score_err: Option<String>,
    heatmap: Option<(TextureHandle, usize)>,
    dirty: bool,
    dist_controls: DistControls,
    /// Sampled score populations, cleared whenever a knob changes (they would
    /// otherwise describe stale params).
    populations: Option<SampledPopulations>,
}

/// Sampled populations plus the seed-pair count they were taken with, so a
/// moved "seeds" slider is detectable as staleness.
struct SampledPopulations {
    n_seed_pairs: usize,
    pops: ScorePopulations,
}

/// Knobs for the score-distribution panel.
struct DistControls {
    /// While on, the populations are resampled whenever they go stale.
    enabled: bool,
    n_seed_pairs: usize,
    show_absent: bool,
}

impl Default for SimApp {
    fn default() -> Self {
        Self {
            params: SimParams::default(),
            data: None,
            score: None,
            score_err: None,
            heatmap: None,
            dirty: true,
            dist_controls: DistControls {
                enabled: false,
                n_seed_pairs: 200,
                show_absent: true,
            },
            populations: None,
        }
    }
}

impl SimApp {
    /// Rebuild the synthetic extraction, rescore it, and refresh the heatmap.
    fn recompute(&mut self, ctx: &egui::Context) {
        let data = sim::build(&self.params);

        let map = self.params.rt_mapper();
        match scorer::run(&data.extraction, &map) {
            Ok(res) => {
                self.score = Some(res);
                self.score_err = None;
            }
            Err(e) => {
                self.score = None;
                self.score_err = Some(e);
            }
        }

        self.heatmap = Some(plots::build_heatmap(ctx, &data));
        self.data = Some(data);
        self.populations = None; // sampled for the previous params -> stale
        self.dirty = false;
    }

    /// Resample the score populations if the ones on hand no longer describe
    /// the current params or seed-pair count. Must run AFTER [`Self::recompute`],
    /// which is what clears them on a param change.
    ///
    /// Skipped while a pointer button is down: a slider drag would otherwise
    /// resample every frame, and one sample is `2 * n_seed_pairs` full
    /// simulations.
    fn resample_if_stale(&mut self, ctx: &egui::Context) {
        if !self.dist_controls.enabled {
            self.populations = None;
            return;
        }
        let n_seed_pairs = self.dist_controls.n_seed_pairs;
        let fresh = self
            .populations
            .as_ref()
            .is_some_and(|p| p.n_seed_pairs == n_seed_pairs);
        if fresh || ctx.input(|i| i.pointer.any_down()) {
            return;
        }
        self.populations = Some(SampledPopulations {
            n_seed_pairs,
            pops: bench::score_populations(&self.params, n_seed_pairs),
        });
    }
}

impl eframe::App for SimApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::SidePanel::left("controls")
            .default_width(320.0)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    self.dirty |= controls(ui, &mut self.params, &mut self.dist_controls);
                    ui.separator();
                    score_readout(ui, self.score.as_ref(), self.score_err.as_deref());
                });
            });

        if self.dirty {
            self.recompute(ctx);
        }
        self.resample_if_stale(ctx);

        egui::CentralPanel::default().show(ctx, |ui| {
            // Every section is collapsible (a collapsed one is not laid out at
            // all, so hiding a plot also skips its rendering cost), and the
            // stack still exceeds most window heights when all are open.
            egui::ScrollArea::vertical().show(ui, |ui| {
                let Some(data) = &self.data else { return };
                let score = self.score.as_ref();
                // The REALIZED apex, not `params.apex_cycle`: the sim always
                // applies a sub-cycle offset, so the configured value is never
                // where the peak actually landed.
                let true_apex = data.realized_apex_cycle;

                section(ui, "Chromatograms", |ui| {
                    plots::apex_caption(ui);
                    plots::chromatograms(ui, data, score, true_apex);
                });

                section(ui, "Heatmap (transitions x cycles)", |ui| {
                    if let Some((tex, _)) = &self.heatmap {
                        let n_cols = data
                            .fragment_rows
                            .first()
                            .map(|r| r.intensities.len())
                            .unwrap_or(1);
                        plots::heatmap(ui, tex, n_cols, score, true_apex);
                    }
                });

                section(ui, "Apex-finder intermediate traces", |ui| {
                    if let Some(score) = score {
                        plots::traces(
                            ui,
                            score,
                            &TRACE_TOGGLES.with(|t| t.borrow().clone()),
                            true_apex,
                        );
                    } else if let Some(err) = &self.score_err {
                        ui.colored_label(egui::Color32::from_rgb(220, 80, 80), err);
                    }
                });

                section(
                    ui,
                    "Main-score across cycles (log10, fraction per bin)",
                    |ui| landscape_histogram(ui, score),
                );

                section(
                    ui,
                    "Main-score across seeds (log10, fraction per bin)",
                    |ui| {
                        seeds_histogram(
                            ui,
                            self.populations.as_ref().map(|p| &p.pops),
                            self.dist_controls.show_absent,
                            score.map(|s| s.pass2.primary.main_score),
                        )
                    },
                );
            });
        });
    }

    /// Dump the final simulation config to stdout on window close, so a run can
    /// be reproduced or pasted into a test.
    fn on_exit(&mut self, _gl: Option<&eframe::glow::Context>) {
        println!("--- apex_sim final config ---\n{:#?}", self.params);
    }
}

// Trace toggles kept in a thread-local so the checkbox row (rendered in the
// central panel header) and the plot share one source of truth without
// threading extra fields through every fn.
thread_local! {
    static TRACE_TOGGLES: std::cell::RefCell<plots::TraceToggles> =
        std::cell::RefCell::new(plots::TraceToggles::default());
}

/// One collapsible plot section, open by default. egui remembers the open/shut
/// state per title for the lifetime of the process.
fn section(ui: &mut egui::Ui, title: &str, body: impl FnOnce(&mut egui::Ui)) {
    ui.separator();
    egui::CollapsingHeader::new(egui::RichText::new(title).heading())
        .id_salt(title)
        .default_open(true)
        .show(ui, body);
}

/// OR a widget's `changed()` into `flag`. Non-capturing so it can be used
/// inside nested `ui.horizontal` closures without a persistent borrow.
fn chg(flag: &mut bool, r: egui::Response) {
    *flag |= r.changed();
}

/// Render every knob. Returns true if any *simulation* param changed (=>
/// recompute needed). The `dist` knobs deliberately do not set it: none of them
/// affect the simulation, and a changed seed-pair count is picked up by
/// [`SimApp::resample_if_stale`] instead.
fn controls(ui: &mut egui::Ui, p: &mut SimParams, dist: &mut DistControls) -> bool {
    let mut changed = false;

    ui.heading("Simulation");

    ui.label("Peak shape");
    chg(
        &mut changed,
        ui.add(
            egui::Slider::new(&mut p.apex_cycle, 0.0..=(p.n_cycles as f32 - 1.0))
                .text("apex cycle"),
        ),
    );
    chg(
        &mut changed,
        ui.add(egui::Slider::new(&mut p.width_sigma, 0.1..=5.0).text("width sigma")),
    );
    chg(
        &mut changed,
        ui.add(
            egui::Slider::new(&mut p.height, 10.0..=10_000.0)
                .logarithmic(true)
                .text("height"),
        ),
    );
    chg(
        &mut changed,
        ui.add(egui::Slider::new(&mut p.n_cycles, 10..=300).text("n cycles")),
    );

    ui.separator();
    ui.label("Noise");
    chg(
        &mut changed,
        ui.add(egui::Slider::new(&mut p.noise_floor, 0.0..=1.0).text("noise floor")),
    );
    if ui
        .add(egui::Slider::new(&mut p.seed, 0..=999).text("seed"))
        .changed()
    {
        changed = true;
    }

    ui.separator();
    ui.label("Precursor");
    chg(
        &mut changed,
        ui.add(egui::Slider::new(&mut p.n_isotopes, 1..=6).text("n isotopes")),
    );
    chg(
        &mut changed,
        ui.add(
            egui::Slider::new(&mut p.precursor_intensity, 0.0..=2.0).text("precursor intensity"),
        ),
    );

    ui.separator();
    ui.horizontal(|ui| {
        ui.label(format!("Real fragments ({})", p.real_fragments.len()));
        if ui.button("+").clicked() {
            let n = p.real_fragments.len();
            p.real_fragments.push(FragmentSpec {
                label: format!("?{n}"),
                theo_intensity: 0.5,
                obs_scale: 1.0,
                noise_mult: 1.0,
            });
            changed = true;
        }
        if ui.button("-").clicked() && p.real_fragments.len() > 1 {
            p.real_fragments.pop();
            changed = true;
        }
    });
    ui.label("theo = library ratio · obs = observed/theory (0 = absent)");
    for f in p.real_fragments.iter_mut() {
        ui.horizontal(|ui| {
            ui.label(&f.label);
            chg(
                &mut changed,
                ui.add(egui::Slider::new(&mut f.theo_intensity, 0.0..=1.0).text("theo")),
            );
            chg(
                &mut changed,
                ui.add(egui::Slider::new(&mut f.obs_scale, 0.0..=2.0).text("obs")),
            );
            chg(
                &mut changed,
                ui.add(egui::Slider::new(&mut f.noise_mult, 0.0..=5.0).text("noise x")),
            );
        });
    }

    ui.separator();
    let rp = &mut p.random_peaks;
    chg(
        &mut changed,
        ui.checkbox(&mut rp.enabled, "Random peak-like objects"),
    );
    if rp.enabled {
        chg(
            &mut changed,
            ui.add(egui::Slider::new(&mut rp.count, 0..=100).text("count")),
        );
        ui.horizontal(|ui| {
            chg(
                &mut changed,
                ui.add(egui::Slider::new(&mut rp.height_frac_min, 0.0..=10.0).text("h min")),
            );
            chg(
                &mut changed,
                ui.add(egui::Slider::new(&mut rp.height_frac_max, 0.0..=10.0).text("h max")),
            );
        });
        ui.label("(width inherited from peptide σ)");
        chg(
            &mut changed,
            ui.checkbox(&mut rp.hit_precursors, "also hit precursors"),
        );
    }

    ui.separator();
    ui.label("Score distribution");
    ui.add(
        egui::Slider::new(&mut dist.n_seed_pairs, 10..=2000)
            .logarithmic(true)
            .text("seeds"),
    );
    ui.checkbox(&mut dist.show_absent, "overlay pure-noise twin");
    ui.checkbox(&mut dist.enabled, "sample across seeds (slow)");

    ui.separator();
    ui.label("Traces to plot");
    TRACE_TOGGLES.with(|t| {
        let mut t = t.borrow_mut();
        ui.checkbox(&mut t.cosine, "cosine");
        ui.checkbox(&mut t.scribe, "scribe");
        ui.checkbox(&mut t.lazyscore, "lazyscore");
        ui.checkbox(&mut t.log_intensity, "log_intensity");
        ui.checkbox(&mut t.ms1_precursor, "ms1_precursor");
        ui.checkbox(&mut t.apex_profile, "apex_profile");
    });

    changed
}

/// The score landscape of the CURRENT run: `main_score` at every cycle of this
/// one extraction, marked with the cycle Pass 2 actually reported.
fn landscape_histogram(ui: &mut egui::Ui, score: Option<&ScoreResult>) {
    ui.label(
        "this run only: main_score at every cycle. Only the 11 features vary with cycle -- \
         the apex evidence is window-global.",
    );
    let Some(s) = score else { return };
    plots::score_histogram(
        ui,
        "hist_landscape",
        &[plots::HistSeries {
            name: "all cycles",
            values: &s.landscape,
            color: plots::PRESENT_COLOR,
        }],
        Some(s.pass2.primary.main_score),
    );
}

/// Present-vs-absent `main_score` populations across seeds, with the ROC-AUC
/// between them. `pops` is `None` until sampling is enabled in the controls.
fn seeds_histogram(
    ui: &mut egui::Ui,
    pops: Option<&ScorePopulations>,
    show_absent: bool,
    marker_score: Option<f32>,
) {
    let Some(pops) = pops else {
        ui.label(
            "enable \"sample across seeds\" in the control panel (resamples on every param \
             change)",
        );
        return;
    };

    let auc = bench::roc_auc(&pops.present, &pops.absent);
    ui.label(format!(
        "{} seeds · ROC-AUC vs pure noise = {auc:.4}",
        pops.present.len()
    ));

    let mut series = vec![plots::HistSeries {
        name: "peptide present",
        values: &pops.present,
        color: plots::PRESENT_COLOR,
    }];
    if show_absent {
        series.push(plots::HistSeries {
            name: "pure noise",
            values: &pops.absent,
            color: plots::ABSENT_COLOR,
        });
    }
    plots::score_histogram(ui, "hist_seeds", &series, marker_score);
}

/// Show the full `ApexBlocks` from pass 2 plus the pass-1 location.
fn score_readout(ui: &mut egui::Ui, score: Option<&ScoreResult>, err: Option<&str>) {
    ui.heading("Score");
    let Some(s) = score else {
        if let Some(e) = err {
            ui.colored_label(egui::Color32::from_rgb(220, 80, 80), e);
        } else {
            ui.label("(no score)");
        }
        return;
    };

    let p2 = &s.pass2;
    egui::Grid::new("score_grid").num_columns(2).show(ui, |ui| {
        let row = |ui: &mut egui::Ui, k: &str, v: String| {
            ui.label(k);
            ui.label(v);
            ui.end_row();
        };
        let ev = &p2.evidence;
        row(ui, "pass1 apex cycle", s.pass1.apex_cycle.to_string());
        row(ui, "pass1 score", format!("{:.3e}", s.pass1.score));
        row(ui, "pass2 joint apex", p2.joint_apex_cycle.to_string());
        row(ui, "pass2 score", format!("{:.3e}", p2.primary.main_score));
        row(ui, "rt (ms)", p2.retention_time_ms.to_string());
        // Score = apex_evidence * feature_product. Decompose so magnitude
        // blowups (evidence ~ intensity^2 via the two area-uniqueness terms)
        // are visible.
        row(ui, "= apex_evidence", format!("{:.3e}", ev.apex_evidence));
        row(
            ui,
            "  x feat_product",
            format!(
                "{:.3e}",
                safe_ratio(p2.primary.main_score, ev.apex_evidence)
            ),
        );
        row(ui, "  cosine_au", format!("{:.3e}", ev.cosine_au));
        row(ui, "  cosine_cg", format!("{:.3}", ev.cosine_cg));
        row(ui, "  scribe_au", format!("{:.3e}", ev.scribe_au));
        row(ui, "  scribe_cg", format!("{:.3}", ev.scribe_cg));
        row(ui, "delta_next", format!("{:.4}", p2.primary.delta_next));
        row(
            ui,
            "delta_2nd_next",
            format!("{:.4}", p2.primary.delta_second_next),
        );
        row(
            ui,
            "lazyscore",
            format!("{:.4}", p2.apex_lazy.apex_lazyscore),
        );
        row(
            ui,
            "lazyscore_z",
            format!("{:.4}", p2.apex_lazy.lazyscore_z),
        );
        row(ui, "npeaks", p2.counts.npeaks.to_string());
        row(
            ui,
            "ms1 summed int",
            format!("{:.1}", p2.intensities.ms1_summed_intensity),
        );
        row(
            ui,
            "ms2 summed int",
            format!("{:.1}", p2.intensities.ms2_summed_intensity),
        );
        row(
            ui,
            "rising/falling",
            format!("{}/{}", p2.counts.rising_cycles, p2.counts.falling_cycles),
        );
    });

    ui.separator();
    ui.label("11 features (each in [0,1])");
    let f = &p2.features;
    egui::Grid::new("feat_grid").num_columns(2).show(ui, |ui| {
        let row = |ui: &mut egui::Ui, k: &str, v: f32| {
            ui.label(k);
            ui.label(format!("{v:.3}"));
            ui.end_row();
        };
        row(ui, "peak_shape", f.peak_shape);
        row(ui, "ratio_cv", f.ratio_cv);
        row(ui, "centered_apex", f.centered_apex);
        row(ui, "precursor_coelution", f.precursor_coelution);
        row(ui, "fragment_coverage", f.fragment_coverage);
        row(ui, "precursor_apex_match", f.precursor_apex_match);
        row(ui, "xic_quality", f.xic_quality);
        row(ui, "fragment_apex_agreement", f.fragment_apex_agreement);
        row(ui, "isotope_correlation", f.isotope_correlation);
        row(ui, "gaussian_correlation", f.gaussian_correlation);
        row(ui, "per_frag_gaussian_corr", f.per_frag_gaussian_corr);
    });
}

/// `a / b`, or 0 when `b` is 0 / non-finite (avoids NaN/inf in the readout).
fn safe_ratio(a: f32, b: f32) -> f32 {
    if b.abs() > 0.0 && b.is_finite() {
        a / b
    } else {
        0.0
    }
}
