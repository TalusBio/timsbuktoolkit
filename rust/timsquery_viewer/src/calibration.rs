//! Viewer calibration state machine and background scoring thread.
//!
//! Drives the live RT calibration panel: shuffles elution groups,
//! scores them on a background thread, feeds calibrant candidates to
//! `calibrt::CalibrationState`, and exposes snapshots to the UI via
//! an MPSC channel.

use std::sync::atomic::{AtomicU8, Ordering};
use std::sync::mpsc::{self, Receiver, SyncSender};
use std::sync::Arc;
use std::thread::JoinHandle;

use eframe::egui;

use calibrt::CalibrationState;
use timscentroid::rt_mapping::{MS1CycleIndex, RTIndex};
use timsquery::models::tolerance::{
    MobilityTolerance, MzTolerance, QuadTolerance, RtTolerance, Tolerance,
};
use timsquery::serde::IndexedPeaksHandle;
use timsseek::scoring::apex_finding::TraceScorer;
use timsseek::scoring::extraction::build_extraction;
use timsseek::scoring::pipeline::{CalibrantCandidate, CalibrantHeap};

use crate::file_loader::ElutionGroupData;

// ---------------------------------------------------------------------------
// Thread-control constants (stored in Arc<AtomicU8>)
// ---------------------------------------------------------------------------

const CONTROL_RUNNING: u8 = 0;
const CONTROL_PAUSED: u8 = 1;
const CONTROL_STOP_REQUESTED: u8 = 2;

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

/// Number of top fragments to keep per elution group during calibration scoring.
const CALIBRATION_TOP_N_FRAGMENTS: usize = 8;

/// How many scored elution groups between channel snapshots.
const SNAPSHOT_INTERVAL: usize = 100;

/// Default CalibrantHeap capacity.
const DEFAULT_HEAP_CAPACITY: usize = 2000;

/// Default calibrt grid size.
const DEFAULT_GRID_SIZE: usize = 100;

/// Default DP lookback for calibrt pathfinding.
const DEFAULT_LOOKBACK: usize = 30;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Phase of the calibration state machine.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CalibrationPhase {
    Idle,
    Running,
    Paused,
    Done,
}

/// Derived tolerance windows (placeholder for downstream use).
#[derive(Debug, Clone)]
pub struct DerivedTolerances {
    pub rt_tolerance_minutes: f32,
    pub wrmse: f64,
}

// Save/load uses shared types from timsseek::rt_calibration:
// SavedCalibration, SavedTolerances, LoadedCalibration, CalibrationSnapshot

/// Messages sent from the background thread to the UI.
#[derive(Debug)]
pub enum CalibrationMessage {
    /// Periodic progress snapshot.
    Snapshot {
        n_scored: usize,
        heap_len: usize,
        /// (library_rt_seconds, apex_rt_seconds, score)
        points: Vec<(f64, f64, f64)>,
    },
    /// Thread completed (all elution groups scored or stopped).
    Done { n_scored: usize },
}

/// Broad tolerance used for calibration extraction queries.
/// Wider than a typical search to ensure calibrants are found even
/// before accurate RT calibration is available.
fn broad_calibration_tolerance() -> Tolerance {
    Tolerance {
        ms: MzTolerance::Ppm((15.0, 15.0)),
        rt: RtTolerance::Unrestricted,
        mobility: MobilityTolerance::Pct((5.0, 5.0)),
        quad: QuadTolerance::Absolute((0.1, 0.1)),
    }
}

// ---------------------------------------------------------------------------
// ViewerCalibrationState
// ---------------------------------------------------------------------------

/// Owns the calibration background thread and collects results.
pub struct ViewerCalibrationState {
    pub phase: CalibrationPhase,
    pub calibration_state: Option<CalibrationState>,
    pub generation: u64,
    pub n_scored: usize,
    pub n_calibrants: usize,
    pub heap_capacity: usize,
    pub elution_group_count: usize,
    pub derived_tolerances: Option<DerivedTolerances>,

    thread_handle: Option<JoinHandle<()>>,
    thread_control: Arc<AtomicU8>,
    receiver: Option<Receiver<CalibrationMessage>>,

    /// Latest calibrant points: (library_rt, apex_rt, score).
    pub snapshot_points: Vec<(f64, f64, f64)>,
}

impl Default for ViewerCalibrationState {
    fn default() -> Self {
        Self {
            phase: CalibrationPhase::Idle,
            calibration_state: None,
            generation: 0,
            n_scored: 0,
            n_calibrants: 0,
            heap_capacity: DEFAULT_HEAP_CAPACITY,
            elution_group_count: 0,
            derived_tolerances: None,
            thread_handle: None,
            thread_control: Arc::new(AtomicU8::new(CONTROL_STOP_REQUESTED)),
            receiver: None,
            snapshot_points: Vec::new(),
        }
    }
}

impl ViewerCalibrationState {
    /// Start the calibration background thread.
    ///
    /// Requires both raw data and elution groups to be loaded.
    /// If already running, this is a no-op.
    pub fn start(
        &mut self,
        index: Arc<IndexedPeaksHandle>,
        elution_groups: Arc<ElutionGroupData>,
    ) {
        if self.phase == CalibrationPhase::Running {
            return;
        }

        // Increment generation to invalidate stale data.
        self.generation += 1;
        self.n_scored = 0;
        self.n_calibrants = 0;
        self.snapshot_points.clear();
        self.elution_group_count = elution_groups.len();

        // Build calibration state with RT range from the data.
        let cycle_mapping = index.ms1_cycle_mapping();
        let (rt_min_ms, rt_max_ms) = cycle_mapping.range_milis();
        let rt_min_sec = rt_min_ms as f64 / 1000.0;
        let rt_max_sec = rt_max_ms as f64 / 1000.0;

        self.calibration_state = CalibrationState::new(
            DEFAULT_GRID_SIZE,
            (rt_min_sec, rt_max_sec),
            (rt_min_sec, rt_max_sec),
            DEFAULT_LOOKBACK,
        )
        .ok();

        // Set up channel and control flag.
        let (tx, rx) = mpsc::sync_channel::<CalibrationMessage>(1);
        self.receiver = Some(rx);

        let control = Arc::new(AtomicU8::new(CONTROL_RUNNING));
        self.thread_control = control.clone();

        let heap_capacity = self.heap_capacity;

        let handle = std::thread::Builder::new()
            .name("calibration-bg".into())
            .spawn(move || {
                Self::background_loop(index, elution_groups, tx, control, heap_capacity);
            })
            .expect("Failed to spawn calibration thread");

        self.thread_handle = Some(handle);
        self.phase = CalibrationPhase::Running;
    }

    /// Pause the background thread (it will park until resumed).
    pub fn pause(&mut self) {
        if self.phase == CalibrationPhase::Running {
            self.thread_control
                .store(CONTROL_PAUSED, Ordering::Release);
            self.phase = CalibrationPhase::Paused;
        }
    }

    /// Resume a paused background thread.
    pub fn resume(&mut self) {
        if self.phase == CalibrationPhase::Paused {
            self.thread_control
                .store(CONTROL_RUNNING, Ordering::Release);
            if let Some(handle) = &self.thread_handle {
                handle.thread().unpark();
            }
            self.phase = CalibrationPhase::Running;
        }
    }

    /// Request the background thread to stop.
    pub fn stop(&mut self) {
        self.thread_control
            .store(CONTROL_STOP_REQUESTED, Ordering::Release);
        // Unpark in case the thread is parked.
        if let Some(handle) = &self.thread_handle {
            handle.thread().unpark();
        }
    }

    /// Stop and reset all state. Returns to Idle.
    pub fn reset(&mut self) {
        self.stop();
        // Wait for the thread to finish (non-blocking check; if it takes
        // too long the Drop impl will also try).
        if let Some(handle) = self.thread_handle.take() {
            let _ = handle.join();
        }
        self.receiver = None;
        self.phase = CalibrationPhase::Idle;
        self.n_scored = 0;
        self.n_calibrants = 0;
        self.snapshot_points.clear();
        self.generation += 1;
        if let Some(cs) = &mut self.calibration_state {
            cs.reset();
        }
    }

    /// Drain the channel and update internal state.
    ///
    /// Returns `true` if any new data was received (caller should
    /// `request_repaint`).
    pub fn poll(&mut self) -> bool {
        let Some(rx) = &self.receiver else {
            return false;
        };

        let mut changed = false;

        loop {
            match rx.try_recv() {
                Ok(CalibrationMessage::Snapshot {
                    n_scored,
                    heap_len,
                    points,
                }) => {
                    self.n_scored = n_scored;
                    self.n_calibrants = heap_len;
                    self.snapshot_points = points;

                    // Feed points into CalibrationState for curve fitting.
                    if let Some(cs) = &mut self.calibration_state {
                        cs.update(
                            self.snapshot_points
                                .iter()
                                .map(|&(lib_rt, apex_rt, score)| (lib_rt, apex_rt, score)),
                        );
                        cs.fit();
                    }
                    changed = true;
                }
                Ok(CalibrationMessage::Done { n_scored }) => {
                    self.n_scored = n_scored;
                    self.phase = CalibrationPhase::Done;
                    // Clean up thread handle.
                    if let Some(handle) = self.thread_handle.take() {
                        let _ = handle.join();
                    }
                    self.receiver = None;
                    changed = true;
                    break;
                }
                Err(mpsc::TryRecvError::Empty) => break,
                Err(mpsc::TryRecvError::Disconnected) => {
                    // Thread exited unexpectedly.
                    self.phase = CalibrationPhase::Done;
                    self.thread_handle = None;
                    self.receiver = None;
                    changed = true;
                    break;
                }
            }
        }

        changed
    }

    // -----------------------------------------------------------------------
    // Background thread
    // -----------------------------------------------------------------------

    fn background_loop(
        index: Arc<IndexedPeaksHandle>,
        elution_groups: Arc<ElutionGroupData>,
        tx: SyncSender<CalibrationMessage>,
        control: Arc<AtomicU8>,
        heap_capacity: usize,
    ) {
        let n_elution_groups = elution_groups.len();
        if n_elution_groups == 0 {
            let _ = tx.send(CalibrationMessage::Done { n_scored: 0 });
            return;
        }

        // Deterministic shuffle via simple LCG (no rand dependency).
        let mut indices: Vec<usize> = (0..n_elution_groups).collect();
        simple_shuffle(&mut indices);

        let tolerance = broad_calibration_tolerance();
        let cycle_mapping = index.ms1_cycle_mapping();

        // Thread-local scorer.
        let n_cycles = cycle_mapping.len();
        let mut scorer = TraceScorer::new(n_cycles);
        let mut heap = CalibrantHeap::new(heap_capacity);

        let mut n_scored: usize = 0;
        let mut last_snapshot_heap_len: usize = 0;

        for &eg_idx in &indices {
            // Check control flag.
            loop {
                let flag = control.load(Ordering::Acquire);
                match flag {
                    CONTROL_RUNNING => break,
                    CONTROL_PAUSED => {
                        std::thread::park();
                        // Re-check after unpark.
                        continue;
                    }
                    _ => {
                        // STOP_REQUESTED or unknown.
                        let _ = tx.send(CalibrationMessage::Done { n_scored });
                        return;
                    }
                }
            }

            // Get elution group data. Skip on error.
            let Ok((elution_group, expected_intensities)) = elution_groups.get_elem(eg_idx) else {
                continue;
            };

            // Build extraction.
            let extraction = match build_extraction(
                &elution_group,
                expected_intensities,
                index.as_ref(),
                &tolerance,
                Some(CALIBRATION_TOP_N_FRAGMENTS),
            ) {
                Ok(ext) => ext,
                Err(_) => continue,
            };

            // Compute traces.
            if scorer.compute_traces(&extraction).is_err() {
                continue;
            }

            // Build RT mapper closure.
            let cycle_offset = extraction.chromatograms.cycle_offset();
            let rt_mapper = |idx: usize| -> u32 {
                cycle_mapping
                    .rt_milis_for_index(&MS1CycleIndex::new((idx + cycle_offset) as u32))
                    .unwrap_or(0)
            };

            // Suggest apex.
            let apex = match scorer.suggest_apex(&rt_mapper, 0) {
                Ok(a) => a,
                Err(_) => continue,
            };

            let candidate = CalibrantCandidate {
                score: apex.score,
                apex_rt_seconds: apex.retention_time_ms as f32 / 1000.0,
                speclib_index: eg_idx,
                library_rt_seconds: elution_group.rt_seconds(),
            };
            heap.push(candidate);
            n_scored += 1;

            // Periodic snapshot.
            if n_scored % SNAPSHOT_INTERVAL == 0 && heap.len() != last_snapshot_heap_len {
                last_snapshot_heap_len = heap.len();
                let points: Vec<(f64, f64, f64)> = heap
                    .iter()
                    .map(|c| {
                        (
                            c.library_rt_seconds as f64,
                            c.apex_rt_seconds as f64,
                            c.score as f64,
                        )
                    })
                    .collect();

                let msg = CalibrationMessage::Snapshot {
                    n_scored,
                    heap_len: heap.len(),
                    points,
                };
                // Use try_send: if the channel is full, skip this snapshot.
                let _ = tx.try_send(msg);
            }
        }

        // Final snapshot with Done marker.
        let _ = tx.send(CalibrationMessage::Done { n_scored });
    }

    // -----------------------------------------------------------------------
    // Save / Load
    // -----------------------------------------------------------------------

    /// Serialize the current calibration state to a JSON v1 file.
    /// Delegates to `CalibrationResult::save_json` format via shared serde types.
    pub fn save_to_file(
        &self,
        path: &std::path::Path,
        rt_range_seconds: [f64; 2],
    ) -> Result<(), String> {
        use timsseek::rt_calibration::{SavedCalibration, SavedTolerances};
        use calibrt::CalibrationSnapshot;

        let tol = self.derived_tolerances.as_ref();
        let saved = SavedCalibration {
            version: "v1".to_string(),
            rt_range_seconds,
            calibration: CalibrationSnapshot {
                points: self.snapshot_points.iter().map(|&(x, y, w)| [x, y, w]).collect(),
                grid_size: DEFAULT_GRID_SIZE,
                lookback: DEFAULT_LOOKBACK,
            },
            tolerances: SavedTolerances {
                rt_minutes: tol.map_or(0.0, |t| t.rt_tolerance_minutes),
                mz_ppm: [0.0, 0.0],
                mobility_pct: [0.0, 0.0],
            },
            n_calibrants: self.n_calibrants,
            n_scored: self.n_scored,
        };
        let json = serde_json::to_string_pretty(&saved).map_err(|e| e.to_string())?;
        std::fs::write(path, json).map_err(|e| e.to_string())
    }

    /// Deserialize calibration state from a JSON v1 file.
    /// Delegates to `CalibrationResult::load_json` + `CalibrationState::from_snapshot`.
    pub fn load_from_file(
        &mut self,
        path: &std::path::Path,
        raw_rt_range: Option<[f64; 2]>,
    ) -> Result<Option<String>, String> {
        use timsseek::rt_calibration::CalibrationResult;

        let loaded = CalibrationResult::load_json(path, raw_rt_range)?;

        // Reconstruct CalibrationState from the snapshot
        if let Ok(cal) = calibrt::CalibrationState::from_snapshot(&loaded.snapshot) {
            self.snapshot_points = loaded.snapshot.points
                .iter()
                .map(|p| (p[0], p[1], p[2]))
                .collect();
            self.calibration_state = Some(cal);
        }

        self.n_calibrants = loaded.n_calibrants;
        self.n_scored = loaded.n_scored;
        self.derived_tolerances = Some(DerivedTolerances {
            rt_tolerance_minutes: loaded.tolerances.rt_minutes,
            wrmse: 0.0, // WRMSE will be recomputed from snapshot on next render
        });
        self.phase = CalibrationPhase::Done;

        Ok(loaded.warning)
    }

    // -----------------------------------------------------------------------
    // UI rendering
    // -----------------------------------------------------------------------

    /// Render the calibration panel inside an egui `Ui`.
    ///
    /// `indexed_data` and `elution_groups` are needed to enable the Start
    /// button (we need both loaded). `tolerance` is written when the user
    /// clicks [Apply].
    pub fn render_panel(
        &mut self,
        ui: &mut egui::Ui,
        indexed_data: &crate::app::IndexedDataState,
        elution_groups: &crate::app::ElutionGroupState,
        tolerance: &mut Tolerance,
    ) {
        // -- Control buttons --------------------------------------------------
        ui.horizontal(|ui| {
            match self.phase {
                CalibrationPhase::Idle => {
                    let both_loaded = matches!(
                        indexed_data,
                        crate::app::IndexedDataState::Loaded { .. }
                    ) && matches!(
                        elution_groups,
                        crate::app::ElutionGroupState::Loaded { .. }
                    );
                    if ui
                        .add_enabled(both_loaded, egui::Button::new("\u{25B6} Start"))
                        .clicked()
                    {
                        // Extract Arc handles from the loaded states.
                        if let (
                            crate::app::IndexedDataState::Loaded { index, .. },
                            crate::app::ElutionGroupState::Loaded { data, .. },
                        ) = (indexed_data, elution_groups)
                        {
                            self.start(Arc::clone(index), Arc::clone(data));
                        }
                    }
                    if ui.button("Load").clicked() {
                        if let Some(path) = rfd::FileDialog::new()
                            .add_filter("JSON", &["json"])
                            .pick_file()
                        {
                            let raw_rt_range = if let crate::app::IndexedDataState::Loaded {
                                index, ..
                            } = indexed_data
                            {
                                let cycle_mapping = index.ms1_cycle_mapping();
                                let (rt_min_ms, rt_max_ms) = cycle_mapping.range_milis();
                                Some([
                                    rt_min_ms as f64 / 1000.0,
                                    rt_max_ms as f64 / 1000.0,
                                ])
                            } else {
                                None
                            };
                            match self.load_from_file(&path, raw_rt_range) {
                                Ok(Some(warning)) => tracing::warn!("{}", warning),
                                Ok(None) => {
                                    tracing::info!("Calibration loaded from {:?}", path)
                                }
                                Err(e) => tracing::error!("Failed to load calibration: {}", e),
                            }
                        }
                    }
                }
                CalibrationPhase::Running => {
                    if ui.button("\u{23F8} Pause").clicked() {
                        self.pause();
                    }
                    if ui.button("\u{23F9} Stop").clicked() {
                        self.stop();
                    }
                }
                CalibrationPhase::Paused => {
                    if ui.button("\u{25B6} Resume").clicked() {
                        self.resume();
                    }
                    if ui.button("\u{23F9} Stop").clicked() {
                        self.stop();
                    }
                    if ui.button("Save").clicked() {
                        if let Some(path) = rfd::FileDialog::new()
                            .set_file_name("calibration.json")
                            .add_filter("JSON", &["json"])
                            .save_file()
                        {
                            let rt_range = if let crate::app::IndexedDataState::Loaded {
                                index, ..
                            } = indexed_data
                            {
                                let cycle_mapping = index.ms1_cycle_mapping();
                                let (rt_min_ms, rt_max_ms) = cycle_mapping.range_milis();
                                [rt_min_ms as f64 / 1000.0, rt_max_ms as f64 / 1000.0]
                            } else {
                                [0.0, 0.0]
                            };
                            if let Err(e) = self.save_to_file(&path, rt_range) {
                                tracing::error!("Failed to save calibration: {}", e);
                            }
                        }
                    }
                }
                CalibrationPhase::Done => {
                    if ui.button("\u{21BA} Reset").clicked() {
                        self.reset();
                    }
                    if ui.button("Save").clicked() {
                        if let Some(path) = rfd::FileDialog::new()
                            .set_file_name("calibration.json")
                            .add_filter("JSON", &["json"])
                            .save_file()
                        {
                            let rt_range = if let crate::app::IndexedDataState::Loaded {
                                index, ..
                            } = indexed_data
                            {
                                let cycle_mapping = index.ms1_cycle_mapping();
                                let (rt_min_ms, rt_max_ms) = cycle_mapping.range_milis();
                                [rt_min_ms as f64 / 1000.0, rt_max_ms as f64 / 1000.0]
                            } else {
                                [0.0, 0.0]
                            };
                            if let Err(e) = self.save_to_file(&path, rt_range) {
                                tracing::error!("Failed to save calibration: {}", e);
                            }
                        }
                    }
                }
            }
        });

        ui.add_space(4.0);

        // -- Progress counters ------------------------------------------------
        ui.horizontal(|ui| {
            let total = if self.elution_group_count > 0 {
                self.elution_group_count
            } else {
                // Fallback: show "?" until we know the total
                0
            };
            ui.label(format!(
                "Scored: {} / {}",
                self.n_scored, total
            ));
            ui.separator();
            ui.label(format!(
                "Calibrants: {} / {}",
                self.n_calibrants, self.heap_capacity
            ));
        });

        ui.add_space(4.0);
        ui.separator();
        ui.add_space(4.0);

        // -- Grid + curve plot ------------------------------------------------
        self.render_calibration_plot(ui);

        ui.add_space(4.0);
        ui.separator();
        ui.add_space(4.0);

        // -- Tolerance suggestion ---------------------------------------------
        self.render_tolerance_suggestion(ui, tolerance);
    }

    /// Render the scatter + curve calibration plot.
    fn render_calibration_plot(&self, ui: &mut egui::Ui) {
        use egui_plot::{Line, Plot, PlotPoints, Points};

        let plot_id = format!("calibration_plot_{}", self.generation);
        let plot = Plot::new(plot_id)
            .height(ui.available_height().min(400.0))
            .x_axis_label("Library RT (s)")
            .y_axis_label("Measured RT (s)")
            .allow_zoom(true)
            .allow_drag(true);

        let cal_state = self.calibration_state.as_ref();

        plot.show(ui, |plot_ui| {
            // Grid cells from CalibrationState (if available).
            if let Some(cs) = cal_state {
                let cells = cs.grid_cells();
                let path_indices = cs.path_indices();

                // Suppressed cells with any weight: small gray dots
                let suppressed_pts: Vec<[f64; 2]> = cells
                    .iter()
                    .filter(|n| n.suppressed && n.center.weight > 0.0)
                    .map(|n| [n.center.x, n.center.y])
                    .collect();

                if !suppressed_pts.is_empty() {
                    plot_ui.points(
                        Points::new(
                            "suppressed",
                            PlotPoints::new(suppressed_pts),
                        )
                        .color(egui::Color32::from_rgb(140, 140, 140))
                        .radius(2.0),
                    );
                }

                // Retained (non-suppressed) cells with weight: larger blue dots
                let retained_pts: Vec<[f64; 2]> = cells
                    .iter()
                    .filter(|n| !n.suppressed && n.center.weight > 0.0)
                    .map(|n| [n.center.x, n.center.y])
                    .collect();

                if !retained_pts.is_empty() {
                    plot_ui.points(
                        Points::new(
                            "retained",
                            PlotPoints::new(retained_pts),
                        )
                        .color(egui::Color32::from_rgb(70, 130, 230))
                        .radius(4.0),
                    );
                }

                // Path nodes: green dots
                let path_pts: Vec<[f64; 2]> = path_indices
                    .iter()
                    .filter_map(|&idx| cells.get(idx))
                    .map(|n| [n.center.x, n.center.y])
                    .collect();

                if !path_pts.is_empty() {
                    plot_ui.points(
                        Points::new(
                            "path",
                            PlotPoints::new(path_pts),
                        )
                        .color(egui::Color32::from_rgb(50, 205, 50))
                        .radius(5.0),
                    );
                }

                // Fitted curve: cyan line sampled at 200 points
                if let Some(curve) = cs.curve() {
                    let curve_points = curve.points();
                    if curve_points.len() >= 2 {
                        let x_min = curve_points.first().unwrap().x;
                        let x_max = curve_points.last().unwrap().x;
                        let n_samples = 200;
                        let step = (x_max - x_min) / n_samples as f64;

                        let line_pts: Vec<[f64; 2]> = (0..=n_samples)
                            .filter_map(|i| {
                                let x = x_min + i as f64 * step;
                                // predict returns Err for out-of-bounds but we stay in range
                                let y = match curve.predict(x) {
                                    Ok(y) => y,
                                    Err(calibrt::CalibRtError::OutOfBounds(y)) => y,
                                    Err(_) => return None,
                                };
                                Some([x, y])
                            })
                            .collect();

                        if !line_pts.is_empty() {
                            plot_ui.line(
                                Line::new(
                                    "fitted curve",
                                    PlotPoints::new(line_pts),
                                )
                                .color(egui::Color32::from_rgb(0, 220, 220))
                                .width(2.0),
                            );
                        }
                    }
                }
            } else if !self.snapshot_points.is_empty() {
                // Fallback: show raw calibrant points if CalibrationState
                // hasn't been built yet (shouldn't normally happen, but
                // keeps the plot populated).
                let raw_pts: Vec<[f64; 2]> = self
                    .snapshot_points
                    .iter()
                    .map(|&(lib_rt, apex_rt, _)| [lib_rt, apex_rt])
                    .collect();

                plot_ui.points(
                    Points::new(
                        "calibrants",
                        PlotPoints::new(raw_pts),
                    )
                    .color(egui::Color32::from_rgb(70, 130, 230))
                    .radius(3.0),
                );
            }
        });
    }

    /// Render tolerance suggestion and Apply button.
    fn render_tolerance_suggestion(&mut self, ui: &mut egui::Ui, tolerance: &mut Tolerance) {
        // Compute WRMSE if we have a curve and snapshot points.
        let wrmse = self
            .calibration_state
            .as_ref()
            .and_then(|cs| {
                let curve = cs.curve()?;
                let points: Vec<calibrt::Point> = self
                    .snapshot_points
                    .iter()
                    .map(|&(lib_rt, apex_rt, score)| calibrt::Point {
                        x: lib_rt,
                        y: apex_rt,
                        weight: score,
                    })
                    .collect();
                let val = curve.wrmse(points.iter());
                if val.is_finite() { Some(val) } else { None }
            });

        // Derive a suggested RT tolerance: 3x WRMSE in minutes.
        let suggested_rt_min = wrmse.map(|w| (w / 60.0) * 3.0);

        if let (Some(rt_min), Some(w)) = (suggested_rt_min, wrmse) {
            self.derived_tolerances = Some(DerivedTolerances {
                rt_tolerance_minutes: rt_min as f32,
                wrmse: w,
            });
        }

        ui.horizontal(|ui| {
            if let (Some(rt_min), Some(w)) = (suggested_rt_min, wrmse) {
                ui.label(format!(
                    "Suggested RT: \u{00B1}{:.2} min   WRMSE: {:.2} s",
                    rt_min, w
                ));
                if ui.button("Apply").clicked() {
                    let rt_tol = rt_min as f32;
                    tolerance.rt = RtTolerance::Minutes((rt_tol, rt_tol));
                }
            } else if self.phase == CalibrationPhase::Idle {
                ui.label("Start calibration to compute RT tolerance suggestion.");
            } else {
                ui.label("Collecting data...");
                ui.spinner();
            }
        });
    }
}

impl Drop for ViewerCalibrationState {
    fn drop(&mut self) {
        self.stop();
        if let Some(handle) = self.thread_handle.take() {
            let _ = handle.join();
        }
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Simple deterministic shuffle using a linear congruential generator.
/// Avoids pulling in the `rand` crate just for this.
fn simple_shuffle(indices: &mut [usize]) {
    let len = indices.len();
    if len <= 1 {
        return;
    }
    // LCG parameters (Numerical Recipes).
    let mut state: u64 = 0xDEAD_BEEF_CAFE_BABE;
    for i in (1..len).rev() {
        state = state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        let j = (state >> 33) as usize % (i + 1);
        indices.swap(i, j);
    }
}
