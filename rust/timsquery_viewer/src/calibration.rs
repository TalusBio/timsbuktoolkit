//! Viewer calibration state machine and background scoring thread.
//!
//! Drives the live RT calibration panel: shuffles elution groups,
//! scores them on a background thread, feeds calibrant candidates to
//! `calibrt::CalibrationState`, and exposes snapshots to the UI via
//! an MPSC channel.

use std::sync::Arc;
use std::sync::atomic::{
    AtomicU8,
    Ordering,
};
use std::sync::mpsc::{
    self,
    Receiver,
    SyncSender,
};
use std::thread::JoinHandle;

use eframe::egui;

use calibrt::{
    CALIBRANT_WEIGHT,
    CalibrationState,
    LibraryRT,
    ObservedRTSeconds,
    RidgeSummary,
};
use timscentroid::rt_mapping::{
    MS1CycleIndex,
    RTIndex,
};
use timsquery::models::tolerance::{
    MobilityTolerance,
    MzTolerance,
    QuadTolerance,
    RtTolerance,
    Tolerance,
};
use timsquery::serde::IndexedPeaksHandle;
use timsseek::rt_calibration::{
    ResidualBlock,
    rt_tolerance_from_ridge,
};
use timsseek::scoring::apex_finding::TraceScorer;
use timsseek::scoring::extraction::build_extraction;
use timsseek::scoring::pipeline::{
    CalibrantCandidate,
    CalibrantHeap,
    CalibrationConfig,
    TOP_N_FRAGMENTS,
};

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

/// How many scored elution groups between channel snapshots. The viewer's own —
/// it governs UI refresh, nothing the search does.
const SNAPSHOT_INTERVAL: usize = 100;

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
}

/// Messages sent from the background thread to the UI.
#[derive(Debug)]
pub enum CalibrationMessage {
    /// Periodic progress snapshot.
    Snapshot {
        n_scored: usize,
        heap_len: usize,
        /// (library_rt, apex_rt)
        points: Vec<(LibraryRT<f64>, ObservedRTSeconds<f64>)>,
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
    pub n_calibrants_found: usize,
    pub heap_capacity: usize,
    pub elution_group_count: usize,
    pub derived_tolerances: Option<DerivedTolerances>,
    /// What the loaded calibration's search measured. The viewer calibrates RT
    /// only and measures no residuals, so it carries this untouched. `None`
    /// whenever the curve on screen was fit here.
    pub residuals: Option<ResidualBlock>,

    thread_handle: Option<JoinHandle<()>>,
    thread_control: Arc<AtomicU8>,
    receiver: Option<Receiver<CalibrationMessage>>,

    /// Latest calibrant points: (library_rt, apex_rt).
    pub snapshot_points: Vec<(LibraryRT<f64>, ObservedRTSeconds<f64>)>,

    /// The search's own calibration settings, so the viewer does not fit on a
    /// grid the search would not use. Defaults: the viewer has no config file.
    search: CalibrationConfig,
}

impl Default for ViewerCalibrationState {
    fn default() -> Self {
        let search = CalibrationConfig::default();
        Self {
            phase: CalibrationPhase::Idle,
            calibration_state: None,
            generation: 0,
            n_scored: 0,
            n_calibrants_found: 0,
            heap_capacity: search.n_calibrants,
            elution_group_count: 0,
            derived_tolerances: None,
            residuals: None,
            thread_handle: None,
            thread_control: Arc::new(AtomicU8::new(CONTROL_STOP_REQUESTED)),
            receiver: None,
            snapshot_points: Vec::new(),
            search,
        }
    }
}

impl ViewerCalibrationState {
    /// Reconstruct from a persisted snapshot (app state restore).
    pub fn from_snapshot(snapshot: Option<calibrt::CalibrationSnapshot>) -> Self {
        let Some(snapshot) = snapshot else {
            return Self::default();
        };
        if snapshot.points.is_empty() {
            return Self::default();
        }

        let snapshot_points: Vec<(LibraryRT<f64>, ObservedRTSeconds<f64>)> = snapshot
            .points
            .iter()
            .map(|p| (LibraryRT(p[0]), ObservedRTSeconds(p[1])))
            .collect();
        let n_calibrants_found = snapshot_points.len();

        let calibration_state = calibrt::CalibrationState::from_snapshot(&snapshot).ok();
        let search = CalibrationConfig::default();

        Self {
            phase: if calibration_state.is_some() {
                CalibrationPhase::Done
            } else {
                CalibrationPhase::Idle
            },
            calibration_state,
            generation: 0,
            n_scored: n_calibrants_found,
            n_calibrants_found,
            heap_capacity: search.n_calibrants,
            elution_group_count: 0,
            derived_tolerances: None,
            residuals: None,
            thread_handle: None,
            thread_control: Arc::new(AtomicU8::new(CONTROL_STOP_REQUESTED)),
            receiver: None,
            snapshot_points,
            search,
        }
    }

    /// Extract snapshot for persistence (returns None if no calibration data).
    pub fn snapshot_for_persistence(&self) -> Option<calibrt::CalibrationSnapshot> {
        if self.snapshot_points.is_empty() {
            return None;
        }
        Some(self.snapshot())
    }

    /// The calibration as a snapshot stores it: points plus grid config.
    fn snapshot(&self) -> calibrt::CalibrationSnapshot {
        calibrt::CalibrationSnapshot {
            points: self.persisted_points(),
            grid_size: self
                .calibration_state
                .as_ref()
                .map_or(self.search.grid_size, CalibrationState::grid_bins),
            lookback: self.search.dp_lookback,
        }
    }

    /// The calibrant points as a snapshot stores them: `[library_rt,
    /// observed_rt, weight]`.
    fn persisted_points(&self) -> Vec<[f64; 3]> {
        self.snapshot_points
            .iter()
            .map(|&(lib, obs)| [lib.0, obs.0, CALIBRANT_WEIGHT])
            .collect()
    }

    /// Start the calibration background thread.
    ///
    /// Requires both raw data and elution groups to be loaded.
    /// If already running, this is a no-op.
    pub fn start(&mut self, index: Arc<IndexedPeaksHandle>, elution_groups: Arc<ElutionGroupData>) {
        if self.phase == CalibrationPhase::Running {
            return;
        }

        // Increment generation to invalidate stale data.
        self.generation += 1;
        self.n_scored = 0;
        self.n_calibrants_found = 0;
        self.snapshot_points.clear();
        self.elution_group_count = elution_groups.len();

        // No calibration state yet: the grid geometry comes from the calibrant
        // points, so it is not knowable until the first snapshot arrives (see
        // `refit`). Dropping any previous run's state also keeps a stale fit
        // from being drawn over the new run's points.
        self.calibration_state = None;

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
            self.thread_control.store(CONTROL_PAUSED, Ordering::Release);
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
        // Drop the receiver BEFORE joining — if the background thread is
        // blocked on tx.send(Done), dropping the receiver unblocks it.
        self.receiver = None;
        if let Some(handle) = self.thread_handle.take() {
            let _ = handle.join();
        }
        self.phase = CalibrationPhase::Idle;
        self.n_scored = 0;
        self.n_calibrants_found = 0;
        self.snapshot_points.clear();
        self.generation += 1;
        if let Some(cs) = &mut self.calibration_state {
            cs.reset();
        }
    }

    /// Drain the channel and update internal state.
    ///
    /// Whether the background thread is still active (Running or Paused).
    pub fn is_active(&self) -> bool {
        matches!(
            self.phase,
            CalibrationPhase::Running | CalibrationPhase::Paused
        )
    }

    /// Returns `true` if any new data was received (caller should
    /// `request_repaint`).
    pub fn poll(&mut self) -> bool {
        let Some(rx) = &self.receiver else {
            return false;
        };

        let mut changed = false;
        let mut new_points = false;

        loop {
            match rx.try_recv() {
                Ok(CalibrationMessage::Snapshot {
                    n_scored,
                    heap_len,
                    points,
                }) => {
                    self.n_scored = n_scored;
                    self.n_calibrants_found = heap_len;
                    self.snapshot_points = points;
                    new_points = true;
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

        // Each snapshot is the full heap, not a delta, so only the last one
        // drained is worth fitting.
        if new_points {
            self.refit();
        }

        changed
    }

    /// Re-fit the curve to `self.snapshot_points`, on a grid the points themselves
    /// span rather than the file's acquisition RT range — which would clamp an
    /// iRT-scaled library, whose RTs fall entirely outside it, into one edge column.
    ///
    /// A geometry the points cannot support leaves the previous fit alone: a later
    /// snapshot with more calibrants may well span a usable range.
    fn refit(&mut self) {
        let lookback = self.search.dp_lookback;
        let bins = self
            .calibration_state
            .as_ref()
            .map_or(self.search.grid_size, CalibrationState::grid_bins);
        let cs = match self.calibration_state.as_mut() {
            Some(cs) => cs,
            None => match CalibrationState::deferred(bins, lookback) {
                Ok(cs) => self.calibration_state.insert(cs),
                Err(e) => {
                    tracing::warn!("Calibration refit skipped: no grid at {bins} bins: {e:?}");
                    return;
                }
            },
        };

        let points = self
            .snapshot_points
            .iter()
            .map(|&(lib_rt, apex_rt)| (lib_rt.0, apex_rt.0));
        let (x_range, y_range) = match cs.refit(bins, points, &mut (), calibrt::ObserveOpts::NONE) {
            Ok(ranges) => ranges,
            Err(e) => {
                let n = self.snapshot_points.len();
                tracing::warn!("Calibration refit skipped over {n} points: {e:?}");
                return;
            }
        };

        let n_retained = cs
            .grid_cells()
            .iter()
            .filter(|n| !n.suppressed && n.center.weight > 0.0)
            .count();
        tracing::info!(
            "Calibration refit: scored={} calibrants={} retained_cells={} path_nodes={} curve={} x={:?} y={:?}",
            self.n_scored,
            self.n_calibrants_found,
            n_retained,
            cs.path_indices().len(),
            cs.curve().is_some(),
            x_range,
            y_range,
        );
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
        use rayon::prelude::*;

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
        let n_cycles = cycle_mapping.len();

        let mut heap = CalibrantHeap::new(heap_capacity);
        let mut n_scored: usize = 0;

        // Process in chunks — each chunk is parallelized via Rayon.
        // Between chunks: merge heaps, send snapshot, check pause/stop.
        for chunk in indices.chunks(SNAPSHOT_INTERVAL) {
            // Check control flag between chunks.
            loop {
                let flag = control.load(Ordering::Acquire);
                match flag {
                    CONTROL_RUNNING => break,
                    CONTROL_PAUSED => {
                        std::thread::park();
                        continue;
                    }
                    _ => {
                        if let Err(e) = tx.try_send(CalibrationMessage::Done { n_scored }) {
                            tracing::warn!("Calibration thread: failed to send Done on stop: {e}");
                        }
                        return;
                    }
                }
            }

            // Score chunk in parallel — per-thread TraceScorer + CalibrantHeap.
            let chunk_heap: CalibrantHeap = chunk
                .par_iter()
                .fold(
                    || {
                        (
                            // Viewer is interactive, not a hot path — a conservative
                            // default capacity is fine; realloc on outliers is free.
                            TraceScorer::new(n_cycles, 16),
                            CalibrantHeap::new(heap_capacity),
                        )
                    },
                    |(mut scorer, mut local_heap), &eg_idx| {
                        let Ok((elution_group, expected_intensities)) =
                            elution_groups.get_elem(eg_idx)
                        else {
                            return (scorer, local_heap);
                        };

                        let extraction = match build_extraction(
                            &elution_group,
                            expected_intensities,
                            index.as_ref(),
                            &tolerance,
                            Some(TOP_N_FRAGMENTS),
                        ) {
                            Ok(ext) => ext,
                            Err(_) => return (scorer, local_heap),
                        };

                        if scorer.compute_traces(&extraction).is_err() {
                            return (scorer, local_heap);
                        }

                        let cycle_offset = extraction.chromatograms.cycle_offset();
                        let rt_mapper = |idx: usize| -> u32 {
                            cycle_mapping
                                .rt_milis_for_index(&MS1CycleIndex::new(
                                    (idx + cycle_offset) as u32,
                                ))
                                .unwrap_or(0)
                        };

                        if let Ok(apex) = scorer.suggest_apex(&rt_mapper, 0) {
                            // Viewer has no skip counter wired up; NaN scores
                            // are rare and silently dropped here.
                            let _ = local_heap.push(CalibrantCandidate {
                                score: apex.score,
                                apex_rt: ObservedRTSeconds(apex.retention_time_ms as f32 / 1000.0),
                                speclib_index: eg_idx,
                                library_rt: LibraryRT(elution_group.rt_seconds()),
                            });
                        }
                        (scorer, local_heap)
                    },
                )
                .map(|(_, local_heap)| local_heap)
                .reduce(|| CalibrantHeap::new(heap_capacity), CalibrantHeap::merge);

            // Merge chunk results into main heap.
            heap = heap.merge(chunk_heap);
            n_scored += chunk.len();

            // Send snapshot.
            let points: Vec<(LibraryRT<f64>, ObservedRTSeconds<f64>)> = heap
                .iter()
                .map(|c| {
                    (
                        LibraryRT(c.library_rt.0 as f64),
                        ObservedRTSeconds(c.apex_rt.0 as f64),
                    )
                })
                .collect();

            let _ = tx.try_send(CalibrationMessage::Snapshot {
                n_scored,
                heap_len: heap.len(),
                points,
            });
        }

        if let Err(e) = tx.try_send(CalibrationMessage::Done { n_scored }) {
            tracing::warn!("Calibration thread: failed to send Done on completion: {e}");
        }
    }

    // -----------------------------------------------------------------------
    // Save / Load
    // -----------------------------------------------------------------------

    /// Write the current calibration where the search's reader can load it.
    ///
    /// The points and their geometry come from the fitted grid, so a reader
    /// refits the same curve and the same ridge widths the viewer is showing.
    /// The residual statistics are the ones this file was loaded with, since the
    /// viewer measures no residuals of its own — see [`Self::residuals`].
    pub fn save_to_file(
        &self,
        path: &std::path::Path,
        rt_range_seconds: [f64; 2],
    ) -> Result<(), String> {
        use timsseek::rt_calibration::SavedCalibration;

        let tol = self.derived_tolerances.as_ref();
        SavedCalibration::new(
            rt_range_seconds,
            self.snapshot(),
            tol.map_or(0.0, |t| t.rt_tolerance_minutes),
            self.residuals.clone(),
            self.n_scored,
        )
        .write(path)
    }

    /// Deserialize calibration state from a saved calibration, returning the
    /// provenance warning if the file may not belong to the loaded raw run.
    pub fn load_from_file(
        &mut self,
        path: &std::path::Path,
        raw_rt_range: Option<[f64; 2]>,
    ) -> Result<Option<String>, String> {
        use timsseek::rt_calibration::SavedCalibration;

        let (saved, warning) = SavedCalibration::read(path, raw_rt_range)?;

        // Reconstruct CalibrationState from the snapshot
        if let Ok(cal) = calibrt::CalibrationState::from_snapshot(&saved.calibration) {
            self.snapshot_points = saved
                .calibration
                .points
                .iter()
                .map(|p| (LibraryRT(p[0]), ObservedRTSeconds(p[1])))
                .collect();
            self.calibration_state = Some(cal);
        }

        self.n_calibrants_found = saved.n_calibrants();
        self.n_scored = saved.n_scored;
        self.derived_tolerances = Some(DerivedTolerances {
            rt_tolerance_minutes: saved.rt_tolerance_minutes,
        });
        self.residuals = saved.residuals;
        self.phase = CalibrationPhase::Done;

        Ok(warning)
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
        selected_library_rt: Option<f64>,
    ) {
        // -- Control buttons --------------------------------------------------
        ui.horizontal(|ui| {
            match self.phase {
                CalibrationPhase::Idle => {
                    let both_loaded =
                        matches!(indexed_data, crate::app::IndexedDataState::Loaded { .. })
                            && matches!(
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
                    if ui.button("Load").clicked()
                        && let Some(path) = rfd::FileDialog::new()
                            .add_filter("JSON", &["json"])
                            .pick_file()
                    {
                        let raw_rt_range = if let crate::app::IndexedDataState::Loaded {
                            index,
                            ..
                        } = indexed_data
                        {
                            let cycle_mapping = index.ms1_cycle_mapping();
                            let (rt_min_ms, rt_max_ms) = cycle_mapping.range_milis();
                            Some([rt_min_ms as f64 / 1000.0, rt_max_ms as f64 / 1000.0])
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
                        save_calibration_dialog(self, indexed_data);
                    }
                }
                CalibrationPhase::Done => {
                    if ui.button("\u{21BA} Reset").clicked() {
                        self.reset();
                    }
                    if ui.button("Save").clicked() {
                        save_calibration_dialog(self, indexed_data);
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
            ui.label(format!("Scored: {} / {}", self.n_scored, total));
            ui.separator();
            ui.label(format!(
                "Calibrants: {} / {}",
                self.n_calibrants_found, self.heap_capacity
            ));
        });

        ui.add_space(4.0);
        ui.separator();
        ui.add_space(4.0);

        // -- Tolerance suggestion (pinned to bottom, reserves its natural height) --
        egui::TopBottomPanel::bottom("calibration_footer").show_inside(ui, |ui| {
            self.render_tolerance_suggestion(ui, tolerance);
        });

        // -- Grid + curve plot (fills remaining space) -------------------------
        self.render_calibration_plot(ui, selected_library_rt);
    }

    /// Render the scatter + curve calibration plot.
    fn render_calibration_plot(&mut self, ui: &mut egui::Ui, selected_library_rt: Option<f64>) {
        use egui_plot::{
            HLine,
            Line,
            Plot,
            PlotPoints,
            Points,
            Polygon,
            VLine,
        };

        let plot_id = format!("calibration_plot_{}", self.generation);
        let plot = Plot::new(plot_id)
            .height(ui.available_height().max(100.0))
            .x_axis_label("Library RT (s)")
            .y_axis_label("Measured RT (s)")
            .allow_zoom(true)
            .allow_drag(true);

        let cal_state = self.calibration_state.as_mut();

        plot.show(ui, |plot_ui| {
            // Grid heatmap from CalibrationState (if available).
            if let Some(cs) = cal_state {
                let cells = cs.grid_cells();
                let path_indices = cs.path_indices();
                let bins = cs.grid_bins();
                let (x_lo, x_hi) = cs.grid_x_range();
                let (y_lo, y_hi) = cs.grid_y_range();
                let cell_w = (x_hi - x_lo) / bins as f64;
                let cell_h = (y_hi - y_lo) / bins as f64;

                // Find max weight for color normalization (log scale)
                let max_weight = cells
                    .iter()
                    .map(|n| n.center.weight)
                    .fold(0.0f64, f64::max)
                    .max(1.0);
                let log_max = max_weight.ln_1p();

                // Draw each non-zero cell as a colored rectangle
                for (i, node) in cells.iter().enumerate() {
                    if node.center.weight <= 0.0 {
                        continue;
                    }

                    let gx = i % bins;
                    let gy = i / bins;
                    let cx = x_lo + (gx as f64 + 0.5) * cell_w;
                    let cy = y_lo + (gy as f64 + 0.5) * cell_h;
                    let hw = cell_w * 0.5;
                    let hh = cell_h * 0.5;

                    // Log-scale color: dark blue → bright yellow
                    let t = (node.center.weight.ln_1p() / log_max) as f32;
                    let color = if node.suppressed {
                        // Suppressed: gray tones
                        let v = (40.0 + t * 80.0) as u8;
                        egui::Color32::from_rgba_unmultiplied(v, v, v, 180)
                    } else {
                        // Retained: blue → cyan → yellow heat
                        let r = (t * 255.0) as u8;
                        let g = (t * 200.0 + 55.0) as u8;
                        let b = ((1.0 - t) * 200.0) as u8;
                        egui::Color32::from_rgba_unmultiplied(r, g, b, 200)
                    };

                    let rect = vec![
                        [cx - hw, cy - hh],
                        [cx + hw, cy - hh],
                        [cx + hw, cy + hh],
                        [cx - hw, cy + hh],
                    ];
                    plot_ui.polygon(
                        Polygon::new(format!("cell_{i}"), PlotPoints::new(rect))
                            .fill_color(color)
                            .stroke(egui::Stroke::new(0.0_f32, egui::Color32::TRANSPARENT)),
                    );
                }

                // Path nodes: bright green dots on top of heatmap
                let path_pts: Vec<[f64; 2]> = path_indices
                    .iter()
                    .filter_map(|&idx| cells.get(idx))
                    .map(|n| [n.center.library, n.center.observed])
                    .collect();

                if !path_pts.is_empty() {
                    plot_ui.points(
                        Points::new("path", PlotPoints::new(path_pts))
                            .color(egui::Color32::from_rgb(50, 255, 50))
                            .radius(5.0_f32),
                    );
                }

                // Fitted curve + extrapolation over full grid range
                if let Some(curve) = cs.curve().cloned() {
                    let curve_points = curve.points();
                    if curve_points.len() >= 2 {
                        let curve_x_min = curve_points.first().unwrap().library;
                        let curve_x_max = curve_points.last().unwrap().library;
                        let (grid_x_min, grid_x_max) = cs.grid_x_range();
                        let n_samples = 200;

                        // Interpolated region (solid cyan) — within curve bounds
                        let interp_step = (curve_x_max - curve_x_min) / n_samples as f64;
                        let interp_pts: Vec<[f64; 2]> = (0..=n_samples)
                            .filter_map(|i| {
                                let x = curve_x_min + i as f64 * interp_step;
                                let y = curve.predict(LibraryRT(x)).ok()?.0;
                                Some([x, y])
                            })
                            .collect();

                        if !interp_pts.is_empty() {
                            plot_ui.line(
                                Line::new("fitted curve", PlotPoints::new(interp_pts))
                                    .color(egui::Color32::from_rgb(0, 220, 220))
                                    .width(2.0_f32),
                            );
                        }

                        // Extrapolated regions (dashed red) — beyond curve bounds
                        // Clamp Y to grid range so extrapolation doesn't fly off
                        let (grid_y_min, grid_y_max) = cs.grid_y_range();
                        let extrap_color = egui::Color32::from_rgb(255, 100, 100);
                        let extrap_predict = |x: f64| -> f64 {
                            let y = match curve.predict(LibraryRT(x)) {
                                Ok(y) => y.0,
                                Err(calibrt::CalibRtError::OutOfBounds(y)) => y,
                                Err(_) => 0.0,
                            };
                            y.clamp(grid_y_min, grid_y_max)
                        };
                        // Left extrapolation
                        if grid_x_min < curve_x_min {
                            let left_step = (curve_x_min - grid_x_min) / 50.0_f64;
                            let left_pts: Vec<[f64; 2]> = (0..=50)
                                .map(|i| {
                                    let x = grid_x_min + i as f64 * left_step;
                                    [x, extrap_predict(x)]
                                })
                                .collect();
                            plot_ui.line(
                                Line::new("extrapolation (left)", PlotPoints::new(left_pts))
                                    .color(extrap_color)
                                    .width(1.5_f32)
                                    .style(egui_plot::LineStyle::dashed_dense()),
                            );
                        }
                        // Right extrapolation
                        if grid_x_max > curve_x_max {
                            let right_step = (grid_x_max - curve_x_max) / 50.0_f64;
                            let right_pts: Vec<[f64; 2]> = (0..=50)
                                .map(|i| {
                                    let x = curve_x_max + i as f64 * right_step;
                                    [x, extrap_predict(x)]
                                })
                                .collect();
                            plot_ui.line(
                                Line::new("extrapolation (right)", PlotPoints::new(right_pts))
                                    .color(extrap_color)
                                    .width(1.5_f32)
                                    .style(egui_plot::LineStyle::dashed_dense()),
                            );
                        }

                        // Ridge envelope: upper and lower boundary lines showing tolerance width
                        let ridge = cs.ridge_widths();
                        if ridge.len() >= 2 {
                            let ridge_color =
                                egui::Color32::from_rgba_unmultiplied(0, 220, 220, 100);

                            let upper: Vec<[f64; 2]> = ridge
                                .iter()
                                .filter_map(|m| {
                                    let y = match curve.predict(m.library) {
                                        Ok(y) => y.0,
                                        Err(calibrt::CalibRtError::OutOfBounds(y)) => y,
                                        Err(_) => return None,
                                    };
                                    Some([m.library.0, y + m.half_width])
                                })
                                .collect();
                            let lower: Vec<[f64; 2]> = ridge
                                .iter()
                                .filter_map(|m| {
                                    let y = match curve.predict(m.library) {
                                        Ok(y) => y.0,
                                        Err(calibrt::CalibRtError::OutOfBounds(y)) => y,
                                        Err(_) => return None,
                                    };
                                    Some([m.library.0, y - m.half_width])
                                })
                                .collect();

                            if upper.len() >= 2 {
                                plot_ui.line(
                                    Line::new("ridge upper", PlotPoints::new(upper))
                                        .color(ridge_color)
                                        .width(1.5_f32)
                                        .style(egui_plot::LineStyle::dashed_dense()),
                                );
                            }
                            if lower.len() >= 2 {
                                plot_ui.line(
                                    Line::new("ridge lower", PlotPoints::new(lower))
                                        .color(ridge_color)
                                        .width(1.5_f32)
                                        .style(egui_plot::LineStyle::dashed_dense()),
                                );
                            }
                        }
                    }
                }

                // Selected peptide overlay: vertical line at library RT,
                // horizontal line at predicted measured RT, tolerance band
                if let Some(lib_rt) = selected_library_rt {
                    // Vertical line: library RT (x-axis)
                    plot_ui.vline(
                        VLine::new("library RT", lib_rt)
                            .color(egui::Color32::from_rgba_unmultiplied(255, 100, 100, 160))
                            .width(1.5_f32),
                    );

                    // If curve is fitted, show predicted RT + tolerance band
                    if let Some(curve) = cs.curve() {
                        let predicted_rt = match curve.predict(LibraryRT(lib_rt)) {
                            Ok(y) => y.0,
                            Err(calibrt::CalibRtError::OutOfBounds(y)) => y,
                            Err(_) => lib_rt,
                        };

                        // Horizontal line: predicted measured RT
                        plot_ui.hline(
                            HLine::new("predicted RT", predicted_rt)
                                .color(egui::Color32::from_rgba_unmultiplied(255, 100, 100, 160))
                                .width(1.5_f32),
                        );

                        // Crosshair point at (lib_rt, predicted_rt)
                        plot_ui.points(
                            Points::new("query", PlotPoints::new(vec![[lib_rt, predicted_rt]]))
                                .color(egui::Color32::from_rgb(255, 80, 80))
                                .radius(6.0_f32),
                        );
                    }
                }
            } else if !self.snapshot_points.is_empty() {
                // Fallback: show raw calibrant points if CalibrationState
                // hasn't been built yet.
                let raw_pts: Vec<[f64; 2]> = self
                    .snapshot_points
                    .iter()
                    .map(|&(lib_rt, apex_rt)| [lib_rt.0, apex_rt.0])
                    .collect();

                plot_ui.points(
                    Points::new("calibrants", PlotPoints::new(raw_pts))
                        .color(egui::Color32::from_rgb(70, 130, 230))
                        .radius(3.0_f32),
                );
            }
        });
    }

    /// Render tolerance suggestion and Apply button.
    fn render_tolerance_suggestion(&mut self, ui: &mut egui::Ui, tolerance: &mut Tolerance) {
        let floor = self.search.min_rt_tolerance_minutes;

        // The weight-averaged half-width gives the global tolerance — heavy
        // columns count more.
        let ridge_stats = self.calibration_state.as_ref().and_then(|cs| {
            cs.curve()?; // ensure curve is fitted
            let summary = RidgeSummary::of(cs.ridge_widths())?;
            summary.weighted_half_width.is_finite().then_some(summary)
        });

        // Through the search's own rule, so the number shown is the window the
        // search would open. The search applies it per query at the interpolated
        // half-width; one number can only carry the weighted average.
        let suggested = ridge_stats.map(|stats| {
            (
                rt_tolerance_from_ridge(stats.weighted_half_width, floor),
                stats,
            )
        });

        if let Some((rt_min, _)) = suggested {
            self.derived_tolerances = Some(DerivedTolerances {
                rt_tolerance_minutes: rt_min,
            });
        }

        ui.vertical(|ui| {
            if let Some((rt_min, stats)) = suggested {
                ui.horizontal(|ui| {
                    ui.label(format!(
                        "RT tol (fallback): \u{00B1}{rt_min:.2} min, uniform over all queries",
                    ));
                    if ui.button("Apply").clicked() {
                        tolerance.rt = RtTolerance::Minutes((rt_min, rt_min));
                    }
                });
                let RidgeSummary {
                    weighted_half_width: hw_s,
                    min_half_width: min_s,
                    max_half_width: max_s,
                    n_columns: n_cols,
                    in_ridge_ratio,
                } = stats;
                let in_ridge_pct = in_ridge_ratio * 100.0;
                ui.label(
                    egui::RichText::new(format!(
                        "A real search interpolates it per query instead: ridge half-width {hw_s:.0} s weighted mean, {min_s:.0}-{max_s:.0} s over {n_cols} cols, {in_ridge_pct:.0}% in-ridge.",
                    ))
                    .weak()
                    .small(),
                );
            } else if self.phase == CalibrationPhase::Idle {
                ui.label("Start calibration to compute RT tolerance suggestion.");
            } else {
                ui.horizontal(|ui| {
                    ui.label("Collecting data...");
                    ui.spinner();
                });
            }
        });
    }
}

impl Drop for ViewerCalibrationState {
    fn drop(&mut self) {
        self.stop();
        // Drop the receiver BEFORE joining — unblocks any tx.send(Done)
        // in the background thread that would otherwise deadlock.
        self.receiver = None;
        if let Some(handle) = self.thread_handle.take() {
            let _ = handle.join();
        }
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Show a file-save dialog and write the current calibration to JSON.
fn save_calibration_dialog(
    state: &mut ViewerCalibrationState,
    indexed_data: &crate::app::IndexedDataState,
) {
    if let Some(path) = rfd::FileDialog::new()
        .set_file_name("calibration.json")
        .add_filter("JSON", &["json"])
        .save_file()
    {
        let rt_range = if let crate::app::IndexedDataState::Loaded { index, .. } = indexed_data {
            let (rt_min_ms, rt_max_ms) = index.ms1_cycle_mapping().range_milis();
            [rt_min_ms as f64 / 1000.0, rt_max_ms as f64 / 1000.0]
        } else {
            [0.0, 0.0]
        };
        if let Err(e) = state.save_to_file(&path, rt_range) {
            tracing::error!("Failed to save calibration: {}", e);
        }
    }
}

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
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        let j = (state >> 33) as usize % (i + 1);
        indices.swap(i, j);
    }
}
