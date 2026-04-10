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
}

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
