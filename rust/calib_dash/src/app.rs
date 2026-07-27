//! Dashboard app state: the tab/stage cursors, the vim-style count prefix,
//! and the pause action a keypress resolves to. [`CalibDash`] is the piece
//! that ties all of it to a running search: [`CalibDash::on_batch`] is
//! called once per Phase 1 scoring chunk.
//!
//! Counts use the vim prefix idiom so there is no modal input: `15n` advances
//! fifteen batches. Without it, "skip N" would need a second text-entry mode
//! with its own key routing, cursor and escape handling, for one integer.
//! The accumulator is a single `Option<u32>` on `App` and serves every
//! motion, not just batch stepping.

use crate::metrics::{
    churn,
    curve_delta,
    weighted_ridge_half_width,
};
use crate::{
    BatchMetrics,
    CalibrantPoint,
    FitRecording,
    FrameStore,
};
use calibrt::{
    CalibrationCurve,
    CalibrationState,
    LibraryRT,
    ObserveOpts,
    ObservedRTSeconds,
};
use ratatui::crossterm::event::{
    self,
    Event,
    KeyCode,
    KeyEvent,
    KeyEventKind,
    KeyModifiers,
};
use std::io::IsTerminal;
use std::ops::Range;

/// Top-level dashboard tabs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tab {
    Fit,
    Convergence,
    Tolerances,
}

impl Tab {
    pub const ALL: [Tab; 3] = [Tab::Fit, Tab::Convergence, Tab::Tolerances];

    pub fn title(self) -> &'static str {
        match self {
            Tab::Fit => "Fit",
            Tab::Convergence => "Convergence",
            Tab::Tolerances => "Tolerances",
        }
    }

    fn index(self) -> usize {
        Tab::ALL
            .iter()
            .position(|t| *t == self)
            .expect("Tab::ALL lists every variant")
    }

    /// One tab forward, wrapping — unlike `Stage::next`, which clamps at the
    /// end of a pipeline that has no "next stage after the last one" to wrap
    /// to. A tab bar has no natural end, so `l` past `Tolerances` is `Fit`
    /// again, the same as vim's own `l`/`h` inside a fixed-width line.
    pub fn next(self) -> Self {
        Tab::ALL[(self.index() + 1) % Tab::ALL.len()]
    }

    /// One tab back, wrapping the other direction.
    pub fn prev(self) -> Self {
        Tab::ALL[(self.index() + Tab::ALL.len() - 1) % Tab::ALL.len()]
    }
}

/// A stage in the fit pipeline the Fit tab steps through with `[` and `]`.
/// Order matches the pipeline itself: a raw grid is suppressed down to
/// row/column maxima, chained into a path, fitted into a curve, and finally
/// measured for ridge width.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stage {
    Grid,
    Suppressed,
    Path,
    Curve,
    Ridge,
}

impl Stage {
    /// Every variant, in pipeline order. Lets a count-driven caller derive
    /// how many steps could possibly matter (`ALL.len() - 1`) instead of
    /// hardcoding the variant count, which would silently go stale the next
    /// time a stage is added.
    pub const ALL: [Stage; 5] = [
        Stage::Grid,
        Stage::Suppressed,
        Stage::Path,
        Stage::Curve,
        Stage::Ridge,
    ];

    /// One step forward, clamped at `Ridge` — stepping past the end of the
    /// pipeline holds rather than wrapping back to `Grid`.
    pub fn next(self) -> Self {
        match self {
            Stage::Grid => Stage::Suppressed,
            Stage::Suppressed => Stage::Path,
            Stage::Path => Stage::Curve,
            Stage::Curve => Stage::Ridge,
            Stage::Ridge => Stage::Ridge,
        }
    }

    /// One step back, clamped at `Grid`.
    pub fn prev(self) -> Self {
        match self {
            Stage::Grid => Stage::Grid,
            Stage::Suppressed => Stage::Grid,
            Stage::Path => Stage::Suppressed,
            Stage::Curve => Stage::Path,
            Stage::Ridge => Stage::Curve,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Stage::Grid => "Grid",
            Stage::Suppressed => "Suppressed",
            Stage::Path => "Path",
            Stage::Curve => "Curve",
            Stage::Ridge => "Ridge",
        }
    }
}

/// Whether the batch loop driving Phase 1 scoring should keep going after a
/// pause. Only `Ctrl-C` ever produces `Abort` — a dashboard failure anywhere
/// else is logged and treated as `Continue`, since a dev tool must never
/// fail a search.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Flow {
    Continue,
    Abort,
}

/// What the user asked for on leaving a pause.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PauseAction {
    /// Stay in the pause — the key changed the view, not the flow.
    Stay,
    /// Advance this many batches before pausing again.
    Next(u32),
    RunToEnd,
    Detach,
    Abort,
}

/// Skip bookkeeping, split out of `App` so the render-or-skip decision is
/// testable without a terminal.
pub struct Stepper {
    /// Batches still to skip before the next pause.
    remaining: u32,
    stopped: bool,
}

impl Default for Stepper {
    fn default() -> Self {
        Self::new()
    }
}

impl Stepper {
    pub fn new() -> Self {
        Self {
            remaining: 0,
            stopped: false,
        }
    }

    /// Consulted by `on_batch` exactly once per chunk. Takes `&mut self`
    /// because asking the question is what draws down the skip counter —
    /// splitting the query from the decrement is how a caller forgets one of
    /// them and the skip never expires.
    pub fn should_render(&mut self) -> bool {
        if self.stopped {
            return false;
        }
        if self.remaining > 0 {
            self.remaining -= 1;
            return false;
        }
        true
    }

    /// Whether a prior pause resolved to `RunToEnd`, `Detach`, or `Abort` —
    /// the user has already asked to stop seeing the dashboard (or stop the
    /// search outright), and nothing later should reopen it. Read-only,
    /// unlike `should_render`: querying this must not draw down the skip
    /// counter, since a caller like `CalibDash::present` isn't asking "should
    /// I render *this* batch" at all, just "has the user already opted out
    /// of every future pause."
    pub fn is_stopped(&self) -> bool {
        self.stopped
    }

    /// Called after a pause with whatever the user chose.
    pub fn apply(&mut self, action: PauseAction) {
        match action {
            PauseAction::Stay => {}
            // `Next(5)` at the chunk-0 pause means "render at chunk 5", so
            // four chunks are skipped, not five.
            PauseAction::Next(n) => self.remaining = n.saturating_sub(1),
            PauseAction::RunToEnd | PauseAction::Detach | PauseAction::Abort => self.stopped = true,
        }
    }
}

/// Dashboard state for one pause: which tab and pipeline stage are showing,
/// the DP pane toggle, the current batch number, the recording of the fit
/// that ran at this batch, and the pending vim-style count prefix.
pub struct App {
    tab: Tab,
    stage: Stage,
    dp_pane: bool,
    batch: u32,
    recording: FitRecording,
    /// The post-Phase-2 recording, set once `calibration.json` is written.
    /// `None` for the whole of Phase 1: the Tolerances tab reads this to
    /// decide whether Step B has run yet at all. `CalibDash::show_final`
    /// wires a real value in here once Step B runs; `ui.rs` only reads it.
    real_fit: Option<FitRecording>,
    /// The Tolerances tab's post-Step-B m/z/mobility/RT summary, set
    /// alongside `real_fit` by `CalibDash::show_final` — the tab's own
    /// `None` arm already explains "Step B hasn't run"; this is the payload
    /// once it has. Lives on `App`, not `CalibDash`, for the same reason
    /// `real_fit` does: `ui.rs` only ever receives `&App`.
    tolerances: Option<ToleranceSummary>,
    /// Every batch's metrics, undecimated — see `metrics.rs`'s module doc for
    /// why the series must have no holes. `CalibDash::on_batch` pushes one
    /// entry every batch, rendered or not; `push_metrics` is also called
    /// directly by `ui.rs`'s fixtures to build a history to render without a
    /// live pipeline.
    metrics: Vec<BatchMetrics>,
    /// The three numbers the Convergence tab's header states about the frame
    /// slab's decimation. `App` itself owns no `FrameStore` — `CalibDash`
    /// does — so these stay plain fields defaulting to "nothing decimated",
    /// and `CalibDash::on_batch`/`finish` are what keep them truthful,
    /// pushing `FrameStore::len`/`keep_every`/`dropped` through
    /// `set_frame_summary` after every recorded frame. That is the
    /// reconciliation: one real source (`FrameStore`), fed through the one
    /// setter, rather than a second, parallel decimation tracker.
    retained_frames: usize,
    frame_stride: usize,
    dropped_frames: usize,
    /// Per-overlay visibility on the Fit tab. Each is independently
    /// toggleable (`set_show_*`), but changing `stage` re-derives all four
    /// from "the set implied by that stage" first — see
    /// `sync_overlays_to_stage` — so stepping `[`/`]` alone still tells the
    /// story, and an independent toggle is a deliberate override on top of
    /// that default rather than a second, disconnected source of truth.
    show_suppressed: bool,
    show_path: bool,
    show_curve: bool,
    show_ridge: bool,
    /// The pending count prefix, e.g. the `15` typed before `n` in `15n`.
    /// Lives for exactly one motion: a motion consumes it via `take_count`,
    /// and any other key clears it directly, so a half-typed number never
    /// leaks into the next command.
    count: Option<u32>,
    /// Which retained frame `<`/`>` have scrubbed the Fit tab to, or `None`
    /// for the live batch. Bounded by `retained_frames` — see `handle_key`'s
    /// `<`/`>` arms — the same way `stage` is bounded by `Stage::ALL`. `App`
    /// only tracks *which* frame; it owns no `FrameStore` to refit from, so
    /// the actual recording lives in `scrub_recording`, kept in lockstep by
    /// `CalibDash::sync_scrub` (the one piece of this that needs a
    /// `FrameStore`).
    scrub_frame: Option<usize>,
    /// The original batch/chunk number of `scrub_frame`, purely for display
    /// (the Fit tab's "you are not looking at the live batch" banner) — a
    /// retained-frame *index* alone ("frame 2 of 5") means little to a user
    /// who thinks in batch numbers.
    scrub_chunk: Option<usize>,
    /// `CalibDash::refit_frame(scrub_frame)`'s result, cloned in by
    /// `CalibDash::sync_scrub` so the Fit tab can render it the same way it
    /// renders the live `recording` — through a plain `&App`, with no
    /// `CalibDash` in reach. `None` whenever `scrub_frame` is `None`, and
    /// also `None` momentarily right after `scrub_frame` changes and before
    /// the next `sync_scrub` call catches up — `active_recording` falls back
    /// to the live view for that gap rather than showing nothing.
    scrub_recording: Option<FitRecording>,
}

impl App {
    pub fn new(bins: usize) -> Self {
        Self {
            tab: Tab::Fit,
            stage: Stage::Grid,
            dp_pane: false,
            batch: 0,
            recording: FitRecording::new(bins),
            real_fit: None,
            tolerances: None,
            metrics: Vec::new(),
            retained_frames: 0,
            frame_stride: 1,
            dropped_frames: 0,
            show_suppressed: false,
            show_path: false,
            show_curve: false,
            show_ridge: false,
            count: None,
            scrub_frame: None,
            scrub_chunk: None,
            scrub_recording: None,
        }
    }

    /// Switches tabs directly, bypassing `handle_key`. `handle_key` does not
    /// yet route a tab-switching key (that lands with the rest of the key
    /// map in the next task), so this is how both tests and, later, that key
    /// handler change the tab.
    pub fn set_tab(&mut self, tab: Tab) {
        self.tab = tab;
    }

    /// Sets the pipeline stage directly, bypassing the `[`/`]` step count,
    /// and re-derives the four overlay toggles from it (see
    /// `sync_overlays_to_stage`). Kept separate from `handle_key`'s stepping
    /// so a test (or a future jump-to-stage key) can land on a specific
    /// stage in one call.
    pub fn set_stage(&mut self, stage: Stage) {
        self.stage = stage;
        self.sync_overlays_to_stage();
    }

    /// Resets all four overlay toggles to the set `self.stage` implies —
    /// cumulative, so `Stage::Ridge` implies all four are on. Called
    /// whenever `stage` changes (`set_stage`, and `handle_key`'s `[`/`]`),
    /// so stepping the stage alone still tells the whole story; an
    /// independent `set_show_*` call after that is an explicit override,
    /// not fighting a second default.
    fn sync_overlays_to_stage(&mut self) {
        self.show_suppressed = matches!(
            self.stage,
            Stage::Suppressed | Stage::Path | Stage::Curve | Stage::Ridge
        );
        self.show_path = matches!(self.stage, Stage::Path | Stage::Curve | Stage::Ridge);
        self.show_curve = matches!(self.stage, Stage::Curve | Stage::Ridge);
        self.show_ridge = matches!(self.stage, Stage::Ridge);
    }

    /// Mutable access to the live recording, so a caller — `CalibDash::on_batch`,
    /// or a test fixture here — can fit directly into the allocation `App`
    /// already owns rather than building a `FitRecording` elsewhere and
    /// having nowhere to put it.
    pub fn recording_mut(&mut self) -> &mut FitRecording {
        &mut self.recording
    }

    /// The post-Phase-2 recording. `None` throughout Phase 1 — the
    /// Tolerances tab's whole reason for existing is to explain that rather
    /// than render an empty panel.
    pub fn real_fit(&self) -> Option<&FitRecording> {
        self.real_fit.as_ref()
    }

    /// The Tolerances tab's post-Step-B summary. `None` exactly when
    /// `real_fit` is — `CalibDash::show_final` always sets both together.
    pub fn tolerances(&self) -> Option<&ToleranceSummary> {
        self.tolerances.as_ref()
    }

    /// Wires a post-Step-B recording and tolerance summary directly onto
    /// `App`. `CalibDash::show_final` is the real path that reaches this
    /// state from a live search; this is the one setter both it and a
    /// `ui.rs` fixture (there is no live `CalibDash` behind a fixture) go
    /// through, the same "one setter for the real caller and tests" shape as
    /// `set_frame_summary`.
    pub fn set_final(&mut self, recording: FitRecording, tolerances: ToleranceSummary) {
        self.real_fit = Some(recording);
        self.tolerances = Some(tolerances);
    }

    pub fn metrics(&self) -> &[BatchMetrics] {
        &self.metrics
    }

    /// Appends one batch's metrics. `CalibDash::on_batch` calls this every
    /// batch, rendered or not; fixtures call it directly here.
    pub fn push_metrics(&mut self, m: BatchMetrics) {
        self.metrics.push(m);
    }

    pub fn retained_frames(&self) -> usize {
        self.retained_frames
    }

    pub fn frame_stride(&self) -> usize {
        self.frame_stride
    }

    pub fn dropped_frames(&self) -> usize {
        self.dropped_frames
    }

    /// Sets the three numbers the Convergence tab's header states about the
    /// frame slab's decimation. `CalibDash::on_batch`/`finish` call this with
    /// the real `FrameStore`'s own numbers every time a frame is recorded;
    /// `ui.rs`'s fixtures call it directly to build a history without a live
    /// `FrameStore` to read from.
    pub fn set_frame_summary(&mut self, retained: usize, stride: usize, dropped: usize) {
        self.retained_frames = retained;
        self.frame_stride = stride;
        self.dropped_frames = dropped;
    }

    pub fn tab(&self) -> Tab {
        self.tab
    }

    pub fn stage(&self) -> Stage {
        self.stage
    }

    pub fn batch(&self) -> u32 {
        self.batch
    }

    pub fn dp_pane(&self) -> bool {
        self.dp_pane
    }

    pub fn show_suppressed(&self) -> bool {
        self.show_suppressed
    }

    pub fn show_path(&self) -> bool {
        self.show_path
    }

    pub fn show_curve(&self) -> bool {
        self.show_curve
    }

    pub fn show_ridge(&self) -> bool {
        self.show_ridge
    }

    pub fn set_show_suppressed(&mut self, on: bool) {
        self.show_suppressed = on;
    }

    pub fn set_show_path(&mut self, on: bool) {
        self.show_path = on;
    }

    pub fn set_show_curve(&mut self, on: bool) {
        self.show_curve = on;
    }

    pub fn set_show_ridge(&mut self, on: bool) {
        self.show_ridge = on;
    }

    pub fn recording(&self) -> &FitRecording {
        &self.recording
    }

    /// What the Fit tab should actually draw: the scrubbed frame's recording
    /// while `scrub_frame` is set and `sync_scrub` has already caught up to
    /// it, the live `recording` otherwise. Never panics and never shows a
    /// blank tab for the gap between `handle_key` moving `scrub_frame` and
    /// the next `sync_scrub` — it just falls back to live for that one
    /// frame.
    pub fn active_recording(&self) -> &FitRecording {
        if self.scrub_frame.is_some()
            && let Some(rec) = self.scrub_recording.as_ref()
        {
            rec
        } else {
            &self.recording
        }
    }

    /// Which retained frame `<`/`>` have scrubbed to, or `None` for the live
    /// batch. `ui.rs` reads this to decide whether to draw the "not live"
    /// banner on the Fit tab.
    pub fn scrub_frame(&self) -> Option<usize> {
        self.scrub_frame
    }

    /// The scrubbed frame's original batch/chunk number, for the banner —
    /// `None` whenever `scrub_frame` is (`CalibDash::sync_scrub` always sets
    /// both together).
    pub fn scrub_chunk(&self) -> Option<usize> {
        self.scrub_chunk
    }

    /// Wires a freshly refit frame into the scrub view. Only ever called by
    /// `CalibDash::sync_scrub`, which is the one piece of this with a
    /// `FrameStore` to refit from — `App` cannot compute this itself.
    pub fn set_scrub_recording(
        &mut self,
        frame: usize,
        chunk: Option<usize>,
        recording: FitRecording,
    ) {
        self.scrub_frame = Some(frame);
        self.scrub_chunk = chunk;
        self.scrub_recording = Some(recording);
    }

    /// Drops back to the live view — either the user scrubbed past the last
    /// retained frame back to "now" (`handle_key`'s `>` arm), or
    /// `CalibDash::sync_scrub` could not refit the frame `scrub_frame` named
    /// (should not happen in practice, since that index is bounded by the
    /// same `FrameStore` `sync_scrub` reads, but falling back to live is
    /// always safe, unlike leaving a stale or absent recording on screen).
    pub fn clear_scrub(&mut self) {
        self.scrub_frame = None;
        self.scrub_chunk = None;
        self.scrub_recording = None;
    }

    pub fn pending_count(&self) -> Option<u32> {
        self.count
    }

    /// Consumes the pending count, defaulting to 1 for a bare motion.
    pub fn take_count(&mut self) -> u32 {
        self.count.take().unwrap_or(1)
    }

    /// Folds one more digit into the pending count. A leading `0` does not
    /// start a count — matching vim, which keeps `0` free as a motion later
    /// — but `0` extends an already-started count like any other digit.
    ///
    /// Saturating, not checked: a held-down digit key is ordinary user input,
    /// not a bug, so the 10th keystroke of a run must not panic (overflow
    /// checks are on in dev/test) or silently wrap (the release behavior,
    /// which is not better — a wrapped count would drive some other,
    /// unrelated number of batches). Clamping at `u32::MAX` is a count no
    /// motion below will ever exhaust anyway.
    fn push_digit(&mut self, d: u32) {
        if d == 0 && self.count.is_none() {
            return;
        }
        self.count = Some(self.count.unwrap_or(0).saturating_mul(10).saturating_add(d));
    }

    /// Consumes the pending count and clamps it to the number of steps that
    /// could possibly move `stage` — `Stage::ALL.len() - 1`, derived rather
    /// than hardcoded so a stage added later cannot silently go uncapped —
    /// so a huge typed count costs no more work than a real motion ever
    /// could.
    fn max_stage_steps(&mut self) -> u32 {
        self.take_count().min(Stage::ALL.len() as u32 - 1)
    }

    /// Consumes the pending count and reduces it mod `Tab::ALL.len()` for
    /// `h`/`l`'s tab cycling. Unlike `max_stage_steps`'s clamp, `Tab` wraps
    /// rather than stopping at an end, so the count that actually matters
    /// isn't "how far can this possibly move" but "where does this land
    /// after wrapping" — `n % 3` answers that directly with no loop at all,
    /// which is a stronger anti-spin guarantee than a clamped loop: the cost
    /// is one division, never proportional to the typed count.
    fn tab_cycle_steps(&mut self) -> u32 {
        self.take_count() % Tab::ALL.len() as u32
    }

    /// Consumes the pending count and reduces it to whether an odd or even
    /// number of presses was requested. A boolean toggle applied twice is a
    /// no-op, so `"10s"` and `"2s"` and a bare `"s"` typed once each mean
    /// exactly one of "flip it" or "leave it" — the count's only observable
    /// effect on a toggle is its parity. Like `tab_cycle_steps`, this is O(1)
    /// regardless of how large the typed count is, so no loop-and-clamp is
    /// needed to keep a huge count from doing proportional work.
    fn toggle_parity(&mut self) -> bool {
        self.take_count() % 2 == 1
    }

    /// Moves the Fit tab's frame-scrub cursor back into history by `n`
    /// retained frames (`n` from `take_count`), or does nothing if no frames
    /// have been retained yet. Saturating arithmetic against `retained_frames`
    /// keeps this O(1) and panic-free for any `n`, including the huge counts
    /// `push_digit` no longer rejects — there is no loop here to spin at all,
    /// which is a stronger bound than `max_stage_steps`'s clamped one.
    /// `CalibDash::sync_scrub` is what actually turns the new `scrub_frame`
    /// into a recording to render; this only moves the cursor.
    fn scrub_back(&mut self, n: u32) {
        let total = self.retained_frames as u32;
        if total == 0 {
            return;
        }
        let pos = self.scrub_frame.map(|i| i as u32).unwrap_or(total);
        self.scrub_frame = Some(pos.saturating_sub(n) as usize);
    }

    /// Moves the cursor forward by `n`. Forward past the last retained frame
    /// returns to the live view (`scrub_frame = None`) rather than clamping
    /// at the last index — "keep pressing `>` and you get back to now" is
    /// the whole point of a cursor that can leave history behind.
    fn scrub_forward(&mut self, n: u32) {
        let Some(i) = self.scrub_frame else {
            return;
        };
        let total = self.retained_frames as u32;
        let new_pos = (i as u32).saturating_add(n);
        self.scrub_frame = if new_pos >= total {
            None
        } else {
            Some(new_pos as usize)
        };
    }

    /// Routes one keypress: digits build the count prefix, motions consume
    /// it via `take_count`, and every other key clears it. Plain-character
    /// arms are guarded on `key.modifiers.is_empty()` so a modified key —
    /// `Ctrl-C` above all — never falls into a bare-letter arm by accident.
    pub fn handle_key(&mut self, key: KeyEvent) -> PauseAction {
        if key.code == KeyCode::Char('c') && key.modifiers.contains(KeyModifiers::CONTROL) {
            self.count = None;
            return PauseAction::Abort;
        }

        match key.code {
            KeyCode::Char(c) if key.modifiers.is_empty() && c.is_ascii_digit() => {
                self.push_digit(c.to_digit(10).expect("guarded by is_ascii_digit"));
                PauseAction::Stay
            }
            // No modifier guard: Esc clears the count and stays either way,
            // which is exactly what the catch-all arm below would also do
            // for a modified Esc — unlike the `Char` arms, there is no
            // distinct action a modified Esc could be mistaken for.
            KeyCode::Esc => {
                self.count = None;
                PauseAction::Stay
            }
            KeyCode::Char('n') if key.modifiers.is_empty() => {
                let n = self.take_count();
                // `n` advances the *live* batch loop, so a scrub cursor left
                // over from browsing history at this pause must not still be
                // showing at the next one — nothing about it would be wrong
                // (the bound is still safe and the banner still names the
                // real batch), but reopening on a replayed frame after the
                // user just asked to advance reads as confusing rather than
                // useful. Advancing back to live is the more expected
                // default; `<`/`>` are still there to scrub again at the
                // next pause if that's what's wanted.
                self.clear_scrub();
                PauseAction::Next(n)
            }
            KeyCode::Char(']') if key.modifiers.is_empty() => {
                // `Stage` has a fixed, tiny state space, so anything past
                // `ALL.len() - 1` steps is wasted work — without this clamp,
                // a large typed count (reachable now that digit overflow no
                // longer panics) would loop synchronously for as long as it
                // takes to freeze the UI.
                for _ in 0..self.max_stage_steps() {
                    self.stage = self.stage.next();
                }
                self.sync_overlays_to_stage();
                PauseAction::Stay
            }
            KeyCode::Char('[') if key.modifiers.is_empty() => {
                for _ in 0..self.max_stage_steps() {
                    self.stage = self.stage.prev();
                }
                self.sync_overlays_to_stage();
                PauseAction::Stay
            }
            KeyCode::Char('d') if key.modifiers.is_empty() => {
                self.dp_pane = !self.dp_pane;
                self.count = None;
                PauseAction::Stay
            }
            // ---- tabs ----
            KeyCode::Char('l') if key.modifiers.is_empty() => {
                for _ in 0..self.tab_cycle_steps() {
                    self.tab = self.tab.next();
                }
                PauseAction::Stay
            }
            KeyCode::Char('h') if key.modifiers.is_empty() => {
                for _ in 0..self.tab_cycle_steps() {
                    self.tab = self.tab.prev();
                }
                PauseAction::Stay
            }
            // ---- Fit-tab batch scrubber ----
            KeyCode::Char('<') if key.modifiers.is_empty() => {
                let n = self.take_count();
                self.scrub_back(n);
                PauseAction::Stay
            }
            KeyCode::Char('>') if key.modifiers.is_empty() => {
                let n = self.take_count();
                self.scrub_forward(n);
                PauseAction::Stay
            }
            // ---- Fit-tab overlay toggles ----
            // Mnemonic letters chosen to avoid every binding above:
            // s(uppressed), p(ath), c(urve), w(idth — the ridge overlay is
            // the ridge *width* measurement everywhere else in this crate,
            // and `r` was already taken by run-to-end).
            KeyCode::Char('s') if key.modifiers.is_empty() => {
                if self.toggle_parity() {
                    self.show_suppressed = !self.show_suppressed;
                }
                PauseAction::Stay
            }
            KeyCode::Char('p') if key.modifiers.is_empty() => {
                if self.toggle_parity() {
                    self.show_path = !self.show_path;
                }
                PauseAction::Stay
            }
            KeyCode::Char('c') if key.modifiers.is_empty() => {
                if self.toggle_parity() {
                    self.show_curve = !self.show_curve;
                }
                PauseAction::Stay
            }
            KeyCode::Char('w') if key.modifiers.is_empty() => {
                if self.toggle_parity() {
                    self.show_ridge = !self.show_ridge;
                }
                PauseAction::Stay
            }
            KeyCode::Char('r') if key.modifiers.is_empty() => {
                self.count = None;
                PauseAction::RunToEnd
            }
            KeyCode::Char('q') if key.modifiers.is_empty() => {
                self.count = None;
                PauseAction::Detach
            }
            _ => {
                self.count = None;
                PauseAction::Stay
            }
        }
    }
}

/// Weight threshold (as a fraction of the path cell's own weight) that
/// `measure_ridge_width_with` expands out to. Matches the value
/// `timsseek_cli::processing` uses for the real Step B ridge measurement
/// (`cal_state.measure_ridge_width(0.1)`), so what a batch's ridge overlay
/// shows here is the same computation Phase 2 will run, not a
/// dashboard-only approximation.
const RIDGE_FRACTION: f64 = 0.1;

/// How many library-RT samples `curve_delta` draws between batches' curves.
/// Not load-bearing for correctness (any `>= 2` produces a defined delta);
/// chosen high enough that a local kink between two batches' curves is
/// unlikely to fall between sample points.
const CURVE_DELTA_SAMPLES: usize = 50;

/// The min/max of `library_rt` and `observed_rt` over a point set, exactly
/// the fold `timsseek_cli::processing` runs over its own `Point` set before
/// building the real `CalibrationState` for a batch. `on_batch`'s live re-fit
/// and `refit_frame`'s replay both call this — the same function, not two
/// copies of the same fold — because a divergence here is exactly what would
/// make a scrubbed frame's heatmap show a fit that never ran.
///
/// `None` when every point shares one coordinate (or there are no points at
/// all): a `CalibrationState` cannot be configured with a zero-width range.
fn point_ranges(points: &[CalibrantPoint]) -> Option<((f64, f64), (f64, f64))> {
    let (min_x, max_x, min_y, max_y) = points.iter().fold(
        (
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
        ),
        |(mnx, mxx, mny, mxy), p| {
            (
                mnx.min(p.library_rt),
                mxx.max(p.library_rt),
                mny.min(p.observed_rt),
                mxy.max(p.observed_rt),
            )
        },
    );
    let usable = min_x.is_finite() && max_x.is_finite() && min_y.is_finite() && max_y.is_finite();
    (usable && min_x < max_x && min_y < max_y).then_some(((min_x, max_x), (min_y, max_y)))
}

/// `CalibrantPoint` as `calibrt::CalibrationState::update`/`CalibrationCurve::wrmse`
/// want it.
fn as_calibrt_tuple(p: &CalibrantPoint) -> (LibraryRT<f64>, ObservedRTSeconds<f64>, f64) {
    (
        LibraryRT(p.library_rt),
        ObservedRTSeconds(p.observed_rt),
        p.score,
    )
}

/// Summary the Tolerances tab renders once Step B has run: the derived m/z,
/// mobility and RT tolerances plus how many calibrants they were estimated
/// from. `CalibDash::show_final` receives one of these alongside the
/// post-Phase-2 `FitRecording`, and `show_final` sets it on `App` (not
/// `CalibDash`) since `ui.rs`'s `draw_tolerances_tab` only ever sees an
/// `App`.
///
/// `mz_ppm` and `mobility_pct` are genuinely asymmetric `(lo, hi)` pairs —
/// Step B can derive a different tolerance on either side of a calibrant.
/// `rt_seconds` is not: `timsseek_cli::processing::calibrate_from_phase1`
/// only ever derives one scalar RT tolerance and duplicates it into both
/// tuple slots so this field's shape matches the other two, not because the
/// two sides differ — see that field's own doc comment.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ToleranceSummary {
    pub mz_ppm: (f64, f64),
    pub mobility_pct: (f64, f64),
    /// Symmetric by construction: both elements are always equal (see the
    /// struct doc comment). Not `f64` alone because the Tolerances tab
    /// renders all three tolerances the same `(lo, hi)` way; a lone scalar
    /// here would need its own special-cased rendering for no real benefit.
    pub rt_seconds: (f64, f64),
    pub n_calibrants: usize,
}

/// Ties the dashboard to a running Phase 1 batch loop.
///
/// `on_batch` is called once per scoring chunk: it records the chunk's
/// calibrant snapshot into the frame slab, re-fits the live calibration
/// curve from exactly those points, computes and records this batch's
/// convergence metrics, and either lets the batch loop continue unattended
/// (`Flow::Continue`) or pauses to render and block on a keypress.
///
/// Steady state allocates nothing new per batch beyond what `calibrt` itself
/// already allocates for every fit (`fit_with` builds a fresh
/// `CalibrationCurve` each call, regardless of caller): `state` and
/// `refit_state` are each a single `CalibrationState` reused via
/// `reconfigure` (`bins` never changes, only the point-derived ranges do),
/// `app`'s live `FitRecording` and `refit_recording` are reused the same way
/// `fit_with` always reuses a `FitRecording`, and `frames`' slab is allocated
/// once at construction. `current_points` is pre-sized to `n_calibrants` at
/// construction and reused via `clear`+`extend` from the very first batch;
/// `prev_points` starts as an empty `Vec` and grows once, on the first
/// `on_batch` call, to that same steady-state size — a one-time asymmetry
/// between the two, not a per-batch cost. `prev_curve` does clone one
/// `CalibrationCurve` per batch to keep a comparison point for `curve_delta`
/// — an additional cost of the same kind `calibrt` already pays per fit, not
/// a new category of allocation this dashboard introduces.
pub struct CalibDash {
    app: App,
    frames: FrameStore,
    stepper: Stepper,
    /// The live re-fit's `CalibrationState`, reused every batch.
    state: CalibrationState,
    /// This batch's points, buffered once so recording, re-fitting and churn
    /// all read the same materialized slice rather than each consuming a
    /// fresh copy of the caller's iterator (which is only handed to
    /// `on_batch` once). Clamped to `n_calibrants` exactly the way
    /// `FrameStore::record` clamps its own copy (`.take(self.n_calibrants)`)
    /// — without matching clamps, an overlong batch would fit the live curve
    /// against every point while the frame slab silently recorded only the
    /// first `n_calibrants` of them, so a later replay of that frame would
    /// disagree with what was actually shown live.
    current_points: Vec<CalibrantPoint>,
    n_calibrants: usize,
    /// The previous batch's points, for `churn`. Empty before the first
    /// batch — `churn`'s own "everything admitted, nothing evicted" case.
    prev_points: Vec<CalibrantPoint>,
    /// The previous batch's fitted curve, for `curve_delta`. `None` before
    /// the first successful fit; that batch's `max_delta`/`mean_delta` are
    /// reported as `NaN` (nothing to compare against yet), the same
    /// "nothing to average" convention `weighted_ridge_half_width` uses.
    prev_curve: Option<CalibrationCurve>,
    /// `bins`, held so `refit_frame` can reconfigure `refit_state` to the
    /// same grid resolution the live fit used without `CalibDash` having to
    /// go back through `app.recording().geom()` for it.
    bins: usize,
    /// A wholly separate `CalibrationState`/`FitRecording` pair for replaying
    /// a recorded frame on demand (scrubbing), so stepping through history
    /// never perturbs `state`/`app`'s live recording — what is currently on
    /// screen for the *live* batch stays exactly as it was fit.
    refit_state: CalibrationState,
    refit_recording: FitRecording,
    /// Counts calls into `render_pause` — test-only instrumentation so
    /// `present`'s stepper guard (see its doc comment) is provable without a
    /// real terminal: under `cargo test`, `render_pause` always detaches
    /// immediately regardless of whether the guard let it through, so a
    /// test asserting `present()` "did nothing" can't tell the two apart by
    /// output alone. This is what it checks instead.
    #[cfg(test)]
    render_pause_calls: u32,
}

impl CalibDash {
    pub fn new(
        n_frames: usize,
        n_calibrants: usize,
        bins: usize,
        lookback: usize,
        budget_bytes: usize,
    ) -> Self {
        // The unit range is only ever a placeholder: the first `on_batch`
        // call reconfigures both states to the batch's own point-derived
        // ranges before anything is fit against them. It is fixed and always
        // nonzero-width, so the only way `CalibrationState::new` can fail
        // here is `bins == 0` — the `.expect` messages below name that, not
        // the placeholder range, so a caller who does pass `bins == 0`
        // (there is no `Result` on this constructor to report it through
        // otherwise — see the standalone binary's `validate_snapshot`, which
        // exists specifically to catch this before it ever reaches here)
        // gets a message that points at the actual cause.
        let placeholder = (0.0, 1.0);
        Self {
            app: App::new(bins),
            frames: FrameStore::new(n_frames, n_calibrants, budget_bytes),
            stepper: Stepper::new(),
            state: CalibrationState::new(bins, placeholder, placeholder, lookback)
                .expect("bins must be nonzero (the placeholder range is always valid)"),
            current_points: Vec::with_capacity(n_calibrants),
            n_calibrants,
            prev_points: Vec::new(),
            prev_curve: None,
            bins,
            refit_state: CalibrationState::new(bins, placeholder, placeholder, lookback)
                .expect("bins must be nonzero (the placeholder range is always valid)"),
            refit_recording: FitRecording::new(bins),
            #[cfg(test)]
            render_pause_calls: 0,
        }
    }

    /// Behavior, in order: record the frame, re-fit the live curve from
    /// those same points, compute and push this batch's metrics (every
    /// batch, rendered or not — see `metrics.rs`'s module doc on why the
    /// series must have no holes), then either skip straight back to
    /// `Flow::Continue` or pause to render.
    ///
    /// Never returns `Err`: a dashboard that cannot fit or cannot draw must
    /// not stop a search. `Flow::Abort` is the one path back to the caller
    /// that actually asks the search to stop, and it only happens because
    /// the user pressed Ctrl-C at a pause — see `render_pause`.
    pub fn on_batch(
        &mut self,
        chunk: usize,
        range: Range<usize>,
        points: impl Iterator<Item = CalibrantPoint>,
    ) -> Flow {
        self.current_points.clear();
        self.current_points.extend(points.take(self.n_calibrants));

        self.frames
            .record(chunk, range, self.current_points.iter().copied());
        self.sync_frame_summary();

        let metrics = self.refit_live(chunk);
        self.app.push_metrics(metrics);
        self.app.batch = chunk as u32;

        if !self.stepper.should_render() {
            return Flow::Continue;
        }

        let action = render_pause(self);
        self.stepper.apply(action);
        match action {
            PauseAction::Abort => Flow::Abort,
            _ => Flow::Continue,
        }
    }

    /// Opens the post-Phase-2 pause: the same `render_pause` `on_batch` uses
    /// for every Phase 1 pause, so the Tolerances tab (and the rest of the
    /// dashboard) is reachable exactly once more after Step B runs — see the
    /// module's "two moments" design. There is no batch loop left to steer
    /// at this point, so every `PauseAction` this can return collapses to
    /// the same thing: return to the caller. `Abort` in particular is a
    /// no-op here rather than propagated — `processing.rs` calls this right
    /// after `show_final`, with nothing left downstream for an abort to
    /// cancel.
    ///
    /// Guarded on `stepper.is_stopped()`: a user who chose `r`/`q`/Ctrl-C at
    /// any Phase 1 pause has already asked to stop seeing the dashboard (or
    /// stop the search outright) — `on_batch` already honors that via
    /// `Stepper::should_render`, and this must too, or that choice stops
    /// meaning anything the moment Phase 2 finishes: `processing.rs` reaches
    /// this call regardless of `Flow::Abort` (Phase 1 breaking its own loop
    /// early does not skip Phase 2), so an unguarded `present` would reopen
    /// the alternate screen and block in `event::read()` for a keypress the
    /// user already told it not to wait for — stalling the run rather than
    /// failing it, but defeating the one documented way to end the search
    /// promptly either way.
    pub fn present(&mut self) {
        if self.stepper.is_stopped() {
            return;
        }
        render_pause(self);
    }

    /// Promotes the reserved final frame (if `finish` hasn't already run for
    /// it) and refreshes the Convergence header's decimation numbers one
    /// last time — `finish` can change `frames.len()` even when no further
    /// batch will call `sync_frame_summary` again.
    pub fn finish(&mut self, last_chunk: usize) {
        self.frames.finish(last_chunk);
        self.sync_frame_summary();
    }

    /// Wires the post-Phase-2 recording and derived tolerances into the
    /// dashboard once Step B has run and `calibration.json` is written. Both
    /// land on `App` (not `CalibDash`) since `ui.rs`'s `draw_tolerances_tab`
    /// only ever receives a plain `&App`. The Tolerances tab reads `real_fit`
    /// as its signal that Step B exists at all — see `App::real_fit`'s doc
    /// comment.
    pub fn show_final(&mut self, recording: FitRecording, tolerances: ToleranceSummary) {
        self.app.set_final(recording, tolerances);
    }

    /// Re-fits a recorded frame's points from scratch, into `refit_state`/
    /// `refit_recording` rather than the live `state`/`app.recording`, using
    /// exactly `point_ranges` — the same fold `refit_live` runs — so what
    /// this reproduces is the batch that actually ran, not an approximation
    /// of it. `None` if the frame index doesn't exist or its points are too
    /// degenerate to configure a grid from (mirrors `on_batch`'s own
    /// skip-the-refit behavior for the live path).
    pub fn refit_frame(&mut self, i: usize) -> Option<&FitRecording> {
        let (_idx, pts) = self.frames.frame(i)?;
        let (x_range, y_range) = point_ranges(pts)?;
        self.refit_state
            .reconfigure(self.bins, x_range, y_range)
            .ok()?;
        self.refit_state
            .update(pts.iter().map(as_calibrt_tuple))
            .ok()?;
        self.refit_state
            .fit_with(&mut self.refit_recording, ObserveOpts { dp_nodes: true });
        self.refit_state
            .measure_ridge_width_with(RIDGE_FRACTION, &mut self.refit_recording);
        Some(&self.refit_recording)
    }

    pub fn metrics(&self) -> &[BatchMetrics] {
        self.app.metrics()
    }

    pub fn frames(&self) -> &FrameStore {
        &self.frames
    }

    /// The live batch's recording — the one `on_batch` just fit into,
    /// distinct from whatever `refit_frame` most recently produced.
    pub fn recording(&self) -> &FitRecording {
        self.app.recording()
    }

    /// Only pins the grid-weights buffer, not the per-DP-node `considered`
    /// `Vec`s `refit_frame`'s `ObserveOpts { dp_nodes: true }` fills:
    /// `FitRecording::on_event` pushes a fresh `considered.to_vec()` into
    /// `self.dp` for every `DpNode` event, and `reset_body_only`'s
    /// `self.dp.clear()` drops those per-node `Vec`s along with the elements
    /// — reusing them would mean `calibrt`'s own `FitRecording` tracking a
    /// pool of them, which is out of this task's scope. See
    /// `a_repeated_refit_does_not_reallocate`'s doc comment.
    #[cfg(test)]
    pub fn refit_capacity(&self) -> usize {
        self.refit_recording.weights_capacity()
    }

    #[cfg(test)]
    pub fn current_points_ptr(&self) -> *const CalibrantPoint {
        self.current_points.as_ptr()
    }

    #[cfg(test)]
    pub fn prev_points_ptr(&self) -> *const CalibrantPoint {
        self.prev_points.as_ptr()
    }

    /// How many times `render_pause` has actually been entered — see the
    /// field's own doc comment for why a test needs this rather than just
    /// reading `present`'s return value (there isn't one) or the terminal
    /// state (there is no terminal under `cargo test`).
    #[cfg(test)]
    pub fn render_pause_calls(&self) -> u32 {
        self.render_pause_calls
    }

    /// Re-fits the live `CalibrationState` from `current_points` and returns
    /// this batch's `BatchMetrics`. See the "Behavior of `on_batch`" list:
    /// step 2 (re-fit) and step 3 (metrics) both live here since neither
    /// makes sense without the other — the metrics *are* a comparison of
    /// this re-fit against the last one.
    fn refit_live(&mut self, chunk: usize) -> BatchMetrics {
        let n_points = self.current_points.len();
        let mut wrmse = f64::NAN;
        let mut max_delta = f64::NAN;
        let mut mean_delta = f64::NAN;
        let mut path_nodes = 0usize;
        let mut ridge_half_width = f64::NAN;

        match point_ranges(&self.current_points) {
            None => {
                tracing::warn!(
                    chunk,
                    "calib_dash: this batch's points span a zero-width range on at least \
                     one axis; skipping the re-fit (metrics for this batch report n_points \
                     only)"
                );
            }
            Some((x_range, y_range)) => {
                if let Err(e) = self.state.reconfigure(self.bins, x_range, y_range) {
                    tracing::warn!(
                        chunk,
                        error = ?e,
                        "calib_dash: failed to reconfigure the live calibration state; \
                         skipping this batch's re-fit"
                    );
                } else if let Err(e) = self
                    .state
                    .update(self.current_points.iter().map(as_calibrt_tuple))
                {
                    tracing::warn!(
                        chunk,
                        error = ?e,
                        "calib_dash: failed to update the live calibration state; skipping \
                         this batch's re-fit"
                    );
                } else {
                    self.state
                        .fit_with(self.app.recording_mut(), ObserveOpts::NONE);
                    self.state
                        .measure_ridge_width_with(RIDGE_FRACTION, self.app.recording_mut());
                    path_nodes = self.app.recording().path_indices().len();
                    ridge_half_width = weighted_ridge_half_width(self.app.recording().ridge());

                    if let Some(curve) = self.state.curve() {
                        wrmse = curve.wrmse(self.current_points.iter().map(as_calibrt_tuple));
                        if let Some(prev) = &self.prev_curve {
                            let (d_max, d_mean) =
                                curve_delta(prev, curve, x_range, CURVE_DELTA_SAMPLES);
                            max_delta = d_max;
                            mean_delta = d_mean;
                        }
                        // Only overwrite the baseline when this batch actually
                        // produced a curve: a batch whose fit fails (e.g.
                        // `suppress_nonmax` finds nothing survives) must not
                        // reset `prev_curve` to `None`, or the *next*
                        // successful batch would report a NaN delta as if
                        // there were no prior curve at all, instead of
                        // comparing against the last real one — a
                        // discontinuity on the Convergence tab that never
                        // actually happened.
                        self.prev_curve = Some(curve.clone());
                    }
                }
            }
        }

        let (admitted, evicted) = churn(&self.prev_points, &self.current_points);
        self.prev_points.clear();
        self.prev_points.extend(self.current_points.iter().copied());

        BatchMetrics {
            chunk,
            n_points,
            wrmse,
            max_delta,
            mean_delta,
            path_nodes,
            ridge_half_width,
            admitted,
            evicted,
        }
    }

    /// Pushes `frames`' real decimation numbers into `app` — the
    /// `retained_frames`/`frame_stride`/`dropped_frames` reconciliation: one
    /// real source (`FrameStore`), fed through `App::set_frame_summary`,
    /// rather than a second tracker duplicating what `FrameStore` already
    /// knows.
    fn sync_frame_summary(&mut self) {
        self.app.set_frame_summary(
            self.frames.len(),
            self.frames.keep_every(),
            self.frames.dropped(),
        );
    }

    /// Keeps the Fit tab's scrub view in lockstep with `App::scrub_frame`.
    /// Called at the top of every `event_loop` iteration, before drawing —
    /// `App` moves `scrub_frame` on `<`/`>` but has no `FrameStore` to refit
    /// from, so this is the other half of that key handling, run once the
    /// dashboard (which does have one) is back in the loop.
    ///
    /// A no-op when `scrub_frame` is `None` (the live view — the common
    /// case, so this costs nothing then). Re-fits and clones unconditionally
    /// otherwise rather than caching whether the index already matches what
    /// was last shown: this runs at most once per keypress inside an
    /// interactive pause, never per batch, so the clone `refit_frame`'s
    /// result needs to hand `App` an owned copy is not on any path
    /// `on_batch`'s steady state cares about.
    fn sync_scrub(&mut self) {
        let Some(i) = self.app.scrub_frame() else {
            return;
        };
        let chunk = self.frames.frame(i).map(|(idx, _)| idx.chunk);
        match self.refit_frame(i) {
            Some(rec) => {
                let cloned = rec.clone();
                self.app.set_scrub_recording(i, chunk, cloned);
            }
            // The index should always be valid — it is bounded by
            // `retained_frames`, set from this same `FrameStore` — but if it
            // ever isn't (or the frame's points are too degenerate to
            // refit), fall back to the live view rather than leave a stale
            // or missing recording on screen.
            None => self.app.clear_scrub(),
        }
    }
}

/// Pauses the batch loop to render one interactive frame and block until the
/// user's keypress resolves to something other than [`PauseAction::Stay`].
/// Takes the whole [`CalibDash`], not just its `App` — the batch scrubber
/// (`<`/`>`) needs `refit_frame`'s `FrameStore` access every time it redraws,
/// which only `CalibDash` has.
///
/// Never fails the caller: when stdout is not a terminal (e.g. under `cargo
/// test`, or a piped/redirected run) this returns [`PauseAction::Detach`]
/// immediately, without touching the terminal at all — the render path must
/// short-circuit rather than block on a keypress that will never arrive. A
/// terminal that *is* present but fails to initialize (`ratatui::try_init`)
/// is logged at `warn` and also treated as `Detach`: a dashboard that cannot
/// draw must not stop a search.
fn render_pause(dash: &mut CalibDash) -> PauseAction {
    #[cfg(test)]
    {
        dash.render_pause_calls += 1;
    }
    if !std::io::stdout().is_terminal() {
        return PauseAction::Detach;
    }

    let mut terminal = match ratatui::try_init() {
        Ok(terminal) => terminal,
        Err(e) => {
            // `try_init` can fail after partially succeeding (raw mode
            // enabled, then entering the alternate screen fails), so restore
            // unconditionally rather than assuming there is nothing to undo.
            tracing::warn!("calib_dash failed to initialize the terminal: {e}");
            ratatui::restore();
            return PauseAction::Detach;
        }
    };
    let action = catch_panics(|| event_loop(&mut terminal, dash));
    ratatui::restore();
    action
}

/// Runs `f`, converting a panic into `PauseAction::Detach` rather than
/// letting it unwind into the search — in builds where panics unwind at all.
/// This workspace's `[profile.release]` sets `panic = "abort"`, so in the
/// shipped release binary a panic inside `f` aborts the process immediately;
/// `catch_unwind` never gets a chance to run and this recovery does not
/// apply there. It is real in dev/test builds (`panic = "unwind"`, the
/// default `profile.dev`), which is what this module's tests exercise, and
/// is kept anyway because unwinding is strictly better than nothing when it
/// is available. Factored out of `render_pause` so the recovery behavior is
/// testable without a real terminal — `event_loop` itself needs one, this
/// wrapper does not.
fn catch_panics(f: impl FnOnce() -> PauseAction) -> PauseAction {
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(f)) {
        Ok(action) => action,
        Err(_) => {
            tracing::warn!(
                "calib_dash dashboard panicked while rendering a pause; detaching so the \
                 search is not aborted"
            );
            PauseAction::Detach
        }
    }
}

/// Draws and reads key events in the alternate screen until a keypress
/// resolves to something other than `Stay`. A draw or read failure is
/// treated the same as `render_pause`'s terminal-setup failure: logged and
/// detached, never propagated to fail the batch loop.
fn event_loop<B: ratatui::backend::Backend>(
    terminal: &mut ratatui::Terminal<B>,
    dash: &mut CalibDash,
) -> PauseAction {
    loop {
        // Before every draw, not just after a scrub keypress: simpler than
        // tracking whether `scrub_frame` actually changed since the last
        // iteration, and a no-op whenever it's `None` (the common, live
        // case) — see `sync_scrub`'s own doc comment.
        dash.sync_scrub();
        if let Err(e) = terminal.draw(|f| crate::ui::draw(f, &mut dash.app)) {
            tracing::warn!("calib_dash failed to draw a frame: {e}");
            return PauseAction::Detach;
        }
        let ev = match event::read() {
            Ok(ev) => ev,
            Err(e) => {
                tracing::warn!("calib_dash failed to read a terminal event: {e}");
                return PauseAction::Detach;
            }
        };
        // Only key *presses*: on Windows crossterm also emits releases,
        // which would double every keystroke.
        if let Event::Key(key) = ev
            && key.kind == KeyEventKind::Press
        {
            let action = dash.app.handle_key(key);
            if !matches!(action, PauseAction::Stay) {
                return action;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ratatui::crossterm::event::{
        KeyCode,
        KeyEvent,
        KeyModifiers,
    };

    fn key(c: char) -> KeyEvent {
        KeyEvent::new(KeyCode::Char(c), KeyModifiers::NONE)
    }

    fn ctrl(c: char) -> KeyEvent {
        KeyEvent::new(KeyCode::Char(c), KeyModifiers::CONTROL)
    }

    fn press(app: &mut App, s: &str) -> PauseAction {
        let mut last = PauseAction::Stay;
        for c in s.chars() {
            last = app.handle_key(key(c));
        }
        last
    }

    // ---- count prefix ----

    #[test]
    fn digits_accumulate_into_a_pending_count() {
        let mut app = App::new(10);
        app.handle_key(key('1'));
        app.handle_key(key('5'));
        assert_eq!(app.pending_count(), Some(15));
    }

    #[test]
    fn a_motion_consumes_the_count() {
        let mut app = App::new(10);
        assert_eq!(press(&mut app, "15n"), PauseAction::Next(15));
        assert_eq!(app.pending_count(), None, "count lives exactly one motion");
    }

    #[test]
    fn a_bare_motion_means_one() {
        let mut app = App::new(10);
        assert_eq!(press(&mut app, "n"), PauseAction::Next(1));
    }

    #[test]
    fn a_leading_zero_does_not_start_a_count() {
        let mut app = App::new(10);
        app.handle_key(key('0'));
        assert_eq!(app.pending_count(), None);
    }

    #[test]
    fn zero_continues_an_existing_count() {
        let mut app = App::new(10);
        press(&mut app, "10");
        assert_eq!(app.pending_count(), Some(10));
    }

    #[test]
    fn a_long_digit_run_saturates_instead_of_overflowing() {
        let mut app = App::new(10);
        // 16 nines: `9_999_999_999 * 10 + 9` overflows a u32 well before the
        // end of this run, so this both must not panic (overflow checks are
        // on in dev/test builds) and must land on a sane clamped value.
        press(&mut app, "9999999999999999");
        assert_eq!(app.pending_count(), Some(u32::MAX));
    }

    #[test]
    fn esc_clears_a_pending_count() {
        let mut app = App::new(10);
        press(&mut app, "42");
        app.handle_key(KeyEvent::new(KeyCode::Esc, KeyModifiers::NONE));
        assert_eq!(app.pending_count(), None);
    }

    #[test]
    fn a_non_motion_key_clears_a_pending_count() {
        let mut app = App::new(10);
        press(&mut app, "42d");
        assert_eq!(
            app.pending_count(),
            None,
            "`d` toggles the DP pane, not a motion"
        );
    }

    #[test]
    fn counts_drive_list_and_scrub_motions_too() {
        let mut app = App::new(10);
        press(&mut app, "3]");
        assert_eq!(app.stage(), Stage::Curve, "Grid + 3 stages");
        assert_eq!(app.pending_count(), None);
    }

    #[test]
    fn a_huge_count_on_a_stage_motion_is_capped_by_the_state_space() {
        // `Stage` has 5 variants, so no count can usefully move it more than
        // 4 steps. Without the cap this loops once per requested step —
        // typed as a huge count, that would freeze the UI for as long as the
        // loop takes, for work none of which has any effect past step 4.
        let mut app = App::new(10);
        press(&mut app, "999999999]");
        assert_eq!(
            app.stage(),
            Stage::Ridge,
            "clamped to the last stage, not looped hundreds of millions of times"
        );
    }

    // ---- pause actions ----

    #[test]
    fn r_runs_to_end_and_q_detaches() {
        let mut app = App::new(10);
        assert_eq!(press(&mut app, "r"), PauseAction::RunToEnd);
        assert_eq!(press(&mut app, "q"), PauseAction::Detach);
    }

    #[test]
    fn ctrl_c_aborts() {
        let mut app = App::new(10);
        assert_eq!(app.handle_key(ctrl('c')), PauseAction::Abort);
    }

    #[test]
    fn ctrl_c_is_not_a_stage_or_overlay_key() {
        let mut app = App::new(10);
        let before = app.stage();
        app.handle_key(ctrl(']'));
        assert_eq!(app.stage(), before, "modified keys are not motions");
    }

    // ---- tabs: h/l ----

    #[test]
    fn l_cycles_forward_through_tab_all() {
        let mut app = App::new(10);
        assert_eq!(app.tab(), Tab::Fit);
        press(&mut app, "l");
        assert_eq!(app.tab(), Tab::Convergence);
        press(&mut app, "l");
        assert_eq!(app.tab(), Tab::Tolerances);
        press(&mut app, "l");
        assert_eq!(app.tab(), Tab::Fit, "wraps back to the first tab");
    }

    #[test]
    fn h_cycles_backward_and_also_wraps() {
        let mut app = App::new(10);
        press(&mut app, "h");
        assert_eq!(
            app.tab(),
            Tab::Tolerances,
            "one step back from the first tab wraps to the last"
        );
    }

    #[test]
    fn a_count_before_l_cycles_that_many_tabs() {
        let mut app = App::new(10);
        press(&mut app, "2l");
        assert_eq!(app.tab(), Tab::Tolerances, "Fit + 2 steps");
        assert_eq!(
            app.pending_count(),
            None,
            "l is a motion: it consumes the count"
        );
    }

    #[test]
    fn a_huge_count_before_l_does_not_spin() {
        // 3 tabs: any count's effect only depends on its value mod 3, and
        // `tab_cycle_steps` reduces it to that before looping at all, so this
        // must return promptly rather than looping ~10^9 times.
        // `1_000_000_001 % 3 == 2`, so this is Fit + 2 steps == Tolerances.
        let mut app = App::new(10);
        press(&mut app, "1000000001l");
        assert_eq!(
            app.tab(),
            Tab::Tolerances,
            "Fit + (1_000_000_001 % 3) steps"
        );
    }

    // ---- Fit-tab batch scrubber: < / > ----

    #[test]
    fn scrub_is_a_no_op_with_no_retained_frames() {
        let mut app = App::new(10);
        press(&mut app, "<");
        assert_eq!(app.scrub_frame(), None, "nothing to scrub back to");
        press(&mut app, ">");
        assert_eq!(app.scrub_frame(), None);
    }

    #[test]
    fn less_than_steps_back_from_live_into_the_most_recent_retained_frame() {
        let mut app = App::new(10);
        app.set_frame_summary(5, 1, 0);
        press(&mut app, "<");
        assert_eq!(
            app.scrub_frame(),
            Some(4),
            "one step back from live lands on the last retained frame (index 4 of 5)"
        );
    }

    #[test]
    fn less_than_with_a_count_steps_back_that_many_frames() {
        let mut app = App::new(10);
        app.set_frame_summary(5, 1, 0);
        press(&mut app, "3<");
        assert_eq!(app.scrub_frame(), Some(2));
        assert_eq!(
            app.pending_count(),
            None,
            "< is a motion: it consumes the count"
        );
    }

    #[test]
    fn less_than_clamps_at_the_first_frame_rather_than_underflowing() {
        let mut app = App::new(10);
        app.set_frame_summary(5, 1, 0);
        press(&mut app, "999999999<");
        assert_eq!(
            app.scrub_frame(),
            Some(0),
            "a huge count clamps at the first frame, not a wrapped/underflowed index"
        );
    }

    #[test]
    fn greater_than_steps_forward_and_returns_to_live_past_the_last_frame() {
        let mut app = App::new(10);
        app.set_frame_summary(5, 1, 0);
        press(&mut app, "<"); // frame 4 (of 5)
        assert_eq!(app.scrub_frame(), Some(4));
        press(&mut app, ">");
        assert_eq!(
            app.scrub_frame(),
            None,
            "> past the last retained frame returns to the live view"
        );
    }

    #[test]
    fn greater_than_is_a_no_op_while_already_live() {
        let mut app = App::new(10);
        app.set_frame_summary(5, 1, 0);
        press(&mut app, ">");
        assert_eq!(app.scrub_frame(), None);
    }

    #[test]
    fn a_huge_count_before_greater_than_does_not_spin() {
        // Saturating arithmetic, not a loop — this must return promptly
        // however large the count is.
        let mut app = App::new(10);
        app.set_frame_summary(5, 1, 0);
        press(&mut app, "<"); // frame 4
        press(&mut app, "999999999>");
        assert_eq!(
            app.scrub_frame(),
            None,
            "saturates past the end, back to live"
        );
    }

    /// `n` advances the *live* batch loop, so a scrub cursor left over from
    /// browsing history at this pause must not silently carry into the next
    /// one — the bound stays safe and the banner would still name the real
    /// batch either way, but reopening on a replayed frame right after the
    /// user asked to advance reads as confusing. `n` returning to the live
    /// view is the more useful default.
    #[test]
    fn n_clears_a_scrub_cursor_left_over_from_this_pause() {
        let mut app = App::new(10);
        app.set_frame_summary(5, 1, 0);
        press(&mut app, "<");
        assert_eq!(
            app.scrub_frame(),
            Some(4),
            "sanity: scrubbed back one frame"
        );
        assert_eq!(press(&mut app, "n"), PauseAction::Next(1));
        assert_eq!(
            app.scrub_frame(),
            None,
            "n must return to the live view, not keep replaying frame 4"
        );
    }

    // ---- Fit-tab overlay toggles: s / p / c / w ----

    #[test]
    fn s_p_c_w_toggle_their_own_overlay_independently() {
        let mut app = App::new(10);
        assert!(
            !app.show_suppressed() && !app.show_path() && !app.show_curve() && !app.show_ridge()
        );
        press(&mut app, "s");
        assert!(app.show_suppressed());
        assert!(!app.show_path(), "s must not touch the other three");
        press(&mut app, "p");
        assert!(app.show_path());
        press(&mut app, "c");
        assert!(app.show_curve());
        press(&mut app, "w");
        assert!(app.show_ridge());
        assert!(
            app.show_suppressed() && app.show_path() && app.show_curve(),
            "each toggle key must leave the others exactly as they were"
        );
    }

    #[test]
    fn an_overlay_toggle_consumes_the_pending_count_by_parity() {
        // Toggling a boolean twice is a no-op — the count's only observable
        // effect is odd vs. even, and `toggle_parity` reduces to that before
        // doing any work at all, so a huge count is not a huge loop either.
        let mut app = App::new(10);
        press(&mut app, "2s");
        assert!(!app.show_suppressed(), "an even count leaves it unchanged");
        assert_eq!(
            app.pending_count(),
            None,
            "s is a motion: it consumes the count"
        );

        press(&mut app, "3s");
        assert!(app.show_suppressed(), "an odd count flips it exactly once");

        // 999_999_999 is odd, so this flips it once more, back to false.
        press(&mut app, "999999999s");
        assert!(!app.show_suppressed());
    }

    // ---- the stepper ----

    #[test]
    fn a_fresh_stepper_renders_every_batch() {
        let mut s = Stepper::new();
        assert!(s.should_render());
        assert!(s.should_render());
    }

    #[test]
    fn next_n_renders_on_the_nth_batch_not_the_n_plus_first() {
        let mut s = Stepper::new();
        // At the chunk-0 pause the user asks for 5.
        s.apply(PauseAction::Next(5));
        // Chunks 1..=4 are skipped; the 5th call renders.
        let rendered: Vec<usize> = (1..=10).filter(|_| s.should_render()).collect();
        assert_eq!(
            rendered.first(),
            Some(&5),
            "chunks 1-4 skipped, render at 5"
        );
    }

    #[test]
    fn skipping_is_consumed_once_rendered() {
        let mut s = Stepper::new();
        s.apply(PauseAction::Next(3));
        assert!(!s.should_render());
        assert!(!s.should_render());
        assert!(s.should_render());
        assert!(
            s.should_render(),
            "back to every batch after the skip is spent"
        );
    }

    #[test]
    fn run_to_end_never_renders_again() {
        let mut s = Stepper::new();
        s.apply(PauseAction::RunToEnd);
        assert!((0..1000).all(|_| !s.should_render()));
    }

    #[test]
    fn detach_never_renders_again() {
        let mut s = Stepper::new();
        s.apply(PauseAction::Detach);
        assert!((0..1000).all(|_| !s.should_render()));
    }

    /// `is_stopped` is what `CalibDash::present` consults instead of
    /// `should_render` — read-only, so checking it must not itself draw down
    /// a pending skip count the way `should_render` deliberately does.
    #[test]
    fn is_stopped_reflects_run_to_end_detach_and_abort_without_side_effects() {
        for action in [
            PauseAction::RunToEnd,
            PauseAction::Detach,
            PauseAction::Abort,
        ] {
            let mut s = Stepper::new();
            assert!(!s.is_stopped(), "a fresh stepper has not stopped");
            s.apply(action);
            assert!(s.is_stopped(), "{action:?} must set stopped");
            // Checking it again must not change anything either.
            assert!(s.is_stopped());
        }
    }

    #[test]
    fn is_stopped_is_false_while_only_stay_or_next_has_been_applied() {
        let mut s = Stepper::new();
        s.apply(PauseAction::Stay);
        assert!(!s.is_stopped());
        s.apply(PauseAction::Next(5));
        assert!(!s.is_stopped(), "Next skips batches, it does not stop");
    }

    // ---- CalibDash::on_batch ----

    #[test]
    fn on_batch_without_a_terminal_never_blocks_and_never_aborts() {
        // stdout is not a tty under `cargo test`, so the render path must
        // short-circuit rather than wait for a keypress.
        let mut d = CalibDash::new(4, 8, 10, 5, 1 << 20);
        for chunk in 0..4 {
            let pts: Vec<_> = (0..8)
                .map(|i| CalibrantPoint {
                    library_rt: i as f64 + 0.5,
                    observed_rt: i as f64 + 0.5,
                    score: 1.0 + i as f64,
                    speclib_index: chunk * 8 + i,
                })
                .collect();
            assert!(matches!(
                d.on_batch(chunk, chunk..chunk + 1, pts.into_iter()),
                Flow::Continue
            ));
        }
        d.finish(3);
        assert_eq!(d.metrics().len(), 4, "a metric per batch, even unrendered");
    }

    #[test]
    fn metrics_are_collected_for_every_batch_including_skipped_ones() {
        let mut d = CalibDash::new(10, 8, 10, 5, 64); // tiny budget: heavy decimation
        for chunk in 0..10 {
            let pts: Vec<_> = (0..8)
                .map(|i| CalibrantPoint {
                    library_rt: i as f64 + 0.5,
                    observed_rt: i as f64 + 0.5,
                    score: 1.0 + i as f64,
                    speclib_index: i,
                })
                .collect();
            d.on_batch(chunk, chunk..chunk + 1, pts.into_iter());
        }
        d.finish(9);
        assert_eq!(d.metrics().len(), 10, "metrics are undecimated");
        assert!(d.frames().len() < 10, "frames are decimated");
    }

    #[test]
    fn refitting_a_recorded_frame_reproduces_the_live_curve() {
        // The app's central claim: what you scrub is what ran.
        let mut d = CalibDash::new(1, 8, 10, 5, 1 << 20);
        let pts: Vec<_> = (0..8)
            .map(|i| CalibrantPoint {
                library_rt: i as f64 + 0.5,
                observed_rt: (i as f64 + 0.5) * 1.5,
                score: 1.0 + i as f64,
                speclib_index: i,
            })
            .collect();
        d.on_batch(0, 0..1, pts.into_iter());
        d.finish(0);

        let live: Vec<_> = d.recording().curve().to_vec();
        let refit = d.refit_frame(0).expect("frame 0 exists");
        assert_eq!(refit.curve().len(), live.len());
        for (a, b) in live.iter().zip(refit.curve()) {
            assert!((a.library - b.library).abs() < 1e-12, "x differs");
            assert!((a.observed - b.observed).abs() < 1e-12, "y differs");
        }
    }

    /// End-to-end version of the `<` scrubber: moving `App::scrub_frame`
    /// with `handle_key` and then letting `CalibDash::sync_scrub` catch up
    /// must leave the Fit tab looking at the earlier batch's real curve, not
    /// the live one — `sync_scrub` is the piece that turns `scrub_frame` (an
    /// index `App` alone cannot refit from) into an actual recording
    /// `App::active_recording` can render.
    #[test]
    fn sync_scrub_reproduces_an_earlier_batchs_curve_on_the_fit_tab() {
        let mut d = CalibDash::new(3, 8, 10, 5, 1 << 20);
        for chunk in 0..3 {
            let pts: Vec<_> = (0..8)
                .map(|i| CalibrantPoint {
                    library_rt: i as f64 + 0.5,
                    observed_rt: (i as f64 + 0.5) * (1.0 + chunk as f64 * 0.1),
                    score: 1.0 + i as f64,
                    speclib_index: chunk * 8 + i,
                })
                .collect();
            d.on_batch(chunk, chunk..chunk + 1, pts.into_iter());
        }
        d.finish(2);

        // Batch 0's curve, refit independently for comparison against what
        // scrubbing shows.
        let batch0_curve: Vec<_> = d.refit_frame(0).expect("frame 0 exists").curve().to_vec();

        // A generous budget retains all 3 batches; scrub 3 back from live
        // (batch 2) lands on frame index 0 (batch 0).
        press(&mut d.app, "3<");
        d.sync_scrub();

        assert_eq!(d.app.scrub_frame(), Some(0));
        assert_eq!(d.app.scrub_chunk(), Some(0));
        let shown = d.app.active_recording().curve().to_vec();
        assert_eq!(shown.len(), batch0_curve.len());
        for (a, b) in shown.iter().zip(&batch0_curve) {
            assert!((a.library - b.library).abs() < 1e-12);
            assert!((a.observed - b.observed).abs() < 1e-12);
        }

        // The live view is unaffected by scrubbing — `recording()` (not
        // `active_recording()`) still reflects batch 2's own fit.
        assert!(
            !d.recording().curve().is_empty(),
            "sanity: batch 2 produced a real curve too"
        );
    }

    /// `present` reuses `render_pause`, which detaches immediately when
    /// stdout is not a terminal — always true under `cargo test` — so this
    /// must return promptly rather than blocking on a keypress that will
    /// never arrive. Mirrors
    /// `on_batch_without_a_terminal_never_blocks_and_never_aborts` for the
    /// post-Phase-2 pause.
    #[test]
    fn present_never_blocks_without_a_terminal() {
        let mut d = CalibDash::new(1, 4, 10, 3, 1 << 20);
        d.present();
    }

    /// `present` must still open the pause when nothing has stopped it —
    /// this is the baseline the next test's guard is checked against, so a
    /// bug that made the guard fire unconditionally wouldn't slip past.
    #[test]
    fn present_renders_when_the_stepper_has_not_stopped() {
        let mut d = CalibDash::new(1, 4, 10, 3, 1 << 20);
        d.present();
        assert_eq!(
            d.render_pause_calls(),
            1,
            "present() must still attempt render_pause when nothing stopped it"
        );
    }

    /// A user who pressed `q` (or `r`, or Ctrl-C) at a Phase 1 pause has
    /// already asked to stop seeing the dashboard — `processing.rs` reaches
    /// `present()` regardless (Phase 1 breaking its own loop early on
    /// `Flow::Abort` does not skip Phase 2), so without this guard the
    /// post-Phase-2 pause would reopen the alternate screen and block in
    /// `event::read()` for a keypress the user already declined to keep
    /// giving — the one documented escape hatch stops working the moment
    /// Phase 2 finishes. `render_pause_calls` is what makes "did not even
    /// attempt to render" provable here: under `cargo test`, `render_pause`
    /// always detaches immediately regardless of this guard (stdout isn't a
    /// terminal), so the *return value* can't tell a fixed bypass from a
    /// working guard — only whether it was entered at all can.
    #[test]
    fn present_does_not_render_once_a_pause_chose_detach() {
        let mut d = CalibDash::new(1, 4, 10, 3, 1 << 20);
        d.stepper.apply(PauseAction::Detach);
        d.present();
        assert_eq!(
            d.render_pause_calls(),
            0,
            "present() must not attempt render_pause once the stepper is stopped"
        );
    }

    /// Same guard, reached via `RunToEnd` (`r`) instead of `Detach` (`q`) —
    /// both, like `Abort`, set `Stepper::stopped`, and `present` must honor
    /// all three identically.
    #[test]
    fn present_does_not_render_once_a_pause_chose_run_to_end() {
        let mut d = CalibDash::new(1, 4, 10, 3, 1 << 20);
        d.stepper.apply(PauseAction::RunToEnd);
        d.present();
        assert_eq!(d.render_pause_calls(), 0);
    }

    /// Same guard, reached via `Abort` (Ctrl-C) — the case the review
    /// specifically called out: a user who explicitly asked to stop must
    /// not have the dashboard reopen on them after Phase 2 regardless.
    #[test]
    fn present_does_not_render_once_a_pause_chose_abort() {
        let mut d = CalibDash::new(1, 4, 10, 3, 1 << 20);
        d.stepper.apply(PauseAction::Abort);
        d.present();
        assert_eq!(d.render_pause_calls(), 0);
    }

    /// This is one of the four tests given verbatim in the task brief, so
    /// its name and body are kept as specified. What it actually pins is
    /// narrower than the name alone suggests: `refit_capacity` only reads
    /// the grid-weights buffer's capacity. `refit_frame` runs with
    /// `ObserveOpts { dp_nodes: true }`, and every `DpNode` event pushes a
    /// fresh `considered.to_vec()` (see `FitRecording::on_event`) that
    /// `reset_body_only`'s `self.dp.clear()` then drops along with the rest
    /// of `self.dp` on the next fit — those per-node allocations are real and
    /// this test does not (and cannot, without changing `calibrt`'s own
    /// `FitRecording` internals) claim they don't happen. See
    /// `refit_capacity`'s doc comment for the same scope note.
    #[test]
    fn a_repeated_refit_does_not_reallocate() {
        let mut d = CalibDash::new(2, 8, 10, 5, 1 << 20);
        for chunk in 0..2 {
            let pts: Vec<_> = (0..8)
                .map(|i| CalibrantPoint {
                    library_rt: i as f64 + 0.5,
                    observed_rt: i as f64 + 0.5,
                    score: 1.0 + i as f64,
                    speclib_index: i,
                })
                .collect();
            d.on_batch(chunk, chunk..chunk + 1, pts.into_iter());
        }
        d.finish(1);
        let cap = d.refit_capacity();
        for _ in 0..5 {
            d.refit_frame(0);
            d.refit_frame(1);
        }
        assert_eq!(d.refit_capacity(), cap);
    }

    /// An overlong batch (more points than `n_calibrants`) must not let the
    /// live re-fit see points the frame slab never recorded — `FrameStore`
    /// already truncates its own copy to `n_calibrants`
    /// (`an_overlong_frame_is_truncated_not_overrun`, frames.rs), and
    /// `on_batch` must clamp identically so a later `refit_frame` of this
    /// same batch reproduces what was actually fit live, not a curve fit
    /// against extra points the recording never kept.
    #[test]
    fn an_overlong_batch_is_clamped_the_same_way_the_frame_slab_clamps_it() {
        let n_calibrants = 4;
        let mut d = CalibDash::new(1, n_calibrants, 10, 3, 1 << 20);
        // Twice as many points as `n_calibrants`; the extra half must be
        // dropped identically by both the live fit and the recorded frame.
        let pts: Vec<_> = (0..8)
            .map(|i| CalibrantPoint {
                library_rt: i as f64 + 0.5,
                observed_rt: i as f64 + 0.5,
                score: 1.0 + i as f64,
                speclib_index: i,
            })
            .collect();
        d.on_batch(0, 0..8, pts.into_iter());
        d.finish(0);

        let (idx, _) = d.frames().frame(0).expect("frame 0 exists");
        assert_eq!(idx.len, n_calibrants, "the slab truncates to n_calibrants");

        let live: Vec<_> = d.recording().curve().to_vec();
        let refit = d.refit_frame(0).expect("frame 0 exists").curve().to_vec();
        assert_eq!(
            live.len(),
            refit.len(),
            "the live fit must have seen exactly the points the slab recorded"
        );
        assert_eq!(
            live.len(),
            n_calibrants,
            "not the full 8-point overlong batch"
        );
        for (a, b) in live.iter().zip(&refit) {
            assert!((a.library - b.library).abs() < 1e-12);
            assert!((a.observed - b.observed).abs() < 1e-12);
        }
    }

    // ---- terminal entry ----

    #[test]
    fn render_pause_detaches_when_stdout_is_not_a_terminal() {
        // Same guarantee `on_batch_without_a_terminal_never_blocks_and_never_aborts`
        // pins end-to-end, exercised directly against `render_pause` itself.
        // `render_pause` takes the whole `CalibDash` (not just `App`) so the
        // scrubber can reach `refit_frame` on every redraw.
        let mut d = CalibDash::new(1, 4, 10, 3, 1 << 20);
        assert_eq!(render_pause(&mut d), PauseAction::Detach);
    }

    /// `catch_panics` is what keeps a panic inside the event loop from
    /// unwinding into the search — in builds where panics unwind, which is
    /// what this test runs under (`panic = "unwind"`, the default dev/test
    /// profile). `[profile.release]` sets `panic = "abort"` workspace-wide,
    /// where this guard cannot run at all; see its doc comment. Exercised
    /// directly here since, unlike `event_loop`, it needs no real terminal.
    ///
    /// `#[cfg(debug_assertions)]`: this deliberately panics-and-catches,
    /// which only unwinds when `panic = "unwind"` (`profile.dev`/`test`'s
    /// default). A CI leg that runs tests against a release-profile build
    /// (`panic = "abort"`, set workspace-wide) would abort the whole test
    /// process on the `panic!` below instead of letting `catch_unwind`
    /// recover it — this guard is what keeps that leg green rather than
    /// weakening what the test proves under a normal dev/test build.
    #[cfg(debug_assertions)]
    #[test]
    fn catch_panics_converts_a_panic_into_detach() {
        let prev_hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(|_| {})); // silence the default panic printout
        let action = catch_panics(|| panic!("simulated event-loop panic"));
        std::panic::set_hook(prev_hook);
        assert_eq!(
            action,
            PauseAction::Detach,
            "a caught panic must not fail the batch loop"
        );
    }

    #[test]
    fn catch_panics_passes_through_a_successful_action() {
        assert_eq!(
            catch_panics(|| PauseAction::RunToEnd),
            PauseAction::RunToEnd
        );
    }

    // ---- no-reallocation, pinned by pointer identity ----
    //
    // Matches the convention used everywhere else in this codebase for a
    // "does not reallocate" claim (`FrameStore`'s `slab_capacity`/
    // `index_capacity`, `calibrt`'s `filtered_ptr`/`path_points_ptr`):
    // capacity equality alone can't distinguish "reused the same allocation"
    // from "reallocated but landed on the same capacity anyway" — pointer
    // identity is the part that actually proves reuse.

    #[test]
    fn current_points_never_reallocates_once_constructed() {
        let n_calibrants = 8;
        let mut d = CalibDash::new(5, n_calibrants, 10, 5, 1 << 20);
        // `current_points` is pre-sized to `n_calibrants` at construction —
        // unlike `prev_points` (see the next test), it should be stable from
        // before the very first batch, not just after it.
        let ptr_before_any_batch = d.current_points_ptr();
        for chunk in 0..5 {
            let pts: Vec<_> = (0..n_calibrants)
                .map(|i| CalibrantPoint {
                    library_rt: i as f64 + 0.5,
                    observed_rt: i as f64 + 0.5,
                    score: 1.0 + i as f64,
                    speclib_index: chunk * n_calibrants + i,
                })
                .collect();
            d.on_batch(chunk, chunk..chunk + 1, pts.into_iter());
            assert_eq!(
                d.current_points_ptr(),
                ptr_before_any_batch,
                "current_points must not reallocate on batch {chunk}"
            );
        }
    }

    #[test]
    fn prev_points_stabilizes_after_growing_once_on_the_first_batch() {
        let n_calibrants = 8;
        let mut d = CalibDash::new(5, n_calibrants, 10, 5, 1 << 20);
        let batch = |chunk: usize| -> Vec<CalibrantPoint> {
            (0..n_calibrants)
                .map(|i| CalibrantPoint {
                    library_rt: i as f64 + 0.5,
                    observed_rt: i as f64 + 0.5,
                    score: 1.0 + i as f64,
                    speclib_index: chunk * n_calibrants + i,
                })
                .collect()
        };
        // `prev_points` starts as an empty `Vec` (no allocation at all), so
        // its first-ever growth happens inside this very call — pointer
        // identity is only meaningful to compare *after* that first batch,
        // not before it.
        d.on_batch(0, 0..n_calibrants, batch(0).into_iter());
        let ptr_after_first_batch = d.prev_points_ptr();
        for chunk in 1..5 {
            d.on_batch(chunk, chunk..chunk + 1, batch(chunk).into_iter());
            assert_eq!(
                d.prev_points_ptr(),
                ptr_after_first_batch,
                "prev_points must not reallocate again after its first growth (batch {chunk})"
            );
        }
    }

    // ---- prev_curve survives a failed batch ----

    #[test]
    fn a_failed_fit_does_not_reset_the_delta_baseline() {
        // A batch whose fit fails (every weight lands below
        // `suppress_nonmax`'s 1.0 seed, so nothing survives suppression and
        // `CalibrationState::curve()` stays `None`) must not wipe out the
        // last real curve as the comparison point for the *next* successful
        // batch's `curve_delta` — otherwise the Convergence tab would show a
        // discontinuity (a sudden NaN-to-real jump) that never actually
        // happened.
        let n_calibrants = 4;
        let mut d = CalibDash::new(3, n_calibrants, 10, 3, 1 << 20);
        let good_points = || -> Vec<CalibrantPoint> {
            (0..n_calibrants)
                .map(|i| CalibrantPoint {
                    library_rt: i as f64 + 0.5,
                    observed_rt: i as f64 + 0.5,
                    score: 2.0 + i as f64, // well above the 1.0 suppression seed
                    speclib_index: i,
                })
                .collect()
        };
        let failing_points = || -> Vec<CalibrantPoint> {
            (0..n_calibrants)
                .map(|i| CalibrantPoint {
                    library_rt: i as f64 + 0.5,
                    observed_rt: i as f64 + 0.5,
                    score: 0.1, // below the seed: everything gets suppressed
                    speclib_index: i,
                })
                .collect()
        };

        d.on_batch(0, 0..n_calibrants, good_points().into_iter());
        assert!(
            d.metrics()[0].wrmse.is_finite(),
            "batch 0 must produce a real curve"
        );

        d.on_batch(1, 0..n_calibrants, failing_points().into_iter());
        assert!(
            d.metrics()[1].wrmse.is_nan(),
            "batch 1's fit must fail outright (every weight is sub-threshold)"
        );

        // Identical points to batch 0: if `prev_curve` correctly survived
        // batch 1's failure, this fits the same curve again and the delta
        // against it is ~0 — not NaN, which is what a reset-to-None baseline
        // would report instead.
        d.on_batch(2, 0..n_calibrants, good_points().into_iter());
        let m2 = d.metrics()[2];
        assert!(
            !m2.max_delta.is_nan() && !m2.mean_delta.is_nan(),
            "batch 2 must compare against batch 0's curve (the last real one, \
             preserved across batch 1's failure), not report NaN as though there \
             were no prior curve at all: {m2:?}"
        );
        assert!(
            m2.max_delta < 1e-6,
            "identical points to batch 0 should reproduce an identical curve, \
             got max_delta={}",
            m2.max_delta
        );
    }
}
