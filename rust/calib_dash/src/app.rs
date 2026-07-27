//! Dashboard app state: the tab cursor, the vim-style count prefix, and the
//! pause action a keypress resolves to. [`CalibDash`] ties all of it to a
//! running search, one [`CalibDash::on_batch`] call per Phase 1 scoring chunk.
//!
//! Counts use the vim prefix idiom so there is no modal input: `15n` advances
//! fifteen batches, off one `Option<u32>` that serves every motion, rather than
//! a second text-entry mode with its own key routing and escape handling.

use crate::frames::{
    CalibrantPoint,
    FrameStore,
    FrameSummary,
};
use crate::metrics::{
    BatchMetrics,
    churn,
    curve_delta,
    weighted_ridge_half_width,
};
use crate::recording::FitRecording;
use calibrt::{
    CalibRtError,
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
}

/// Which mark overlay the Fit heatmap draws on top of the density field,
/// cycled by `m`/`M`. Exactly one is active at a time: a layer is a *view*,
/// not a set of flags, so there is no priority order to maintain between
/// marks that can never appear on screen together.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Layer {
    /// Bare heatmap — the raw density the fit works from, undisturbed by
    /// anything derived from it.
    None,
    /// The DP chain and the greedily attached tails, glyph-distinguished within
    /// the layer (`O`/`X` — see `ui::mark_path`) since they answer one question
    /// together: is a bad edge the DP's own choice or a tail grafted on after.
    Path,
    Curve,
    Ridge,
    Suppressed,
}

impl Layer {
    pub const ALL: [Layer; 5] = [
        Layer::None,
        Layer::Path,
        Layer::Curve,
        Layer::Ridge,
        Layer::Suppressed,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Layer::None => "none",
            Layer::Path => "path",
            Layer::Curve => "curve",
            Layer::Ridge => "ridge",
            Layer::Suppressed => "suppressed",
        }
    }
}

/// Steps `cur` through `all` by `steps` stops, wrapping. `Tab` and `Layer` are
/// both fixed stop lists with no natural end to clamp at, so `l` past the last
/// tab is the first again, like vim's own `l`/`h` inside a fixed-width line.
///
/// `steps` is reduced mod the stop count before anything moves, so the count a
/// user can type (`push_digit` saturates at `u32::MAX`) costs one division
/// rather than that many iterations.
fn cycle<T: PartialEq + Copy, const N: usize>(
    all: &[T; N],
    cur: T,
    steps: u32,
    forward: bool,
) -> T {
    let n = N as u32;
    let at = all
        .iter()
        .position(|x| *x == cur)
        .expect("ALL lists every variant") as u32;
    let steps = steps % n;
    let offset = if forward { steps } else { n - steps };
    all[((at + offset) % n) as usize]
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

/// Skip bookkeeping, so the render-or-skip decision is testable without a
/// terminal.
struct Stepper {
    /// Batches still to skip before the next pause.
    remaining: u32,
    stopped: bool,
}

impl Stepper {
    fn new() -> Self {
        Self {
            remaining: 0,
            stopped: false,
        }
    }

    /// Consulted by `on_batch` exactly once per chunk. Takes `&mut self`
    /// because asking the question is what draws down the skip counter —
    /// splitting the query from the decrement is how a caller forgets one of
    /// them and the skip never expires.
    fn should_render(&mut self) -> bool {
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
    /// the user has already asked to stop seeing the dashboard, and nothing
    /// later should reopen it. Read-only, unlike `should_render`: the question
    /// is "has the user opted out of every future pause", not "should I render
    /// *this* batch", so it must not draw down the skip counter.
    fn is_stopped(&self) -> bool {
        self.stopped
    }

    /// Called after a pause with whatever the user chose.
    fn apply(&mut self, action: PauseAction) {
        match action {
            PauseAction::Stay => {}
            // `Next(5)` at the chunk-0 pause means "render at chunk 5", so
            // four chunks are skipped, not five.
            PauseAction::Next(n) => self.remaining = n.saturating_sub(1),
            PauseAction::RunToEnd | PauseAction::Detach | PauseAction::Abort => self.stopped = true,
        }
    }
}

/// Dashboard state for one pause: which tab is showing, the DP pane toggle,
/// the current batch number, the recording of the fit that ran at this
/// batch, and the pending vim-style count prefix.
pub struct App {
    tab: Tab,
    dp_pane: bool,
    batch: u32,
    recording: FitRecording,
    /// The post-Phase-2 recording, set once `calibration.json` is written.
    /// `None` for the whole of Phase 1: the Tolerances tab reads this to
    /// decide whether Step B has run yet at all.
    real_fit: Option<FitRecording>,
    /// The Tolerances tab's post-Step-B m/z/mobility/RT summary, set alongside
    /// `real_fit`. Lives on `App`, not `CalibDash`, for the same reason
    /// `real_fit` does: `ui.rs` only ever receives `&App`.
    tolerances: Option<ToleranceSummary>,
    /// Every batch's metrics, undecimated — see `metrics.rs`'s module doc for
    /// why the series must have no holes.
    metrics: Vec<BatchMetrics>,
    /// The three numbers the Convergence tab's header states about the frame
    /// slab's decimation. `App` owns no `FrameStore` — `CalibDash` does — so
    /// these are plain fields defaulting to "nothing decimated" and refreshed
    /// from the store through `set_frame_summary`.
    retained_frames: usize,
    frame_stride: usize,
    dropped_frames: usize,
    layer: Layer,
    /// The pending count prefix, e.g. the `15` typed before `n` in `15n`.
    /// Lives for exactly one motion: a motion consumes it via `take_count`,
    /// and any other key clears it directly, so a half-typed number never
    /// leaks into the next command.
    count: Option<u32>,
    /// Which retained frame `<`/`>` have scrubbed the Fit tab to, or `None`
    /// for the live batch. Bounded by `retained_frames`. `App` tracks only
    /// *which* frame; the recording itself lives in `scrub_recording`.
    scrub_frame: Option<usize>,
    /// The original batch/chunk number of `scrub_frame`, purely for the Fit
    /// tab's "not live" banner — a retained-frame *index* alone ("frame 2 of
    /// 5") means little to a user who thinks in batch numbers.
    scrub_chunk: Option<usize>,
    /// The scrubbed frame's refit recording. `None` whenever `scrub_frame` is,
    /// and also `None` momentarily after `scrub_frame` moves and before the
    /// next `sync_scrub` catches up — `active_recording` falls back to the live
    /// view for that gap rather than showing nothing.
    scrub_recording: Option<FitRecording>,
    /// The live batch's on-demand, `dp_nodes`-observed recording — what the DP
    /// pane reads while the Fit tab is showing the live batch. `None` when the
    /// pane is off, when the Fit tab is scrubbed, or when the live batch's
    /// points were too degenerate to refit at all.
    live_dp_recording: Option<FitRecording>,
    /// Whether the `?` key-map overlay is open. Modal in the simplest sense:
    /// while it is showing, every keypress (digits and `Ctrl-C` included) only
    /// dismisses it.
    show_keys: bool,
}

impl App {
    pub(crate) fn new(bins: usize) -> Self {
        Self {
            tab: Tab::Fit,
            dp_pane: false,
            batch: 0,
            recording: FitRecording::new(bins),
            real_fit: None,
            tolerances: None,
            metrics: Vec::new(),
            retained_frames: 0,
            frame_stride: 1,
            dropped_frames: 0,
            layer: Layer::None,
            count: None,
            scrub_frame: None,
            scrub_chunk: None,
            scrub_recording: None,
            live_dp_recording: None,
            show_keys: false,
        }
    }

    /// Mutable access to the live recording, so a caller can fit directly into
    /// the allocation `App` already owns rather than building a `FitRecording`
    /// elsewhere and having nowhere to put it.
    pub(crate) fn recording_mut(&mut self) -> &mut FitRecording {
        &mut self.recording
    }

    pub fn real_fit(&self) -> Option<&FitRecording> {
        self.real_fit.as_ref()
    }

    /// `None` exactly when `real_fit` is — `set_final` sets both together.
    pub fn tolerances(&self) -> Option<&ToleranceSummary> {
        self.tolerances.as_ref()
    }

    pub(crate) fn set_final(&mut self, recording: FitRecording, tolerances: ToleranceSummary) {
        self.real_fit = Some(recording);
        self.tolerances = Some(tolerances);
    }

    pub fn metrics(&self) -> &[BatchMetrics] {
        &self.metrics
    }

    /// Appends one batch's metrics — called every batch, rendered or not.
    pub(crate) fn push_metrics(&mut self, m: BatchMetrics) {
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

    /// Refreshes what the Convergence tab's header states about the frame
    /// slab's decimation, from the `FrameStore`'s own numbers.
    pub(crate) fn set_frame_summary(&mut self, summary: FrameSummary) {
        self.retained_frames = summary.retained;
        self.frame_stride = summary.keep_every;
        self.dropped_frames = summary.dropped;
    }

    pub fn tab(&self) -> Tab {
        self.tab
    }

    pub fn batch(&self) -> u32 {
        self.batch
    }

    pub fn dp_pane(&self) -> bool {
        self.dp_pane
    }

    pub fn layer(&self) -> Layer {
        self.layer
    }

    pub(crate) fn recording(&self) -> &FitRecording {
        &self.recording
    }

    /// What the Fit tab should actually draw: the scrubbed frame's recording
    /// while `scrub_frame` is set and `sync_scrub` has already caught up to it,
    /// the live `recording` otherwise. Never shows a blank tab for the gap
    /// between `scrub_frame` moving and the next `sync_scrub`.
    pub fn active_recording(&self) -> &FitRecording {
        if self.scrub_frame.is_some()
            && let Some(rec) = self.scrub_recording.as_ref()
        {
            rec
        } else {
            &self.recording
        }
    }

    pub fn scrub_frame(&self) -> Option<usize> {
        self.scrub_frame
    }

    pub fn scrub_chunk(&self) -> Option<usize> {
        self.scrub_chunk
    }

    /// Wires a freshly refit frame into the scrub view. The caller needs a
    /// `FrameStore` to refit from, which `App` does not have.
    pub(crate) fn set_scrub_recording(
        &mut self,
        frame: usize,
        chunk: Option<usize>,
        recording: FitRecording,
    ) {
        self.scrub_frame = Some(frame);
        self.scrub_chunk = chunk;
        self.scrub_recording = Some(recording);
    }

    /// Drops back to the live view — either the user scrubbed forward past the
    /// last retained frame, or a frame could not be refit and falling back to
    /// live beats leaving a stale or absent recording on screen.
    pub(crate) fn clear_scrub(&mut self) {
        self.scrub_frame = None;
        self.scrub_chunk = None;
        self.scrub_recording = None;
    }

    /// The Fit tab falls back to `active_recording()` (which renders "no DP
    /// trace" on its own) when this is `None`, rather than needing a
    /// special-cased blank pane.
    pub fn live_dp_recording(&self) -> Option<&FitRecording> {
        self.live_dp_recording.as_ref()
    }

    pub fn pending_count(&self) -> Option<u32> {
        self.count
    }

    pub fn show_keys(&self) -> bool {
        self.show_keys
    }

    /// Consumes the pending count, defaulting to 1 for a bare motion.
    fn take_count(&mut self) -> u32 {
        self.count.take().unwrap_or(1)
    }

    /// Folds one more digit into the pending count. A leading `0` does not
    /// start a count — matching vim, which keeps `0` free as a motion later —
    /// but `0` extends an already-started count like any other digit.
    ///
    /// Saturating, not checked: a held-down digit key is ordinary user input,
    /// so a long run must not panic (overflow checks are on in dev/test) or
    /// wrap into some other, unrelated number of batches.
    fn push_digit(&mut self, d: u32) {
        if d == 0 && self.count.is_none() {
            return;
        }
        self.count = Some(self.count.unwrap_or(0).saturating_mul(10).saturating_add(d));
    }

    /// Moves the Fit tab's frame-scrub cursor back into history by `n` retained
    /// frames, or does nothing if no frames have been retained yet. Saturating
    /// arithmetic keeps this O(1) and panic-free for any `n`, including the
    /// saturated counts `push_digit` can produce.
    fn scrub_back(&mut self, n: u32) {
        let total = self.retained_frames as u32;
        if total == 0 {
            return;
        }
        let pos = self.scrub_frame.map(|i| i as u32).unwrap_or(total);
        self.scrub_frame = Some(pos.saturating_sub(n) as usize);
    }

    /// Moves the cursor forward by `n`. Forward past the last retained frame
    /// returns to the live view rather than clamping at the last index — "keep
    /// pressing `>` and you get back to now".
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

    /// Routes one keypress. Every bare-character binding lives in
    /// `handle_plain_char`, behind a single `modifiers.is_empty()` guard, so a
    /// modified key can never fall into a bare-letter arm by accident. Only the
    /// three deliberately lax bindings stay here.
    pub(crate) fn handle_key(&mut self, key: KeyEvent) -> PauseAction {
        // The `?` overlay is modal: any keypress only dismisses it, so a user
        // who presses `n` to dismiss it does not also skip a batch.
        // Deliberately ahead of the `Ctrl-C` check below: while the overlay is
        // up, dismissing it is the whole of what *any* key does, abort
        // included.
        if self.show_keys {
            self.show_keys = false;
            return PauseAction::Stay;
        }
        if key.code == KeyCode::Char('c') && key.modifiers.contains(KeyModifiers::CONTROL) {
            self.count = None;
            return PauseAction::Abort;
        }

        match key.code {
            // `?` and `M` are shifted characters and terminals disagree on
            // whether they report `Shift` on the resulting `Char`, so these two
            // guard only against `Control`. They precede the plain-character
            // arm so an unshifted-looking `?`/`M` still lands here.
            KeyCode::Char('?') if !key.modifiers.contains(KeyModifiers::CONTROL) => {
                self.show_keys = true;
                self.count = None;
                PauseAction::Stay
            }
            KeyCode::Char('M') if !key.modifiers.contains(KeyModifiers::CONTROL) => {
                self.layer = cycle(&Layer::ALL, self.layer, self.take_count(), false);
                PauseAction::Stay
            }
            KeyCode::Char(c) if key.modifiers.is_empty() => self.handle_plain_char(c),
            // Esc lands here, and needs no modifier guard of its own: it
            // clears the count and stays either way, exactly what every other
            // unbound key does.
            _ => {
                self.count = None;
                PauseAction::Stay
            }
        }
    }

    /// The unmodified-character keymap: digits build the count prefix, motions
    /// consume it via `take_count`, and every other key clears it.
    fn handle_plain_char(&mut self, c: char) -> PauseAction {
        if let Some(d) = c.to_digit(10) {
            self.push_digit(d);
            return PauseAction::Stay;
        }
        match c {
            'n' => {
                let n = self.take_count();
                // `n` advances the *live* batch loop, so a scrub cursor left
                // over from browsing history at this pause must not still be
                // showing at the next one: reopening on a replayed frame right
                // after the user asked to advance reads as confusing.
                self.clear_scrub();
                PauseAction::Next(n)
            }
            'd' => {
                self.dp_pane = !self.dp_pane;
                if !self.dp_pane {
                    // Closing the pane leaves the last recording with no
                    // reader; drop it rather than let `sync_dp` see a stale one
                    // if the pane is reopened on a later batch.
                    self.live_dp_recording = None;
                }
                self.count = None;
                PauseAction::Stay
            }
            'l' => {
                self.tab = cycle(&Tab::ALL, self.tab, self.take_count(), true);
                PauseAction::Stay
            }
            'h' => {
                self.tab = cycle(&Tab::ALL, self.tab, self.take_count(), false);
                PauseAction::Stay
            }
            // `m` for "mark", cycling forward; shifted (`M`, in `handle_key`)
            // cycles back, the same pairing `l`/`h` uses for tabs.
            'm' => {
                self.layer = cycle(&Layer::ALL, self.layer, self.take_count(), true);
                PauseAction::Stay
            }
            '<' => {
                let n = self.take_count();
                self.scrub_back(n);
                PauseAction::Stay
            }
            '>' => {
                let n = self.take_count();
                self.scrub_forward(n);
                PauseAction::Stay
            }
            'r' => {
                self.count = None;
                PauseAction::RunToEnd
            }
            'q' => {
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

/// A batch's ridge overlay is the same computation Phase 2 runs, not a
/// dashboard-only approximation, so the threshold comes from calibrt rather
/// than being spelled here.
const RIDGE_FRACTION: f64 = calibrt::DEFAULT_RIDGE_FRACTION;

/// How many library-RT samples `curve_delta` draws between batches' curves. Any
/// `>= 2` produces a defined delta; this is high enough that a local kink is
/// unlikely to fall between sample points.
const CURVE_DELTA_SAMPLES: usize = 50;

/// The `(x_range, y_range)` pair `point_ranges` derives and
/// `CalibrationState::reconfigure` is configured with.
type GridRanges = ((f64, f64), (f64, f64));

/// The min/max of `library_rt` and `observed_rt` over a point set, the same fold
/// `timsseek_cli::processing` runs before building the real `CalibrationState`.
///
/// `None` when every point shares one coordinate (or there are no points at
/// all): a `CalibrationState` cannot be configured with a zero-width range.
fn point_ranges(points: &[CalibrantPoint]) -> Option<GridRanges> {
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

/// Why a re-fit produced nothing. Only `refit_live` reports these (through
/// `Display`); the other callers treat any of them as "no recording".
#[derive(Debug)]
enum RefitSkipped {
    /// The points span a zero-width range on at least one axis, so there is no
    /// grid to configure — see [`point_ranges`].
    DegenerateRange,
    /// The named [`fit_points`] stage was refused by `calibrt`.
    Failed {
        stage: &'static str,
        error: CalibRtError,
    },
}

impl std::fmt::Display for RefitSkipped {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RefitSkipped::DegenerateRange => {
                f.write_str("the points span a zero-width range on at least one axis")
            }
            // `CalibRtError` has no `Display`.
            RefitSkipped::Failed { stage, error } => {
                write!(f, "{stage} failed: {error:?}")
            }
        }
    }
}

/// The whole re-fit sequence, in the one place every caller shares:
/// `reconfigure` → `update` → `fit_with` → `measure_ridge_width_with`, against
/// whichever `CalibrationState`/`FitRecording` pair the caller owns.
///
/// The grid ranges are derived *here*, by [`point_ranges`], rather than passed
/// in: the live fit and a scrubbed frame's replay agreeing on that fold is what
/// makes a replayed frame show the fit that actually ran. They are returned for
/// callers that need them afterwards (`refit_live` compares curves over the x
/// range it just fit).
fn fit_points(
    state: &mut CalibrationState,
    recording: &mut FitRecording,
    bins: usize,
    points: &[CalibrantPoint],
    opts: ObserveOpts,
) -> Result<GridRanges, RefitSkipped> {
    let (x_range, y_range) = point_ranges(points).ok_or(RefitSkipped::DegenerateRange)?;
    state
        .reconfigure(bins, x_range, y_range)
        .map_err(|error| RefitSkipped::Failed {
            stage: "reconfiguring the calibration state",
            error,
        })?;
    state
        .update(points.iter().map(as_calibrt_tuple))
        .map_err(|error| RefitSkipped::Failed {
            stage: "updating the calibration state",
            error,
        })?;
    state.fit_with(recording, opts);
    state.measure_ridge_width_with(RIDGE_FRACTION, recording);
    Ok((x_range, y_range))
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
/// from.
///
/// `mz_ppm` and `mobility_pct` are signed `(lo, hi)` windows around zero and are
/// genuinely asymmetric — Step B can derive a different tolerance on either side
/// of a calibrant. `rt_seconds` is symmetric by construction, and is a bare
/// half-width rather than a third `(lo, hi)`: a tuple would carry the same number
/// twice under a type that says the two sides may differ, and any renderer
/// reading only `.0` would silently drop the other.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ToleranceSummary {
    pub mz_ppm: (f64, f64),
    pub mobility_pct: (f64, f64),
    /// Half-width: the tolerance is `±rt_seconds`.
    pub rt_seconds: f64,
    pub n_calibrants: usize,
}

/// Ties the dashboard to a running Phase 1 batch loop, through `on_batch` once
/// per scoring chunk.
///
/// Every buffer here — the `CalibrationState`s (reused via `reconfigure`, since
/// `bins` never changes), their `FitRecording`s, the frame slab,
/// `current_points`, `prev_points` — is allocated at startup and reused, so
/// steady state adds nothing per batch beyond the one `CalibrationCurve`
/// `calibrt` allocates for every fit anyway (plus the `prev_curve` clone
/// `curve_delta` compares against).
pub struct CalibDash {
    app: App,
    frames: FrameStore,
    stepper: Stepper,
    /// The live re-fit's `CalibrationState`, reused every batch.
    state: CalibrationState,
    /// This batch's points, buffered once so recording, re-fitting and churn all
    /// read the same materialized slice rather than each consuming a fresh copy
    /// of the caller's iterator. Clamped to `n_calibrants` exactly the way
    /// `FrameStore::record` clamps its own copy: without matching clamps, an
    /// overlong batch would fit the live curve against points the frame slab
    /// never recorded, so replaying that frame would disagree with what was
    /// shown live.
    current_points: Vec<CalibrantPoint>,
    n_calibrants: usize,
    /// The previous batch's points, for `churn`. Empty before the first batch —
    /// `churn`'s own "everything admitted, nothing evicted" case.
    prev_points: Vec<CalibrantPoint>,
    /// The previous batch's fitted curve, for `curve_delta`. `None` before the
    /// first *successful* fit: a batch whose fit fails must not reset this, or
    /// the next successful batch would report a NaN delta as though there had
    /// never been a curve at all.
    prev_curve: Option<CalibrationCurve>,
    /// Held so `refit_frame` can reconfigure `refit_state` to the same grid
    /// resolution the live fit used.
    bins: usize,
    /// A separate `CalibrationState`/`FitRecording` pair for replaying a
    /// recorded frame on demand (scrubbing), so stepping through history never
    /// perturbs the live recording.
    refit_state: CalibrationState,
    refit_recording: FitRecording,
    /// A third pair, for the on-demand `dp_nodes: true` re-fit of the *live*
    /// batch that pressing `d` triggers. Kept separate from the scrubber's
    /// `refit_state`/`refit_recording` so the two features cannot clobber each
    /// other's buffer within one draw: `sync_dp` already skips while scrubbing,
    /// so today they never overlap, and the separate buffer is what keeps that a
    /// property of the code rather than of call order.
    dp_state: CalibrationState,
    dp_recording: FitRecording,
    /// Which batch `dp_recording` was last computed for, so `sync_dp` can skip
    /// the refit on every keystroke of the same pause.
    dp_recording_batch: Option<u32>,
    /// Counts calls into `render_pause` — test-only instrumentation so
    /// `present`'s stepper guard is provable without a real terminal: under
    /// `cargo test`, `render_pause` detaches immediately whether or not the
    /// guard let it through, so a test asserting `present()` "did nothing"
    /// cannot tell the two apart by output alone.
    #[cfg(test)]
    render_pause_calls: u32,
    /// Counts `dp_refit_current` calls, for the same reason: `sync_dp`'s
    /// `dp_recording_batch` cache has no other observable — a skipped refit and
    /// a repeated one leave `live_dp_recording` looking identical.
    #[cfg(test)]
    dp_refit_calls: u32,
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
        // reconfigures every state to the batch's own point-derived ranges
        // before anything is fit against them. Being fixed and nonzero-width, it
        // leaves `bins == 0` as the only way `CalibrationState::new` can fail
        // here. This constructor has no `Result` to report that through; the
        // standalone binary's `validate_snapshot` catches it first.
        let placeholder = (0.0, 1.0);
        let new_state = || {
            CalibrationState::new(bins, placeholder, placeholder, lookback)
                .expect("bins must be nonzero (the placeholder range is always valid)")
        };
        Self {
            app: App::new(bins),
            frames: FrameStore::new(n_frames, n_calibrants, budget_bytes),
            stepper: Stepper::new(),
            state: new_state(),
            current_points: Vec::with_capacity(n_calibrants),
            n_calibrants,
            prev_points: Vec::new(),
            prev_curve: None,
            bins,
            refit_state: new_state(),
            refit_recording: FitRecording::new(bins),
            dp_state: new_state(),
            dp_recording: FitRecording::new(bins),
            dp_recording_batch: None,
            #[cfg(test)]
            render_pause_calls: 0,
            #[cfg(test)]
            dp_refit_calls: 0,
        }
    }

    /// In order: record the frame, re-fit the live curve from those same points,
    /// push this batch's metrics (every batch, rendered or not — see
    /// `metrics.rs`'s module doc), then either continue or pause to render.
    ///
    /// `Flow::Abort` is the one path that asks the search to stop, and it happens
    /// only because the user pressed Ctrl-C at a pause: a dashboard that cannot
    /// fit or cannot draw never fails a search.
    pub fn on_batch(&mut self, chunk: usize, points: impl Iterator<Item = CalibrantPoint>) -> Flow {
        self.current_points.clear();
        self.current_points.extend(points.take(self.n_calibrants));

        self.frames
            .record(chunk, self.current_points.iter().copied());
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

    /// Opens the post-Phase-2 pause, so the Tolerances tab is reachable once
    /// more after Step B runs. There is no batch loop left to steer, so every
    /// `PauseAction` collapses to "return to the caller" — `Abort` included.
    ///
    /// Guarded on `stepper.is_stopped()`: a user who chose `r`/`q`/Ctrl-C at any
    /// Phase 1 pause has already asked to stop seeing the dashboard, and the
    /// caller reaches this call regardless of `Flow::Abort`. Unguarded, it would
    /// reopen the alternate screen and block on a keypress that never comes.
    pub fn present(&mut self) {
        if self.stepper.is_stopped() {
            return;
        }
        render_pause(self);
    }

    /// Promotes the reserved final frame and refreshes the Convergence header's
    /// decimation numbers one last time — the promotion can change them with no
    /// further batch left to call `sync_frame_summary`.
    pub fn finish(&mut self) {
        self.frames.finish();
        self.sync_frame_summary();
    }

    /// Wires the post-Phase-2 recording and derived tolerances into the
    /// dashboard once Step B has run. Both land on `App` (not `CalibDash`) since
    /// `ui.rs` only ever receives a plain `&App`.
    pub fn show_final(&mut self, recording: FitRecording, tolerances: ToleranceSummary) {
        self.app.set_final(recording, tolerances);
    }

    /// Re-fits a recorded frame's points from scratch, into `refit_state`/
    /// `refit_recording` rather than the live pair, through the same
    /// `fit_points` sequence `refit_live` runs — so what this reproduces is the
    /// batch that actually ran, not an approximation of it. `None` if the frame
    /// index doesn't exist or the re-fit was skipped.
    fn refit_frame(&mut self, i: usize) -> Option<&FitRecording> {
        let (_chunk, pts) = self.frames.frame(i)?;
        fit_points(
            &mut self.refit_state,
            &mut self.refit_recording,
            self.bins,
            pts,
            ObserveOpts { dp_nodes: true },
        )
        .ok()?;
        Some(&self.refit_recording)
    }

    /// Re-fits the live `CalibrationState` from `current_points` and returns
    /// this batch's `BatchMetrics`. The re-fit and the metrics live together
    /// because the metrics *are* a comparison of this re-fit against the last
    /// one.
    fn refit_live(&mut self, chunk: usize) -> BatchMetrics {
        let n_points = self.current_points.len();
        let mut wrmse = f64::NAN;
        let mut max_delta = f64::NAN;
        let mut mean_delta = f64::NAN;
        let mut path_nodes = 0usize;
        let mut ridge_half_width = f64::NAN;

        match fit_points(
            &mut self.state,
            self.app.recording_mut(),
            self.bins,
            &self.current_points,
            ObserveOpts::NONE,
        ) {
            Err(skipped) => {
                tracing::warn!(
                    chunk,
                    "calib_dash: skipping this batch's re-fit ({skipped}); its metrics \
                     report n_points only"
                );
            }
            Ok((x_range, _y_range)) => {
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
                    // produced a curve: a failed fit must not reset
                    // `prev_curve` to `None`, or the *next* successful batch
                    // would report a NaN delta as if there were no prior curve
                    // at all — a discontinuity on the Convergence tab that
                    // never actually happened.
                    self.prev_curve = Some(curve.clone());
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

    /// Pushes `frames`' real decimation numbers into `app`, so
    /// `retained_frames`/`frame_stride`/`dropped_frames` have one source
    /// (`FrameStore`) rather than a second tracker duplicating it.
    fn sync_frame_summary(&mut self) {
        self.app.set_frame_summary(self.frames.summary());
    }

    /// Keeps the Fit tab's scrub view in lockstep with `App::scrub_frame`.
    /// `App` moves `scrub_frame` on `<`/`>` but has no `FrameStore` to refit
    /// from, so this is the other half of that key handling.
    ///
    /// A no-op when `scrub_frame` is `None` (the live view, the common case).
    /// Otherwise re-fits and clones unconditionally rather than caching whether
    /// the index already matches what was last shown: this runs at most once per
    /// keypress inside an interactive pause, never per batch, so the clone is
    /// off `on_batch`'s steady-state path entirely.
    fn sync_scrub(&mut self) {
        let Some(i) = self.app.scrub_frame() else {
            return;
        };
        let chunk = self.frames.frame(i).map(|(chunk, _)| chunk);
        match self.refit_frame(i) {
            Some(rec) => {
                let cloned = rec.clone();
                self.app.set_scrub_recording(i, chunk, cloned);
            }
            // The index is bounded by `retained_frames`, set from this same
            // `FrameStore`, so this only fires on a frame too degenerate to
            // refit — fall back to live rather than leave a stale recording up.
            None => self.app.clear_scrub(),
        }
    }

    /// Re-fits `current_points` into `dp_state`/`dp_recording` with
    /// `ObserveOpts { dp_nodes: true }`, so the DP pane has decisions to show.
    /// Only ever called from `sync_dp`, never from `on_batch`: that separation is
    /// what keeps the per-DP-node observation cost off the per-batch hot path.
    /// `None` on the same degenerate-range skip `refit_live` makes.
    fn dp_refit_current(&mut self) -> Option<&FitRecording> {
        #[cfg(test)]
        {
            self.dp_refit_calls += 1;
        }
        fit_points(
            &mut self.dp_state,
            &mut self.dp_recording,
            self.bins,
            &self.current_points,
            ObserveOpts { dp_nodes: true },
        )
        .ok()?;
        Some(&self.dp_recording)
    }

    /// Keeps `App::live_dp_recording` in sync with `App::dp_pane`. A no-op
    /// whenever the pane is off, or the Fit tab is showing a scrubbed frame:
    /// `refit_frame` already observes DP nodes for every scrubbed frame, so
    /// recomputing here would at best duplicate that work and at worst show the
    /// *live* batch's DP trace on a Fit tab showing an older, scrubbed one.
    ///
    /// Otherwise re-fits once per batch, not once per draw.
    fn sync_dp(&mut self) {
        if !self.app.dp_pane() || self.app.scrub_frame().is_some() {
            return;
        }
        // The cache is only good while the recording it names is still around:
        // closing the pane drops the recording, so a reopen within the same
        // batch has to refit rather than serve a cache entry pointing at
        // nothing.
        if self.dp_recording_batch == Some(self.app.batch())
            && self.app.live_dp_recording().is_some()
        {
            return;
        }
        let refit = self.dp_refit_current().cloned();
        self.app.live_dp_recording = refit;
        self.dp_recording_batch = Some(self.app.batch());
    }
}

/// Pauses the batch loop to render one interactive frame and block until the
/// user's keypress resolves to something other than [`PauseAction::Stay`]. Takes
/// the whole [`CalibDash`] because the batch scrubber needs `refit_frame`'s
/// `FrameStore` access on every redraw.
///
/// Never fails the caller: with stdout not a terminal (under `cargo test`, or a
/// piped run) this returns [`PauseAction::Detach`] immediately rather than block
/// on a keypress that never arrives, and a terminal that fails to initialize is
/// logged at `warn` and treated the same way.
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

/// Runs `f`, converting a panic into `PauseAction::Detach` rather than letting
/// it unwind into the search — in builds where panics unwind at all. This
/// workspace's `[profile.release]` sets `panic = "abort"`, so the shipped binary
/// aborts before `catch_unwind` gets a chance; the recovery is real in dev/test.
/// Factored out of `render_pause` so it is testable without a real terminal.
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

/// Draws and reads key events in the alternate screen until a keypress resolves
/// to something other than `Stay`. A draw or read failure is treated the same as
/// `render_pause`'s terminal-setup failure: logged and detached, never propagated
/// to fail the batch loop.
fn event_loop<B: ratatui::backend::Backend>(
    terminal: &mut ratatui::Terminal<B>,
    dash: &mut CalibDash,
) -> PauseAction {
    loop {
        // Before every draw, not just after a scrub keypress: simpler than
        // tracking whether `scrub_frame` changed, and a no-op when it's `None`.
        dash.sync_scrub();
        dash.sync_dp();
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
    use calibrt::Point;
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

    /// A dashboard with a budget generous enough to retain every frame. Tests
    /// that want decimation spell `CalibDash::new` out with a tiny budget.
    fn dash(n_frames: usize, n_calibrants: usize, bins: usize, lookback: usize) -> CalibDash {
        CalibDash::new(n_frames, n_calibrants, bins, lookback, 1 << 20)
    }

    /// `n` points on the line `observed_rt = library_rt * slope`, weighted
    /// `1.0 + i` — comfortably above `suppress_nonmax`'s 1.0 seed, so the fit
    /// succeeds.
    fn points(n: usize, slope: f64, chunk: usize) -> Vec<CalibrantPoint> {
        scored_points(n, slope, chunk, |i| 1.0 + i as f64)
    }

    /// `points`, with the score picked per point — the only knob the
    /// suppression-threshold tests need beyond the line itself.
    fn scored_points(
        n: usize,
        slope: f64,
        chunk: usize,
        score: impl Fn(usize) -> f64,
    ) -> Vec<CalibrantPoint> {
        (0..n)
            .map(|i| CalibrantPoint {
                library_rt: i as f64 + 0.5,
                observed_rt: (i as f64 + 0.5) * slope,
                score: score(i),
                speclib_index: chunk * n + i,
            })
            .collect()
    }

    /// Points at explicit `speclib_index`es, so a test can choose the two sets
    /// `churn` diffs independently of how many batches have gone by. On the
    /// identity line with weights above `suppress_nonmax`'s 1.0 seed, so the fit
    /// succeeds like `points`'.
    fn indexed_points(indices: &[usize]) -> Vec<CalibrantPoint> {
        indices
            .iter()
            .enumerate()
            .map(|(i, &idx)| CalibrantPoint {
                library_rt: i as f64 + 0.5,
                observed_rt: i as f64 + 0.5,
                score: 1.0 + i as f64,
                speclib_index: idx,
            })
            .collect()
    }

    /// `n` retained frames, nothing decimated — what the scrubber tests need
    /// from a `FrameStore` they do not have.
    fn retained(n: usize) -> FrameSummary {
        FrameSummary {
            retained: n,
            keep_every: 1,
            dropped: 0,
        }
    }

    /// Pointwise curve equality — the assertion behind every "what you scrub is
    /// what ran" test.
    fn assert_same_curve(live: &[Point], refit: &[Point]) {
        assert_eq!(live.len(), refit.len(), "curve lengths differ");
        for (a, b) in live.iter().zip(refit) {
            assert!((a.library - b.library).abs() < 1e-12, "x differs");
            assert!((a.observed - b.observed).abs() < 1e-12, "y differs");
        }
    }

    // ---- count prefix ----

    #[test]
    fn a_motion_consumes_the_count() {
        let mut app = App::new(10);
        assert_eq!(press(&mut app, "15n"), PauseAction::Next(15));
        assert_eq!(app.pending_count(), None, "count lives exactly one motion");
    }

    /// `0` is the one digit whose meaning depends on what came before it: it
    /// cannot *start* a count (vim keeps `0` free as a motion) but extends one
    /// like any other digit.
    #[test]
    fn a_leading_zero_does_not_start_a_count_but_a_later_zero_extends_one() {
        let mut app = App::new(10);
        app.handle_key(key('0'));
        assert_eq!(app.pending_count(), None, "a leading 0 starts nothing");

        press(&mut app, "10");
        assert_eq!(
            app.pending_count(),
            Some(10),
            "0 after a 1 extends the count rather than being swallowed"
        );
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

    /// Both halves of "only a motion consumes a count": a bound non-motion
    /// (`d`, which reaches `handle_plain_char`) and an unbound key (`Esc`,
    /// which falls through `handle_key`'s catch-all arm) each clear it.
    #[test]
    fn a_non_motion_key_clears_a_pending_count() {
        let mut app = App::new(10);
        press(&mut app, "42d");
        assert_eq!(
            app.pending_count(),
            None,
            "`d` toggles the DP pane, not a motion"
        );

        press(&mut app, "42");
        assert_eq!(app.pending_count(), Some(42), "sanity: the count is back");
        app.handle_key(KeyEvent::new(KeyCode::Esc, KeyModifiers::NONE));
        assert_eq!(
            app.pending_count(),
            None,
            "Esc is unbound, and every unbound key clears the count too"
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
    fn a_modified_key_is_not_a_motion() {
        let mut app = App::new(10);
        let before = app.tab();
        app.handle_key(ctrl('l'));
        assert_eq!(app.tab(), before, "modified keys are not motions");
    }

    // ---- the `?` key-map overlay ----

    #[test]
    fn question_mark_opens_the_overlay_and_clears_a_pending_count() {
        let mut app = App::new(10);
        press(&mut app, "15");
        assert_eq!(app.pending_count(), Some(15));
        let action = app.handle_key(key('?'));
        assert!(app.show_keys());
        assert_eq!(action, PauseAction::Stay);
        assert_eq!(
            app.pending_count(),
            None,
            "`?` is a non-motion key, like `d`/`r`/`q`, so it clears a pending count"
        );
    }

    #[test]
    fn any_key_dismisses_the_overlay_without_also_acting_on_it() {
        let mut app = App::new(10);
        app.handle_key(key('?'));
        assert!(app.show_keys());

        // `n` would normally advance the batch loop; while the overlay is
        // open it must only dismiss the overlay.
        let action = app.handle_key(key('n'));
        assert!(!app.show_keys(), "the overlay must close on the next key");
        assert_eq!(
            action,
            PauseAction::Stay,
            "the dismissing key must not also perform its usual action"
        );

        // A digit takes a different route through `handle_plain_char` than a
        // motion does, so it needs its own case: `5` must not start a pending
        // count on its way out of the overlay either.
        app.handle_key(key('?'));
        let action = app.handle_key(key('5'));
        assert!(!app.show_keys(), "a digit dismisses the overlay too");
        assert_eq!(action, PauseAction::Stay);
        assert_eq!(
            app.pending_count(),
            None,
            "the dismissing digit must not start a count"
        );
    }

    #[test]
    fn ctrl_c_while_the_overlay_is_open_only_dismisses_it() {
        let mut app = App::new(10);
        app.handle_key(key('?'));
        assert_eq!(app.handle_key(ctrl('c')), PauseAction::Stay);
        assert!(!app.show_keys());
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
        app.set_frame_summary(retained(5));
        press(&mut app, "<");
        assert_eq!(
            app.scrub_frame(),
            Some(4),
            "one step back from live lands on the last retained frame (index 4 of 5)"
        );
    }

    #[test]
    fn less_than_clamps_at_the_first_frame_rather_than_underflowing() {
        let mut app = App::new(10);
        app.set_frame_summary(retained(5));
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
        app.set_frame_summary(retained(5));
        press(&mut app, "3<"); // frame 2 (of 5)
        assert_eq!(app.scrub_frame(), Some(2));
        press(&mut app, ">");
        assert_eq!(
            app.scrub_frame(),
            Some(3),
            "a bare > moves forward exactly one retained frame"
        );
        press(&mut app, "2>");
        assert_eq!(
            app.scrub_frame(),
            None,
            "> past the last retained frame returns to the live view"
        );

        // A count large enough to overflow the `+ n` must saturate rather than
        // wrap into some other index — the same run must also return promptly,
        // this being arithmetic and not a loop.
        press(&mut app, "<"); // back to frame 4
        assert_eq!(app.scrub_frame(), Some(4), "sanity");
        press(&mut app, "4294967295>");
        assert_eq!(
            app.scrub_frame(),
            None,
            "a saturated count lands past the end, back to live"
        );
    }

    #[test]
    fn greater_than_is_a_no_op_while_already_live() {
        let mut app = App::new(10);
        app.set_frame_summary(retained(5));
        press(&mut app, ">");
        assert_eq!(app.scrub_frame(), None);
    }

    #[test]
    fn n_clears_a_scrub_cursor_left_over_from_this_pause() {
        let mut app = App::new(10);
        app.set_frame_summary(retained(5));
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

    // ---- Fit-tab mark layer: m / M ----

    #[test]
    fn m_cycles_the_mark_layer_forward_through_all_five_stops_and_wraps() {
        let mut app = App::new(10);
        assert_eq!(app.layer(), Layer::None);
        press(&mut app, "m");
        assert_eq!(app.layer(), Layer::Path);
        press(&mut app, "m");
        assert_eq!(app.layer(), Layer::Curve);
        press(&mut app, "m");
        assert_eq!(app.layer(), Layer::Ridge);
        press(&mut app, "m");
        assert_eq!(app.layer(), Layer::Suppressed);
        press(&mut app, "m");
        assert_eq!(app.layer(), Layer::None, "wraps back to the first stop");
    }

    #[test]
    fn capital_m_cycles_the_mark_layer_backward_and_wraps() {
        let mut app = App::new(10);
        press(&mut app, "M");
        assert_eq!(
            app.layer(),
            Layer::Suppressed,
            "one step back from the first stop wraps to the last"
        );
    }

    /// `cycle` reduces the count mod the number of stops before anything moves,
    /// so a ten-digit count returns promptly rather than driving ~10^9
    /// iterations on the input thread. Both counts below land two stops from the
    /// start (`1_000_000_001 % 3 == 2`, `1_000_000_002 % 5 == 2`).
    ///
    /// The backward keys are the ones that make the reduction load-bearing
    /// rather than merely fast: `cycle` walks back by `n - steps`, which
    /// underflows for any count past the stop list's length — and a count past
    /// 5 is an ordinary keystroke, not a stress case.
    #[test]
    fn a_huge_count_before_a_cycling_key_does_not_spin_or_underflow() {
        let mut app = App::new(10);
        press(&mut app, "1000000001l");
        assert_eq!(app.tab(), Tab::Tolerances, "Fit + 2 of 3 tab stops");
        press(&mut app, "1000000001h");
        assert_eq!(app.tab(), Tab::Fit, "and 2 of 3 stops back again");

        let mut app = App::new(10);
        press(&mut app, "1000000002m");
        assert_eq!(app.layer(), Layer::Curve, "None + 2 of 5 layer stops");
        press(&mut app, "1000000002M");
        assert_eq!(app.layer(), Layer::None, "and 2 of 5 stops back again");
    }

    // ---- the stepper ----

    #[test]
    fn next_n_renders_on_the_nth_batch_not_the_n_plus_first() {
        let mut s = Stepper::new();
        // At the chunk-0 pause the user asks for 5.
        s.apply(PauseAction::Next(5));
        // Chunks 1..=4 are skipped, the 5th call renders, and the skip is spent
        // once rather than re-armed — so every later chunk renders too.
        let rendered: Vec<usize> = (1..=10).filter(|_| s.should_render()).collect();
        assert_eq!(
            rendered,
            vec![5, 6, 7, 8, 9, 10],
            "chunks 1-4 skipped, render at 5, then every batch again"
        );
    }

    #[test]
    fn a_stopping_action_never_renders_again() {
        for action in [
            PauseAction::RunToEnd,
            PauseAction::Detach,
            PauseAction::Abort,
        ] {
            let mut s = Stepper::new();
            s.apply(action);
            assert!((0..1000).all(|_| !s.should_render()), "{action:?}");
        }
    }

    // ---- CalibDash::on_batch ----

    #[test]
    fn metrics_are_collected_for_every_batch_including_skipped_ones() {
        let mut d = CalibDash::new(10, 8, 10, 5, 64); // tiny budget: heavy decimation
        for chunk in 0..10 {
            // stdout is not a tty under `cargo test`, so the render path must
            // short-circuit rather than wait for a keypress — and must never
            // abort the search on its own.
            let flow = d.on_batch(chunk, points(8, 1.0, chunk).into_iter());
            assert!(matches!(flow, Flow::Continue), "chunk {chunk}");
        }
        d.finish();
        assert_eq!(d.app.metrics().len(), 10, "metrics are undecimated");
        assert!(d.frames.summary().retained < 10, "frames are decimated");
    }

    /// `BatchMetrics.admitted`/`evicted` are the only two fields nothing but
    /// `refit_live` ever fills, so the wiring needs its own test: `metrics.rs`
    /// exercises `churn` directly and `ui.rs` hand-builds `BatchMetrics`.
    ///
    /// The second batch is deliberately asymmetric (one index in, two out), so
    /// swapping the two fields at the construction site cannot pass. It also
    /// fails if the `prev_points` refresh moves above the `churn` call, which
    /// would diff `current_points` against itself and report `(0, 0)` forever.
    #[test]
    fn on_batch_reports_this_batchs_admissions_and_evictions() {
        let mut d = dash(2, 8, 10, 5);

        d.on_batch(0, indexed_points(&[0, 1, 2, 3]).into_iter());
        let m0 = d.app.metrics()[0];
        assert_eq!(
            (m0.admitted, m0.evicted),
            (4, 0),
            "the first batch admits all four against an empty baseline: {m0:?}"
        );

        d.on_batch(1, indexed_points(&[2, 3, 4]).into_iter());
        let m1 = d.app.metrics()[1];
        assert_eq!(
            (m1.admitted, m1.evicted),
            (1, 2),
            "index 4 admitted, indices 0 and 1 evicted: {m1:?}"
        );
    }

    /// End-to-end version of the `<` scrubber: moving `App::scrub_frame` with
    /// `handle_key` and then letting `CalibDash::sync_scrub` catch up must leave
    /// the Fit tab looking at the earlier batch's real curve, not the live one.
    #[test]
    fn sync_scrub_reproduces_an_earlier_batchs_curve_on_the_fit_tab() {
        let mut d = dash(3, 8, 10, 5);
        for chunk in 0..3 {
            let slope = 1.0 + chunk as f64 * 0.1;
            d.on_batch(chunk, points(8, slope, chunk).into_iter());
        }
        d.finish();

        // Batch 0's curve, refit independently for comparison against what
        // scrubbing shows.
        let batch0_curve = d.refit_frame(0).expect("frame 0 exists").curve().to_vec();

        // A generous budget retains all 3 batches; scrub 3 back from live
        // (batch 2) lands on frame index 0 (batch 0).
        press(&mut d.app, "3<");
        d.sync_scrub();

        assert_eq!(d.app.scrub_frame(), Some(0));
        assert_eq!(d.app.scrub_chunk(), Some(0));
        assert_same_curve(d.app.active_recording().curve(), &batch0_curve);

        // The live view is unaffected by scrubbing — `recording()` (not
        // `active_recording()`) still reflects batch 2's own fit.
        assert!(
            !d.app.recording().curve().is_empty(),
            "sanity: batch 2 produced a real curve too"
        );
    }

    /// A *live* `on_batch` followed by pressing `d` must leave the pane with
    /// decisions to show. The live re-fit deliberately uses `ObserveOpts::NONE`,
    /// so `recording().dp()` staying empty is expected; the DP data has to come
    /// from `sync_dp`'s separate refit, called directly here since there is no
    /// terminal to run `event_loop` under `cargo test`.
    #[test]
    fn d_on_a_live_batch_refits_with_dp_nodes_so_the_pane_has_decisions() {
        let mut d = dash(2, 8, 10, 5);
        d.on_batch(0, points(8, 1.0, 0).into_iter());

        assert!(
            d.app.recording().dp().is_empty(),
            "sanity: the live per-batch fit must not pay the DP-node cost"
        );
        assert!(
            d.app.live_dp_recording().is_none(),
            "sanity: nothing has refit on demand yet"
        );

        assert_eq!(press(&mut d.app, "d"), PauseAction::Stay);
        assert!(d.app.dp_pane(), "d must turn the pane on");

        d.sync_dp();

        let rec = d
            .app
            .live_dp_recording()
            .expect("d on the live batch must have triggered an on-demand refit");
        assert!(
            !rec.dp().is_empty(),
            "the DP pane must have decisions to show after pressing d on a live batch, \
             not an empty recording"
        );
    }

    /// Pressing `d` while scrubbed to an earlier frame must not trigger the
    /// on-demand live refit at all — `active_recording()` already carries DP
    /// data in that case, and `sync_dp` recomputing it from `current_points`
    /// would show the *live* batch's decisions on a Fit tab that is actually
    /// showing a different, scrubbed one.
    #[test]
    fn d_while_scrubbed_does_not_touch_the_live_on_demand_recording() {
        let mut d = dash(2, 8, 10, 5);
        for chunk in 0..2 {
            let slope = 1.0 + chunk as f64 * 0.1;
            d.on_batch(chunk, points(8, slope, chunk).into_iter());
        }
        d.finish();

        press(&mut d.app, "<");
        d.sync_scrub();
        assert_eq!(
            d.app.scrub_frame(),
            Some(1),
            "sanity: one step back from live lands on the last retained frame (index 1 of 2)"
        );
        // `active_recording()` already has DP data (`refit_frame` always
        // observes it) without anyone pressing `d` at all.
        assert!(!d.app.active_recording().dp().is_empty());

        press(&mut d.app, "d");
        d.sync_dp();
        assert!(
            d.app.live_dp_recording().is_none(),
            "sync_dp must stay a no-op while scrubbed, not refit the live batch"
        );
    }

    /// `sync_dp` runs before every draw, i.e. once per keystroke of a pause, but
    /// must refit at most once per batch — and the cache that makes that true
    /// must not also serve batch N's DP trace at batch N+1.
    #[test]
    fn sync_dp_refits_once_per_batch_and_never_serves_a_stale_trace() {
        let mut d = dash(3, 8, 10, 5);
        d.on_batch(0, points(8, 1.0, 0).into_iter());
        press(&mut d.app, "d");

        for _ in 0..3 {
            d.sync_dp();
        }
        assert_eq!(
            d.dp_refit_calls, 1,
            "three draws within one pause must share one refit"
        );

        // A new batch invalidates the cache: the pane must show *this* batch's
        // decisions, which for a different slope is a visibly different curve.
        d.on_batch(1, points(8, 2.0, 1).into_iter());
        d.sync_dp();
        assert_eq!(d.dp_refit_calls, 2, "a new batch must refit again");
        assert_same_curve(
            d.app.live_dp_recording().expect("batch 1 refit").curve(),
            d.app.recording().curve(),
        );

        // Turning the pane off drops the recording so nothing can read a
        // recording the pane is no longer showing.
        press(&mut d.app, "d");
        assert!(!d.app.dp_pane());
        assert!(d.app.live_dp_recording().is_none());
        d.sync_dp();
        assert_eq!(
            d.dp_refit_calls, 2,
            "sync_dp must not refit for a pane that is off"
        );

        // Reopening within the same batch: the cache still names batch 1, but
        // the recording it named was dropped on the way out, so serving the
        // cache would leave the pane on with nothing in it.
        press(&mut d.app, "d");
        d.sync_dp();
        assert_eq!(
            d.dp_refit_calls, 3,
            "a reopened pane must refit rather than trust a cache whose \
             recording was dropped"
        );
        assert!(
            d.app.live_dp_recording().is_some(),
            "a reopened pane must have a DP trace to show"
        );
    }

    /// A user who pressed `r`/`q`/Ctrl-C at a Phase 1 pause has already asked to
    /// stop seeing the dashboard — see `present`'s doc comment for why an
    /// unguarded post-Phase-2 pause would stall the run. `render_pause_calls` is
    /// what makes "did not even attempt to render" provable: under `cargo test`
    /// `render_pause` detaches immediately whether or not the guard fired, so
    /// the return value cannot tell a bypassed guard from a working one.
    ///
    /// The `Stay`/`Next` rows are the baseline that keeps an inverted or
    /// always-firing guard from passing: `present` must still open the pause
    /// when nothing has asked it not to.
    #[test]
    fn present_does_not_render_once_a_pause_chose_to_stop() {
        // (action applied at the last pause, expected render_pause calls)
        let cases = [
            (PauseAction::Stay, 1),
            (PauseAction::Next(5), 1),
            (PauseAction::RunToEnd, 0),
            (PauseAction::Detach, 0),
            (PauseAction::Abort, 0),
        ];
        for (action, expected) in cases {
            let mut d = dash(1, 4, 10, 3);
            d.stepper.apply(action);
            d.present();
            assert_eq!(
                d.render_pause_calls, expected,
                "{action:?} must leave present() with {expected} render_pause call(s)"
            );
        }
    }

    /// An overlong batch (more points than `n_calibrants`) must not let the live
    /// re-fit see points the frame slab never recorded — `FrameStore` truncates
    /// its own copy, and `on_batch` must clamp identically so a later
    /// `refit_frame` reproduces what was actually fit live.
    #[test]
    fn an_overlong_batch_is_clamped_the_same_way_the_frame_slab_clamps_it() {
        let n_calibrants = 4;
        let mut d = dash(1, n_calibrants, 10, 3);
        d.on_batch(0, points(2 * n_calibrants, 1.0, 0).into_iter());
        d.finish();

        let (_chunk, pts) = d.frames.frame(0).expect("frame 0 exists");
        assert_eq!(
            pts.len(),
            n_calibrants,
            "the slab truncates to n_calibrants"
        );

        let live = d.app.recording().curve().to_vec();
        let refit = d.refit_frame(0).expect("frame 0 exists").curve().to_vec();
        assert_eq!(
            live.len(),
            n_calibrants,
            "the live fit saw exactly the points the slab recorded, not the full \
             overlong batch"
        );
        assert_same_curve(&live, &refit);
    }

    // ---- terminal entry ----

    #[test]
    fn render_pause_detaches_when_stdout_is_not_a_terminal() {
        let mut d = dash(1, 4, 10, 3);
        assert_eq!(render_pause(&mut d), PauseAction::Detach);
    }

    /// `catch_panics` is what keeps a panic inside the event loop from unwinding
    /// into the search. `#[cfg(debug_assertions)]` because catching only works
    /// under `panic = "unwind"` (`profile.dev`/`test`'s default): a CI leg
    /// running tests against a release-profile build has `panic = "abort"` set
    /// workspace-wide and would abort the whole test process below.
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

    // ---- no-reallocation, pinned by pointer identity ----
    //
    // The convention this codebase uses for every "does not reallocate" claim
    // (`frames.rs`'s `the_slab_is_allocated_once`, `calibrt`'s `filtered_ptr`):
    // capacity equality alone cannot distinguish "reused the same allocation"
    // from "reallocated but landed on the same capacity anyway".

    #[test]
    fn the_point_buffers_never_reallocate_once_constructed() {
        let n_calibrants = 8;
        let mut d = dash(5, n_calibrants, 10, 5);
        // `current_points` is pre-sized at construction, so it is stable from
        // before the very first batch. `prev_points` starts as an empty `Vec`
        // (no allocation at all), so its first-ever growth happens inside
        // batch 0 and its pointer is only meaningful to compare afterwards.
        let current_ptr = d.current_points.as_ptr();
        d.on_batch(0, points(n_calibrants, 1.0, 0).into_iter());
        assert_eq!(
            d.current_points.as_ptr(),
            current_ptr,
            "current_points must not reallocate on batch 0"
        );
        let prev_ptr = d.prev_points.as_ptr();

        for chunk in 1..5 {
            d.on_batch(chunk, points(n_calibrants, 1.0, chunk).into_iter());
            assert_eq!(
                d.current_points.as_ptr(),
                current_ptr,
                "current_points must not reallocate on batch {chunk}"
            );
            assert_eq!(
                d.prev_points.as_ptr(),
                prev_ptr,
                "prev_points must not reallocate again after its first growth (batch {chunk})"
            );
        }
    }

    // ---- prev_curve survives a failed batch ----

    #[test]
    fn a_failed_fit_does_not_reset_the_delta_baseline() {
        // A batch whose fit fails (every weight below `suppress_nonmax`'s 1.0
        // seed, so nothing survives suppression and `curve()` stays `None`) must
        // not wipe out the last real curve as the comparison point for the *next*
        // successful batch's `curve_delta` — otherwise the Convergence tab would
        // show a NaN-to-real discontinuity that never actually happened.
        let n = 4;
        let mut d = dash(3, n, 10, 3);
        // Well above the suppression seed, versus below it for every point.
        let good_points = || scored_points(n, 1.0, 0, |i| 2.0 + i as f64);
        let failing_points = || scored_points(n, 1.0, 0, |_| 0.1);

        d.on_batch(0, good_points().into_iter());
        assert!(
            d.app.metrics()[0].wrmse.is_finite(),
            "batch 0 must produce a real curve"
        );

        d.on_batch(1, failing_points().into_iter());
        assert!(
            d.app.metrics()[1].wrmse.is_nan(),
            "batch 1's fit must fail outright (every weight is sub-threshold)"
        );

        // Identical points to batch 0: if `prev_curve` correctly survived
        // batch 1's failure, this fits the same curve again and the delta
        // against it is ~0 — not NaN, which is what a reset-to-None baseline
        // would report instead.
        d.on_batch(2, good_points().into_iter());
        let m2 = d.app.metrics()[2];
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
