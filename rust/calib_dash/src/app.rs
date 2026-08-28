//! Dashboard app state: the tab cursor, the vim-style count prefix, and the
//! pause action a keypress resolves to. [`CalibDash`] ties all of it to a
//! running search, one [`CalibDash::on_batch`] call per Phase 1 scoring chunk.
//!
//! Counts use the vim prefix idiom: `15n` advances fifteen batches.

use crate::frames::{
    CalibrantPoint,
    FrameStore,
    FrameSummary,
};
use crate::metrics::{
    BatchMetrics,
    churn,
    curve_delta,
};
use crate::recording::FitRecording;
use calibrt::{
    CalibrationCurve,
    CalibrationState,
    LibraryRT,
    ObservedRTSeconds,
    RidgeSummary,
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
pub(crate) enum Tab {
    Fit,
    Convergence,
    Tolerances,
}

impl Tab {
    pub(crate) const ALL: [Tab; 3] = [Tab::Fit, Tab::Convergence, Tab::Tolerances];

    pub(crate) fn title(self) -> &'static str {
        match self {
            Tab::Fit => "Fit",
            Tab::Convergence => "Convergence",
            Tab::Tolerances => "Tolerances",
        }
    }
}

/// Which mark overlay the Fit heatmap draws on top of the density field,
/// cycled by `m`/`M`. Exactly one is active at a time.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Layer {
    /// Bare heatmap — the raw density the fit works from.
    None,
    /// The DP chain and the greedily attached tails, glyph-distinguished within the
    /// layer (`O`/`X`): is a bad edge the DP's choice or a tail grafted on after?
    Path,
    Curve,
    Ridge,
    Suppressed,
}

impl Layer {
    pub(crate) const ALL: [Layer; 5] = [
        Layer::None,
        Layer::Path,
        Layer::Curve,
        Layer::Ridge,
        Layer::Suppressed,
    ];

    pub(crate) fn label(self) -> &'static str {
        match self {
            Layer::None => "none",
            Layer::Path => "path",
            Layer::Curve => "curve",
            Layer::Ridge => "ridge",
            Layer::Suppressed => "suppressed",
        }
    }
}

/// Steps `cur` through `all` by `steps` stops, wrapping. `Tab` and `Layer` are both
/// fixed stop lists with no natural end to clamp at, so `l` past the last tab is the
/// first again.
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

/// Whether the batch loop driving Phase 1 scoring should keep going after a pause.
/// Only `Ctrl-C` ever produces `Abort` — a dashboard failure anywhere else is
/// logged and treated as `Continue`, since a dev tool must never fail a search.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Flow {
    Continue,
    Abort,
}

/// What the user asked for on leaving a pause.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PauseAction {
    /// Stay in the pause — the key changed the view, not the flow.
    Stay,
    /// Advance this many batches before pausing again.
    Next(u32),
    /// Stop showing the dashboard and let the pipeline run on. Both `r` (run to
    /// end) and `q` (detach) land here — they differ only in what the user meant,
    /// not in what happens.
    Detach,
    Abort,
}

/// Skip bookkeeping, so the render-or-skip decision is testable without a terminal.
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

    /// Consulted by `on_batch` exactly once per chunk: asking the question is
    /// what draws down the skip counter.
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

    /// Whether a prior pause resolved to `Detach` or `Abort` — the user has already
    /// asked to stop seeing the dashboard, and nothing should reopen it.
    fn is_stopped(&self) -> bool {
        self.stopped
    }

    /// Called after a pause with whatever the user chose.
    fn apply(&mut self, action: PauseAction) {
        match action {
            PauseAction::Stay => {}
            // `Next(5)` at the chunk-0 pause means "render at chunk 5": four skipped.
            PauseAction::Next(n) => self.remaining = n.saturating_sub(1),
            PauseAction::Detach | PauseAction::Abort => self.stopped = true,
        }
    }
}

/// Dashboard state for one pause: the tab, the batch number, the recording of the
/// fit that ran at it, and the pending count prefix.
pub(crate) struct App {
    tab: Tab,
    batch: u32,
    recording: FitRecording,
    /// The post-Phase-2 recording, `None` for the whole of Phase 1: the Tolerances
    /// tab reads this to decide whether Step B has run yet at all.
    real_fit: Option<FitRecording>,
    /// The Tolerances tab's post-Step-B m/z/mobility/RT summary, set alongside
    /// `real_fit`.
    tolerances: Option<ToleranceSummary>,
    /// Every batch's metrics, undecimated.
    metrics: Vec<BatchMetrics>,
    /// What the Convergence tab's header states about the frame slab's decimation.
    /// Defaults to "nothing decimated", refreshed through `set_frame_summary`.
    frames: FrameSummary,
    layer: Layer,
    /// The pending count prefix, e.g. the `15` typed before `n` in `15n`. Lives for
    /// exactly one motion: a motion consumes it via `take_count` and any other key
    /// clears it, so a half-typed number never leaks into the next command.
    count: Option<u32>,
    /// Which retained frame `<`/`>` have scrubbed the Fit tab to, `None` for live.
    /// Bounded by `frames.retained`; the recording lives in `scrub_recording`.
    scrub_frame: Option<usize>,
    /// The original batch/chunk number of `scrub_frame`, for the Fit tab's "not live"
    /// banner — "frame 2 of 5" alone means little to a user who thinks in batches.
    scrub_chunk: Option<usize>,
    /// The scrubbed frame's refit recording. `None` whenever `scrub_frame` is, and
    /// also momentarily after `scrub_frame` moves and before the next `sync_scrub`
    /// catches up — `active_recording` falls back to the live view for that gap.
    scrub_recording: Option<FitRecording>,
    /// Whether the `?` key-map overlay is open. Modal: while it shows, every keypress
    /// (digits and `Ctrl-C` included) only dismisses it.
    show_keys: bool,
}

impl App {
    pub(crate) fn new(bins: usize) -> Self {
        Self {
            tab: Tab::Fit,
            batch: 0,
            recording: FitRecording::new(bins),
            real_fit: None,
            tolerances: None,
            metrics: Vec::new(),
            frames: FrameSummary {
                retained: 0,
                keep_every: 1,
                dropped: 0,
            },
            layer: Layer::None,
            count: None,
            scrub_frame: None,
            scrub_chunk: None,
            scrub_recording: None,
            show_keys: false,
        }
    }

    /// Mutable access to the live recording, so a caller can fit directly into the
    /// allocation `App` already owns.
    pub(crate) fn recording_mut(&mut self) -> &mut FitRecording {
        &mut self.recording
    }

    pub(crate) fn real_fit(&self) -> Option<&FitRecording> {
        self.real_fit.as_ref()
    }

    pub(crate) fn tolerances(&self) -> Option<&ToleranceSummary> {
        self.tolerances.as_ref()
    }

    pub(crate) fn set_final(&mut self, recording: FitRecording, tolerances: ToleranceSummary) {
        self.real_fit = Some(recording);
        self.tolerances = Some(tolerances);
    }

    pub(crate) fn metrics(&self) -> &[BatchMetrics] {
        &self.metrics
    }

    pub(crate) fn push_metrics(&mut self, m: BatchMetrics) {
        self.metrics.push(m);
    }

    pub(crate) fn frames(&self) -> FrameSummary {
        self.frames
    }

    pub(crate) fn set_frame_summary(&mut self, summary: FrameSummary) {
        self.frames = summary;
    }

    pub(crate) fn tab(&self) -> Tab {
        self.tab
    }

    pub(crate) fn batch(&self) -> u32 {
        self.batch
    }

    pub(crate) fn layer(&self) -> Layer {
        self.layer
    }

    pub(crate) fn recording(&self) -> &FitRecording {
        &self.recording
    }

    /// What the Fit tab should actually draw: the scrubbed frame's recording while
    /// `scrub_frame` is set and `sync_scrub` has caught up to it, the live `recording`
    /// otherwise — never a blank tab for the gap between the two.
    pub(crate) fn active_recording(&self) -> &FitRecording {
        if self.scrub_frame.is_some()
            && let Some(rec) = self.scrub_recording.as_ref()
        {
            rec
        } else {
            &self.recording
        }
    }

    pub(crate) fn scrub_frame(&self) -> Option<usize> {
        self.scrub_frame
    }

    pub(crate) fn scrub_chunk(&self) -> Option<usize> {
        self.scrub_chunk
    }

    /// Wires a freshly refit frame into the scrub view (`App` has no `FrameStore`).
    pub(crate) fn set_scrub_recording(
        &mut self,
        frame: usize,
        chunk: usize,
        recording: FitRecording,
    ) {
        self.scrub_frame = Some(frame);
        self.scrub_chunk = Some(chunk);
        self.scrub_recording = Some(recording);
    }

    /// Drops back to the live view — the user scrubbed past the last retained frame,
    /// or a frame could not be refit and live beats a stale or absent recording.
    pub(crate) fn clear_scrub(&mut self) {
        self.scrub_frame = None;
        self.scrub_chunk = None;
        self.scrub_recording = None;
    }

    pub(crate) fn pending_count(&self) -> Option<u32> {
        self.count
    }

    pub(crate) fn show_keys(&self) -> bool {
        self.show_keys
    }

    /// Consumes the pending count, defaulting to 1 for a bare motion.
    fn take_count(&mut self) -> u32 {
        self.count.take().unwrap_or(1)
    }

    /// Folds one more digit into the pending count. A leading `0` does not start a
    /// count — matching vim, which keeps `0` free as a motion — but extends one.
    fn push_digit(&mut self, d: u32) {
        if d == 0 && self.count.is_none() {
            return;
        }
        self.count = Some(self.count.unwrap_or(0).saturating_mul(10).saturating_add(d));
    }

    /// Moves the Fit tab's frame-scrub cursor back by `n` retained frames, or does
    /// nothing if none have been retained. Saturating, so any `n` is safe.
    fn scrub_back(&mut self, n: u32) {
        let total = self.frames.retained as u32;
        if total == 0 {
            return;
        }
        let pos = self.scrub_frame.map(|i| i as u32).unwrap_or(total);
        self.scrub_frame = Some(pos.saturating_sub(n) as usize);
    }

    /// Moves the cursor forward by `n`. Past the last retained frame this returns to
    /// the live view rather than clamping — "keep pressing `>` and you get to now".
    fn scrub_forward(&mut self, n: u32) {
        let Some(i) = self.scrub_frame else {
            return;
        };
        let total = self.frames.retained as u32;
        let new_pos = (i as u32).saturating_add(n);
        self.scrub_frame = if new_pos >= total {
            None
        } else {
            Some(new_pos as usize)
        };
    }

    /// Routes one keypress. Every bare-character binding lives in
    /// `handle_plain_char`, behind a single `modifiers.is_empty()` guard, so a
    /// modified key can never fall into a bare-letter arm. The two shifted
    /// characters below are the exception, for the reason given at their arm.
    pub(crate) fn handle_key(&mut self, key: KeyEvent) -> PauseAction {
        // The `?` overlay is modal: any keypress only dismisses it, so `n` to dismiss
        // does not also skip a batch. Deliberately ahead of the `Ctrl-C` check: while
        // the overlay is up, dismissing it is the whole of what *any* key does.
        if self.show_keys {
            self.show_keys = false;
            return PauseAction::Stay;
        }
        if key.code == KeyCode::Char('c') && key.modifiers.contains(KeyModifiers::CONTROL) {
            self.count = None;
            return PauseAction::Abort;
        }

        match key.code {
            // `?` and `M` are shifted characters and terminals disagree on whether
            // they report `Shift` on the resulting `Char`, so these two guard only
            // against `Control`. They precede the plain-character arm.
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
            // Esc lands here, and needs no modifier guard: it clears the count and
            // stays either way, exactly what every other unbound key does.
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
                // `n` advances the *live* batch loop, so a scrub cursor left over from
                // browsing history must not still be showing at the next pause.
                self.clear_scrub();
                PauseAction::Next(n)
            }
            'l' => {
                self.tab = cycle(&Tab::ALL, self.tab, self.take_count(), true);
                PauseAction::Stay
            }
            'h' => {
                self.tab = cycle(&Tab::ALL, self.tab, self.take_count(), false);
                PauseAction::Stay
            }
            // `m` for "mark"; shifted `M` (in `handle_key`) cycles back, as `h` does.
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
            'r' | 'q' => {
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

/// How many library-RT samples `curve_delta` draws between batches' curves. Any
/// `>= 2` is defined; this is high enough that a local kink is unlikely to hide.
const CURVE_DELTA_SAMPLES: usize = 50;

/// [`calibrt::CalibrationState::refit`], replacing `recording` with what the fit
/// left behind, against whichever state/recording pair the caller owns. The live
/// fit and a scrubbed frame's replay sharing this is what makes a replayed frame
/// show the fit that actually ran.
///
/// `None` when `calibrt` refuses, logged at `warn` here: every caller treats a
/// refusal the same way, as "no recording".
fn fit_points(
    state: &mut CalibrationState,
    recording: &mut FitRecording,
    bins: usize,
    points: &[CalibrantPoint],
) -> Option<(f64, f64)> {
    let (x_range, _) = state
        .refit(bins, points.iter().map(|p| (p.library_rt, p.observed_rt)))
        // `CalibRtError` has no `Display`.
        .inspect_err(|e| tracing::warn!("calib_dash: skipping this re-fit: {e:?}"))
        .ok()?;
    *recording = FitRecording::from_state(state);
    Some(x_range)
}

/// `CalibrantPoint` as `CalibrationCurve::wrmse` wants it, at the same
/// [`calibrt::CALIBRANT_WEIGHT`] the fit itself used.
fn as_calibrt_tuple(p: &CalibrantPoint) -> (LibraryRT<f64>, ObservedRTSeconds<f64>, f64) {
    (
        LibraryRT(p.library_rt),
        ObservedRTSeconds(p.observed_rt),
        calibrt::CALIBRANT_WEIGHT,
    )
}

/// Summary the Tolerances tab renders once Step B has run: the derived m/z, mobility
/// and RT tolerances plus how many calibrants they came from.
///
/// `mz_ppm` and `mobility_pct` are signed `(lo, hi)` windows around zero and are
/// genuinely asymmetric — Step B can derive a different tolerance on either side of a
/// calibrant. `rt_seconds` is symmetric by construction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ToleranceSummary {
    pub mz_ppm: (f64, f64),
    pub mobility_pct: (f64, f64),
    /// Half-width: the tolerance is `±rt_seconds`.
    ///
    /// A *fallback*: Phase 3 prefers the ridge half-width interpolated at each
    /// query's own library RT (`CalibrationResult::get_tolerance`) and only falls
    /// back to this MAD-derived scalar where the ridge has no measurement.
    pub rt_seconds: f64,
    pub n_calibrants: usize,
}

/// Ties the dashboard to a running Phase 1 batch loop, through `on_batch` once per
/// scoring chunk.
///
/// Every buffer here is allocated at startup and reused, so steady state adds nothing
/// per batch beyond the `CalibrationCurve` `calibrt` allocates for every fit anyway
/// (plus the `prev_curve` clone `curve_delta` compares against).
pub struct CalibDash {
    app: App,
    frames: FrameStore,
    stepper: Stepper,
    /// The live re-fit's `CalibrationState`, reused every batch.
    state: CalibrationState,
    /// This batch's points, buffered once so recording, re-fitting and churn all read
    /// the same materialized slice. Clamped to `n_calibrants` exactly the way
    /// `FrameStore::record` clamps its own copy: without matching clamps, an overlong
    /// batch would fit the live curve against points the frame slab never recorded, so
    /// replaying that frame would disagree with what was shown live.
    current_points: Vec<CalibrantPoint>,
    n_calibrants: usize,
    /// The previous batch's points, for `churn`. Empty before the first batch —
    /// `churn`'s own "everything admitted, nothing evicted" case.
    prev_points: Vec<CalibrantPoint>,
    /// The previous batch's fitted curve, for `curve_delta`. `None` before the first
    /// *successful* fit: a failed batch must not reset this, or the next success would
    /// report a NaN delta as though there had never been a curve at all.
    prev_curve: Option<CalibrationCurve>,
    /// Held so `refit_frame` can reconfigure `refit_state` to the live grid size.
    bins: usize,
    /// A separate `CalibrationState`/`FitRecording` pair for a scrubbed frame's
    /// replay, so replaying does not perturb the live recording.
    refit_state: CalibrationState,
    refit_recording: FitRecording,
}

impl CalibDash {
    pub fn new(
        n_frames: usize,
        n_calibrants: usize,
        bins: usize,
        lookback: usize,
        budget_bytes: usize,
    ) -> Self {
        // The unit range is only ever a placeholder: the first `on_batch` reconfigures
        // every state to the batch's own point-derived ranges before any fit.
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
        }
    }

    /// In order: record the frame, re-fit the live curve from those same points, push
    /// this batch's metrics (every batch, rendered or not — see `metrics.rs`), then
    /// either continue or pause to render.
    ///
    /// `Flow::Abort` happens only because the user pressed Ctrl-C at a pause: a
    /// dashboard that cannot fit or cannot draw never fails a search.
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

    /// Opens the post-Phase-2 pause, so the Tolerances tab is reachable once Step B
    /// runs. There is no batch loop left to steer, so every `PauseAction` collapses to
    /// "return to the caller", `Abort` included.
    ///
    /// Guarded on `stepper.is_stopped()`: a user who chose `r`/`q`/Ctrl-C at any Phase
    /// 1 pause has already asked to stop seeing the dashboard, and the caller reaches
    /// this call regardless of `Flow::Abort`. Unguarded, it would reopen the alternate
    /// screen and block on a keypress that never comes.
    pub fn present(&mut self) {
        if self.stepper.is_stopped() {
            return;
        }
        render_pause(self);
    }

    /// Promotes the reserved final frame and refreshes the Convergence header's
    /// decimation numbers one last time — the promotion can change them.
    pub fn finish(&mut self) {
        self.frames.finish();
        self.sync_frame_summary();
    }

    /// Wires the post-Phase-2 recording and derived tolerances in once Step B ran.
    pub fn show_final(&mut self, recording: FitRecording, tolerances: ToleranceSummary) {
        self.app.set_final(recording, tolerances);
    }

    /// Re-fits a recorded frame's points from scratch, into `refit_state`/
    /// `refit_recording` rather than the live pair, through the same `fit_points`
    /// sequence `refit_live` runs — so this reproduces the batch that actually ran.
    /// `None` if the frame index doesn't exist or the re-fit was skipped.
    fn refit_frame(&mut self, i: usize) -> Option<(usize, &FitRecording)> {
        let (chunk, pts) = self.frames.frame(i)?;
        fit_points(
            &mut self.refit_state,
            &mut self.refit_recording,
            self.bins,
            pts,
        )?;
        Some((chunk, &self.refit_recording))
    }

    /// Re-fits the live `CalibrationState` from `current_points` and returns this
    /// batch's `BatchMetrics`. The two live together because the metrics *are* a
    /// comparison of this re-fit against the last one.
    fn refit_live(&mut self, chunk: usize) -> BatchMetrics {
        let n_points = self.current_points.len();
        let mut wrmse = f64::NAN;
        let mut max_delta = f64::NAN;
        let mut mean_delta = f64::NAN;
        let mut path_nodes = 0usize;
        let mut ridge_half_width = f64::NAN;

        if let Some(x_range) = fit_points(
            &mut self.state,
            self.app.recording_mut(),
            self.bins,
            &self.current_points,
        ) {
            path_nodes = self.app.recording().path_indices().len();
            ridge_half_width = RidgeSummary::of(self.app.recording().ridge())
                .map_or(f64::NAN, |s| s.weighted_half_width);

            if let Some(curve) = self.state.curve() {
                wrmse = curve.wrmse(self.current_points.iter().map(as_calibrt_tuple));
                if let Some(prev) = &self.prev_curve {
                    let (d_max, d_mean) = curve_delta(prev, curve, x_range, CURVE_DELTA_SAMPLES);
                    max_delta = d_max;
                    mean_delta = d_mean;
                }
                self.prev_curve = Some(curve.clone());
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

    /// Pushes `frames`' real decimation numbers into `app`.
    fn sync_frame_summary(&mut self) {
        self.app.set_frame_summary(self.frames.summary());
    }

    /// Keeps the Fit tab's scrub view in lockstep with `App::scrub_frame`. `App` moves
    /// `scrub_frame` on `<`/`>` but has no `FrameStore` to refit from, so this is the
    /// other half of that key handling. A no-op when `scrub_frame` is `None`.
    ///
    /// Otherwise re-fits and clones unconditionally: this runs at most once per
    /// keypress inside an interactive pause, never per batch, so the clone is off
    /// `on_batch`'s steady-state path entirely.
    fn sync_scrub(&mut self) {
        let Some(i) = self.app.scrub_frame() else {
            return;
        };
        match self.refit_frame(i) {
            Some((chunk, rec)) => {
                let cloned = rec.clone();
                self.app.set_scrub_recording(i, chunk, cloned);
            }
            // The index is bounded by `frames.retained`, set from this same
            // `FrameStore`, so this only fires on a frame too degenerate to refit.
            None => self.app.clear_scrub(),
        }
    }
}

/// Pauses the batch loop to render one interactive frame and block until the user's
/// keypress resolves to something other than [`PauseAction::Stay`]. Takes the whole
/// [`CalibDash`] because the batch scrubber needs `refit_frame`'s `FrameStore` access
/// on every redraw.
///
/// Never fails the caller: with stdout not a terminal (under `cargo test`, or a piped
/// run) this returns [`PauseAction::Detach`] immediately, and a terminal that fails to
/// initialize is logged at `warn` and treated the same way.
fn render_pause(dash: &mut CalibDash) -> PauseAction {
    if !std::io::stdout().is_terminal() {
        return PauseAction::Detach;
    }

    let mut terminal = match ratatui::try_init() {
        Ok(terminal) => terminal,
        Err(e) => {
            // `try_init` can fail after partially succeeding (raw mode enabled, then
            // entering the alternate screen fails), so restore unconditionally.
            tracing::warn!("calib_dash failed to initialize the terminal: {e}");
            ratatui::restore();
            return PauseAction::Detach;
        }
    };
    let action = event_loop(&mut terminal, dash);
    ratatui::restore();
    action
}

/// Draws and reads key events in the alternate screen until a keypress resolves to
/// something other than `Stay`. A draw or read failure is logged and detached, never
/// propagated to fail the batch loop.
fn event_loop<B: ratatui::backend::Backend>(
    terminal: &mut ratatui::Terminal<B>,
    dash: &mut CalibDash,
) -> PauseAction {
    loop {
        // Before every draw, not just after a scrub keypress: a no-op when `None`.
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
        // Only key *presses*: on Windows crossterm also emits releases, which would
        // double every keystroke.
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

    /// A dashboard with a budget generous enough to retain every frame.
    fn dash(n_frames: usize, n_calibrants: usize, bins: usize, lookback: usize) -> CalibDash {
        CalibDash::new(n_frames, n_calibrants, bins, lookback, 1 << 20)
    }

    /// `n` points on the line `observed_rt = library_rt * slope`. A negative slope is
    /// how a test asks for a fit that fails: the DP only chains pairs that increase in
    /// both axes, so an anti-correlated batch yields a one-node path and no curve.
    fn points(n: usize, slope: f64, chunk: usize) -> Vec<CalibrantPoint> {
        (0..n)
            .map(|i| CalibrantPoint {
                library_rt: i as f64 + 0.5,
                observed_rt: (i as f64 + 0.5) * slope,
                library_id: (chunk * n + i) as u64,
            })
            .collect()
    }

    /// Points at explicit `library_id`es, so a test can choose the two sets `churn`
    /// diffs. On the identity line, so the fit succeeds.
    fn indexed_points(indices: &[usize]) -> Vec<CalibrantPoint> {
        indices
            .iter()
            .enumerate()
            .map(|(i, &idx)| CalibrantPoint {
                library_rt: i as f64 + 0.5,
                observed_rt: i as f64 + 0.5,
                library_id: idx as u64,
            })
            .collect()
    }

    /// `n` retained frames, nothing decimated — what the scrubber tests need.
    fn retained(n: usize) -> FrameSummary {
        FrameSummary {
            retained: n,
            keep_every: 1,
            dropped: 0,
        }
    }

    /// A recording's curve points, empty when it holds no curve.
    fn curve_of(rec: &FitRecording) -> Vec<Point> {
        rec.curve().map_or_else(Vec::new, |c| c.points().to_vec())
    }

    /// Pointwise curve equality — behind every "what you scrub is what ran" test.
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

    /// `0` is the one digit whose meaning depends on what came before it: it cannot
    /// *start* a count (vim keeps `0` free as a motion) but extends one.
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
        // 16 nines: `9_999_999_999 * 10 + 9` overflows a u32 well before the end of
        // this run, so this must not panic (dev/test builds check overflow).
        press(&mut app, "9999999999999999");
        assert_eq!(app.pending_count(), Some(u32::MAX));
    }

    /// Both halves of "only a motion consumes a count": a bound non-motion (`r`) and
    /// an unbound key (`Esc`, `handle_key`'s catch-all arm) each clear it.
    #[test]
    fn a_non_motion_key_clears_a_pending_count() {
        let mut app = App::new(10);
        press(&mut app, "42r");
        assert_eq!(
            app.pending_count(),
            None,
            "`r` detaches, it is not a motion"
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
    fn both_r_and_q_detach() {
        let mut app = App::new(10);
        assert_eq!(press(&mut app, "r"), PauseAction::Detach);
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
            "`?` is a non-motion key, like `r`/`q`, so it clears a pending count"
        );
    }

    #[test]
    fn any_key_dismisses_the_overlay_without_also_acting_on_it() {
        let mut app = App::new(10);
        app.handle_key(key('?'));
        assert!(app.show_keys());

        let action = app.handle_key(key('n'));
        assert!(!app.show_keys(), "the overlay must close on the next key");
        assert_eq!(
            action,
            PauseAction::Stay,
            "the dismissing key must not also perform its usual action"
        );

        // A digit takes a different route through `handle_plain_char` than a motion,
        // so it needs its own case.
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

        // A count large enough to overflow the `+ n` must saturate rather than wrap
        // into some other index.
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

    /// `cycle` reduces the count mod the number of stops before anything moves, so a
    /// ten-digit count returns promptly rather than driving ~10^9 iterations on the
    /// input thread. Both counts below land two stops from the start
    /// (`1_000_000_001 % 3 == 2`, `1_000_000_002 % 5 == 2`).
    ///
    /// The backward keys make the reduction load-bearing rather than merely fast:
    /// `cycle` walks back by `n - steps`, which underflows for any count past the stop
    /// list's length — and a count past 5 is an ordinary keystroke.
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
        let rendered: Vec<usize> = (1..=10).filter(|_| s.should_render()).collect();
        assert_eq!(
            rendered,
            vec![5, 6, 7, 8, 9, 10],
            "chunks 1-4 skipped, render at 5, then every batch again"
        );
    }

    #[test]
    fn a_stopping_action_never_renders_again() {
        for action in [PauseAction::Detach, PauseAction::Abort] {
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
            // short-circuit rather than wait for a keypress.
            let flow = d.on_batch(chunk, points(8, 1.0, chunk).into_iter());
            assert!(matches!(flow, Flow::Continue), "chunk {chunk}");
        }
        d.finish();
        assert_eq!(d.app.metrics().len(), 10, "metrics are undecimated");
        assert!(d.frames.summary().retained < 10, "frames are decimated");
    }

    /// `BatchMetrics.admitted`/`evicted` are the only two fields nothing but
    /// `refit_live` ever fills, so the wiring needs its own test.
    ///
    /// The second batch is deliberately asymmetric (one index in, two out), so
    /// swapping the two fields at the construction site cannot pass. It also fails if
    /// the `prev_points` refresh moves above the `churn` call, which would diff
    /// `current_points` against itself and report `(0, 0)` forever.
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

    /// End-to-end version of the `<` scrubber: moving `App::scrub_frame` and letting
    /// `sync_scrub` catch up must leave the Fit tab on the earlier batch's real curve.
    #[test]
    fn sync_scrub_reproduces_an_earlier_batchs_curve_on_the_fit_tab() {
        let mut d = dash(3, 8, 10, 5);
        for chunk in 0..3 {
            let slope = 1.0 + chunk as f64 * 0.1;
            d.on_batch(chunk, points(8, slope, chunk).into_iter());
        }
        d.finish();

        // Batch 0's curve, refit independently for comparison.
        let batch0_curve = curve_of(d.refit_frame(0).expect("frame 0 exists").1);

        // A generous budget retains all 3 batches; scrub 3 back from live
        // (batch 2) lands on frame index 0 (batch 0).
        press(&mut d.app, "3<");
        d.sync_scrub();

        assert_eq!(d.app.scrub_frame(), Some(0));
        assert_eq!(d.app.scrub_chunk(), Some(0));
        assert_same_curve(&curve_of(d.app.active_recording()), &batch0_curve);

        // `recording()`, unlike `active_recording()`, still reflects batch 2's fit.
        assert!(
            !curve_of(d.app.recording()).is_empty(),
            "sanity: batch 2 produced a real curve too"
        );
    }

    /// An overlong batch (more points than `n_calibrants`) must not let the live re-fit
    /// see points the frame slab never recorded — `FrameStore` truncates its own copy,
    /// and `on_batch` must clamp identically.
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

        let live = curve_of(d.app.recording());
        let refit = curve_of(d.refit_frame(0).expect("frame 0 exists").1);
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

    // ---- prev_curve survives a failed batch ----

    #[test]
    fn a_failed_fit_does_not_reset_the_delta_baseline() {
        // Anti-correlated calibrants fail the fit: no monotonic chain, no curve.
        let n = 4;
        let mut d = dash(3, n, 10, 3);
        let good_points = || points(n, 1.0, 0);
        let failing_points = || points(n, -1.0, 0);

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

        // Identical points to batch 0: if `prev_curve` survived batch 1's failure this
        // fits the same curve again, and the delta against it is ~0.
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
