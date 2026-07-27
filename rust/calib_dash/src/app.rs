//! Dashboard app state: the tab/stage cursors, the vim-style count prefix,
//! and the pause action a keypress resolves to.
//!
//! Counts use the vim prefix idiom so there is no modal input: `15n` advances
//! fifteen batches. Without it, "skip N" would need a second text-entry mode
//! with its own key routing, cursor and escape handling, for one integer.
//! The accumulator is a single `Option<u32>` on `App` and serves every
//! motion, not just batch stepping.

use crate::FitRecording;
use ratatui::crossterm::event::{
    KeyCode,
    KeyEvent,
    KeyModifiers,
};

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
    /// The pending count prefix, e.g. the `15` typed before `n` in `15n`.
    /// Lives for exactly one motion: a motion consumes it via `take_count`,
    /// and any other key clears it directly, so a half-typed number never
    /// leaks into the next command.
    count: Option<u32>,
}

impl App {
    pub fn new(bins: usize) -> Self {
        Self {
            tab: Tab::Fit,
            stage: Stage::Grid,
            dp_pane: false,
            batch: 0,
            recording: FitRecording::new(bins),
            count: None,
        }
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

    pub fn recording(&self) -> &FitRecording {
        &self.recording
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
            KeyCode::Char('n') if key.modifiers.is_empty() => PauseAction::Next(self.take_count()),
            KeyCode::Char(']') if key.modifiers.is_empty() => {
                // `Stage` has a fixed, tiny state space, so anything past
                // `ALL.len() - 1` steps is wasted work — without this clamp,
                // a large typed count (reachable now that digit overflow no
                // longer panics) would loop synchronously for as long as it
                // takes to freeze the UI.
                for _ in 0..self.max_stage_steps() {
                    self.stage = self.stage.next();
                }
                PauseAction::Stay
            }
            KeyCode::Char('[') if key.modifiers.is_empty() => {
                for _ in 0..self.max_stage_steps() {
                    self.stage = self.stage.prev();
                }
                PauseAction::Stay
            }
            KeyCode::Char('d') if key.modifiers.is_empty() => {
                self.dp_pane = !self.dp_pane;
                self.count = None;
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
}
