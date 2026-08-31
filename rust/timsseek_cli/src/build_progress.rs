//! Rendering for msspeculator's build progress callback.
//!
//! A proteome build runs for minutes in three phases that count different
//! things, so the callback hands over a [`Phase`] with every update rather than
//! a bare fraction. All this module decides is where the update goes.
//!
//! Two renderings, chosen once from whether stderr is a terminal, because they
//! are not the same thing narrowed: an in-place bar written to a redirected
//! stream is one enormous line, and a line every few seconds on a terminal
//! scrolls the bar away. The same split `make_progress_bar` already makes for
//! the search phases, so a build looks like the rest of the tool.

use std::io::IsTerminal;
use std::sync::Mutex;
use std::time::{
    Duration,
    Instant,
};

use indicatif::ProgressBar;
use msspeculator_inference::{
    Exactness,
    Phase,
    Progress,
};
use timsseek::scoring::timings::make_progress_bar;
use tracing::info;

/// Line rate for the logged rendering. A ten-minute build leaves twenty lines.
const LOG_INTERVAL: Duration = Duration::from_secs(30);

/// The phase currently being rendered, and whatever that rendering needs to
/// carry between updates.
enum Current {
    /// A terminal: one bar per phase, replaced when the phase changes.
    Bar(ProgressBar),
    /// Not a terminal: when a line last went out, so the next one can wait.
    /// `None` until the first one, which is never made to wait -- a phase whose
    /// only update is also its last (a small FASTA digests in one read) would
    /// otherwise report nothing at all.
    Logged(Option<Instant>),
}

struct State {
    phase: Phase,
    current: Current,
}

/// Renders [`Progress`] updates for one build.
///
/// Reports through `&self`: the callback is borrowed by the build, which reports
/// from its producer thread, so the interior mutability has to live here. The
/// `Mutex` is uncontended in practice for that reason -- it is what makes the
/// closure `Sync`, not a concurrency mechanism.
pub struct BuildProgress {
    state: Mutex<Option<State>>,
    terminal: bool,
}

impl BuildProgress {
    pub fn new() -> Self {
        Self {
            state: Mutex::new(None),
            terminal: std::io::stderr().is_terminal(),
        }
    }

    /// Hand this to `StreamOptions::progress`.
    pub fn callback(&self) -> impl Fn(Progress) + Send + Sync + '_ {
        move |progress| self.report(progress)
    }

    fn report(&self, progress: Progress) {
        let mut state = self.state.lock().expect("progress state is never poisoned");
        // A phase change is itself the news that the previous phase finished --
        // the callback contract says phases do not interleave -- so the previous
        // rendering is closed here rather than waiting for a final update that
        // `Loading` never sends.
        if state.as_ref().is_none_or(|s| s.phase != progress.phase) {
            if let Some(State {
                current: Current::Bar(bar),
                ..
            }) = state.take()
            {
                close(bar);
            }
            *state = Some(State {
                phase: progress.phase,
                current: self.begin(progress),
            });
        }
        let Some(state) = state.as_mut() else {
            return;
        };
        match &mut state.current {
            Current::Bar(bar) => bar.set_position(progress.done),
            Current::Logged(last) => {
                // `done >= total` is the closing update; letting the interval
                // swallow it would leave the last line reading part-done.
                let closing = progress.total > 0 && progress.done >= progress.total;
                let due = last.is_none_or(|last| last.elapsed() >= LOG_INTERVAL);
                if closing || due {
                    *last = Some(Instant::now());
                    log_line(progress);
                }
            }
        }
    }

    /// Choose this phase's rendering. Emits nothing: the update that opened the
    /// phase then goes through the same path as every later one, so a phase
    /// cannot report its first update twice.
    ///
    /// A phase that measures nothing gets no bar even on a terminal --
    /// `make_progress_bar` would draw `0/0` and an ETA it cannot have. That is
    /// the phase which exists to explain a pause, and one line does that.
    fn begin(&self, progress: Progress) -> Current {
        if !self.terminal || progress.total == 0 {
            return Current::Logged(None);
        }
        Current::Bar(make_progress_bar(
            progress.total,
            &capitalized(progress.phase),
        ))
    }

    /// Close the phase still being rendered.
    ///
    /// Has to be called before anything else writes to the terminal: a bar left
    /// open is a line another writer appends to, so the build's summary lands on
    /// the end of the progress bar. Idempotent, and also runs on drop, which is
    /// what covers the error path -- `write_library` returning `Err` skips
    /// everything after it, and an error message is exactly the thing worst
    /// spliced onto a half-drawn bar.
    pub fn finish(&self) {
        let taken = self
            .state
            .lock()
            .expect("progress state is never poisoned")
            .take();
        if let Some(State {
            current: Current::Bar(bar),
            ..
        }) = taken
        {
            close(bar);
        }
    }
}

impl Drop for BuildProgress {
    fn drop(&mut self) {
        self.finish();
    }
}

impl Default for BuildProgress {
    fn default() -> Self {
        Self::new()
    }
}

/// Take a finished bar off the screen.
///
/// Cleared rather than left, unlike the search phases: `finish` leaves the
/// cursor on the bar's own line, so the next thing written -- the next phase's
/// line, the build summary, an error -- continues it. What a completed phase
/// counted and how long it took is in the `.config.json` sidecar either way.
fn close(bar: ProgressBar) {
    bar.finish_and_clear();
}

/// The phase label as a sentence opener: "Digesting", not "digesting".
fn capitalized(phase: Phase) -> String {
    let label = phase.label();
    let mut out = label.to_uppercase();
    out.truncate(1);
    out.push_str(&label[1..]);
    out
}

/// One line for a log, hedged where the phase says its own numbers are.
fn log_line(progress: Progress) {
    let phase = progress.phase;
    if progress.total == 0 {
        info!("{}", capitalized(phase));
        return;
    }
    let about = match phase.exactness() {
        Exactness::Exact => "",
        Exactness::Approximate => "about ",
    };
    info!(
        "{}: {about}{}/{} {} ({:.0}%)",
        capitalized(phase),
        progress.done,
        progress.total,
        phase.unit(),
        progress.fraction() * 100.0,
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    fn at(phase: Phase, done: u64, total: u64) -> Progress {
        Progress { phase, done, total }
    }

    /// Every phase's label reads as a sentence opener, including one this
    /// version has never seen -- the point of `Phase::label` is that a renderer
    /// does not enumerate phases.
    #[test]
    fn phase_labels_are_capitalized_without_being_enumerated() {
        assert_eq!(capitalized(Phase::Digesting), "Digesting");
        assert_eq!(capitalized(Phase::Loading), "Loading");
        assert_eq!(capitalized(Phase::Predicting), "Predicting");
    }

    /// Reporting is `&self` because the build borrows the callback, so the
    /// renderer has to be usable from the closure it hands over. If this stops
    /// compiling the callback cannot be passed at all.
    #[test]
    fn the_callback_can_be_handed_to_a_build() {
        let progress = BuildProgress::new();
        let callback = progress.callback();
        fn takes_progress_fn(_: &(dyn Fn(Progress) + Send + Sync + '_)) {}
        takes_progress_fn(&callback);

        // A whole build's shape, in order, including the phase that measures
        // nothing and the one that restarts `done` at zero.
        for update in [
            at(Phase::Digesting, 0, 4),
            at(Phase::Digesting, 4, 4),
            at(Phase::Loading, 0, 0),
            at(Phase::Predicting, 0, 9),
            at(Phase::Predicting, 9, 9),
        ] {
            callback(update);
        }
    }
}
