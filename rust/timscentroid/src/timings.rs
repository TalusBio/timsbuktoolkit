//! Progressive CLI timing output.
//!
//! [`TimedStep`] prints a dot-padded label immediately, opens a tracing span,
//! then appends elapsed time when the work finishes.

use std::fmt;
use std::time::{
    Duration,
    Instant,
};

/// A timed step that prints a dot-padded label immediately, opens a tracing
/// span, and appends elapsed time on finish.
///
/// ```ignore
/// let step = TimedStep::begin("Loading speclib");
/// let speclib = load_speclib()?;
/// let elapsed = step.finish_with(format_args!("{} entries", speclib.len()));
/// // terminal: "Loading speclib .......... 834.567ms (225178 entries)"
/// ```
pub struct TimedStep {
    start: Instant,
    stderr: bool,
    _span: tracing::span::EnteredSpan,
}

/// Column width for dot-padded labels on stdout.
const LABEL_WIDTH: usize = 26;

impl TimedStep {
    /// Dot-pad `label` to stdout, open a tracing span, flush, start clock.
    pub fn begin(label: impl fmt::Display) -> Self {
        let label = label.to_string();
        let span = tracing::info_span!("step", label = label.as_str());
        let dots = LABEL_WIDTH.saturating_sub(label.len() + 1);
        if dots > 0 {
            print!("{label} {:.<width$}", "", width = dots);
        } else {
            print!("{label}");
        }
        std::io::Write::flush(&mut std::io::stdout()).ok();
        Self {
            start: Instant::now(),
            stderr: false,
            _span: span.entered(),
        }
    }

    /// Print `label` to stderr (no dot-padding), open a tracing span, start clock.
    pub fn begin_stderr(label: impl fmt::Display) -> Self {
        let label = label.to_string();
        let span = tracing::info_span!("step", label = label.as_str());
        eprint!("{label}");
        Self {
            start: Instant::now(),
            stderr: true,
            _span: span.entered(),
        }
    }

    /// Print ` {elapsed:?}\n`, return Duration.
    pub fn finish(self) -> Duration {
        let d = self.start.elapsed();
        self.emit(format_args!(" {:?}", d));
        d
    }

    /// Print ` {elapsed:?} ({detail})\n`, return Duration.
    pub fn finish_with(self, detail: impl fmt::Display) -> Duration {
        let d = self.start.elapsed();
        self.emit(format_args!(" {:?} ({})", d, detail));
        d
    }

    fn emit(&self, msg: fmt::Arguments<'_>) {
        if self.stderr {
            eprintln!("{msg}");
        } else {
            println!("{msg}");
        }
    }
}
