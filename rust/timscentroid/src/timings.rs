//! Progressive CLI timing output.
//!
//! [`TimedStep`] prints a dot-padded label immediately, opens a tracing span,
//! then appends elapsed time when the work finishes.

use std::fmt;
use std::io::IsTerminal;
use std::time::{
    Duration,
    Instant,
};

use indicatif::{
    ProgressBar,
    ProgressFinish,
    ProgressStyle,
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

/// Create a progress bar on stderr when it is interactive, or a hidden no-op
/// bar when output is captured or redirected. Keeping the TTY decision here
/// gives library work and CLI phases the same non-spamming behavior.
pub fn make_progress_bar(len: u64, label: &str) -> ProgressBar {
    if !std::io::stderr().is_terminal() {
        return ProgressBar::hidden();
    }
    let style = ProgressStyle::with_template(&format!(
        "{{spinner:.green}} {} [{{elapsed_precise}}] [{{wide_bar:.cyan/blue}}] {{pos}}/{{len}} ({{eta}})",
        label
    ))
    .expect("the built-in progress template must be valid");
    ProgressBar::new(len)
        .with_style(style)
        .with_finish(ProgressFinish::AndLeave)
}

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
