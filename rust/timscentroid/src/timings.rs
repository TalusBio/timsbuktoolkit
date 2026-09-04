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
    MultiProgress,
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
    /// Held back rather than printed at `begin`, for a step that cannot own the
    /// terminal while it runs. `None` once the label is already out.
    deferred_label: Option<String>,
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
    styled_progress_bar(len, label)
}

fn styled_progress_bar(len: u64, label: &str) -> ProgressBar {
    let style = ProgressStyle::with_template(&format!(
        "{{spinner:.green}} {} [{{elapsed_precise}}] [{{wide_bar:.cyan/blue}}] {{pos}}/{{len}} ({{eta}})",
        label
    ))
    .expect("the built-in progress template must be valid");
    ProgressBar::new(len)
        .with_style(style)
        .with_finish(ProgressFinish::AndLeave)
}

/// A TTY-aware set of independently updating progress bars. On captured or
/// redirected stderr, [`Self::add`] returns hidden no-op bars so callers can keep
/// one progress-reporting code path without emitting redraw noise.
pub struct ProgressGroup {
    inner: Option<MultiProgress>,
}

impl ProgressGroup {
    pub fn new() -> Self {
        Self {
            inner: std::io::stderr().is_terminal().then(MultiProgress::new),
        }
    }

    pub fn is_hidden(&self) -> bool {
        self.inner.is_none()
    }

    pub fn add(&self, len: u64, label: &str) -> ProgressBar {
        match &self.inner {
            Some(group) => group.add(styled_progress_bar(len, label)),
            None => ProgressBar::hidden(),
        }
    }
}

impl Default for ProgressGroup {
    fn default() -> Self {
        Self::new()
    }
}

impl TimedStep {
    /// Dot-pad `label` to stdout, open a tracing span, flush, start clock.
    pub fn begin(label: impl fmt::Display) -> Self {
        let label = label.to_string();
        let span = tracing::info_span!("step", label = label.as_str());
        print!("{}", padded(&label));
        std::io::Write::flush(&mut std::io::stdout()).ok();
        Self {
            start: Instant::now(),
            stderr: false,
            deferred_label: None,
            _span: span.entered(),
        }
    }

    /// As [`Self::begin`], but the label waits until `finish`.
    ///
    /// For a step that draws its own progress on stderr. `begin` leaves an open
    /// line on stdout, and a bar drawn into it starts mid-line and then clears a
    /// span of the terminal it does not own -- the label and the bar cannot
    /// share the line, so the label goes last and the whole step is one line
    /// either way.
    pub fn begin_deferred(label: impl fmt::Display) -> Self {
        let label = label.to_string();
        let span = tracing::info_span!("step", label = label.as_str());
        Self {
            start: Instant::now(),
            stderr: false,
            deferred_label: Some(padded(&label)),
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
            deferred_label: None,
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
        let label = self.deferred_label.as_deref().unwrap_or("");
        if self.stderr {
            eprintln!("{label}{msg}");
        } else {
            println!("{label}{msg}");
        }
    }
}

/// `label` followed by dots out to [`LABEL_WIDTH`], so every step's elapsed time
/// lands in the same column.
fn padded(label: &str) -> String {
    let dots = LABEL_WIDTH.saturating_sub(label.len() + 1);
    if dots > 0 {
        format!("{label} {:.<width$}", "", width = dots)
    } else {
        label.to_string()
    }
}
