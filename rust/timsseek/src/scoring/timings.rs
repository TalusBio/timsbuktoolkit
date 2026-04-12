//! Timing instrumentation for the scoring pipeline.
//!
//! This module provides timing measurement for each stage of the scoring process.
//! The timings are aggregated across parallel scoring operations to provide overall
//! performance metrics.
//!
//! Two helpers reduce timing boilerplate:
//!
//! - [`timed!`] — times a block and accumulates elapsed into a `Duration` field.
//! - [`TimedStep`] — progressive CLI output: prints a label immediately, then
//!   appends elapsed time when the work finishes.

use serde::Serialize;
use std::fmt;
use std::time::{Duration, Instant};

/// Time a block, accumulate elapsed into `$target`, return the block's value.
///
/// ```ignore
/// let ctx = timed!(timings.extraction, {
///     build_extraction(&query, &index, &tolerance)
/// })?;
/// ```
#[doc(hidden)]
#[macro_export]
macro_rules! timed {
    ($target:expr, $body:expr) => {{
        let __start = std::time::Instant::now();
        let __result = $body;
        $target += __start.elapsed();
        __result
    }};
}

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
        Self { start: Instant::now(), stderr: false, _span: span.entered() }
    }

    /// Print `label` to stderr (no dot-padding), open a tracing span, start clock.
    pub fn begin_stderr(label: impl fmt::Display) -> Self {
        let label = label.to_string();
        let span = tracing::info_span!("step", label = label.as_str());
        eprint!("{label}");
        Self { start: Instant::now(), stderr: true, _span: span.entered() }
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

/// Accumulated timing measurements for the four scoring stages.
///
/// Each field tracks the total time spent in that stage across all queries.
/// When scoring in parallel, timings from all threads are aggregated.
///
/// # Example
///
/// ```ignore
/// let (results, timings) = scorer.score_iter(&queries);
/// println!("Extraction: {}ms", timings.extraction.as_millis());
/// println!("Scoring: {}ms", timings.scoring.as_millis());
/// ```
#[derive(Debug, Default)]
pub struct ScoreTimings {
    /// Time spent collecting chromatographic data (Stage 1: Extraction).
    pub extraction: Duration,

    /// Time spent finding peak apex (Stage 2: Scoring).
    /// This is typically the bottleneck (~62% of total time).
    pub scoring: Duration,

    /// Time spent refining search at detected apex (Stage 3: Spectral Query).
    pub spectral_query: Duration,

    /// Time spent assembling final results (Stage 4: Assembly).
    pub assembly: Duration,
}

impl Serialize for ScoreTimings {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        use serde::ser::SerializeStruct;
        let mut state = serializer.serialize_struct("ScoreTimings", 4)?;
        state.serialize_field("extraction_ms", &self.extraction.as_millis())?;
        state.serialize_field("scoring_ms", &self.scoring.as_millis())?;
        state.serialize_field("spectral_query_ms", &self.spectral_query.as_millis())?;
        state.serialize_field("assembly_ms", &self.assembly.as_millis())?;
        state.end()
    }
}

impl std::ops::AddAssign for ScoreTimings {
    fn add_assign(&mut self, rhs: Self) {
        self.extraction += rhs.extraction;
        self.scoring += rhs.scoring;
        self.spectral_query += rhs.spectral_query;
        self.assembly += rhs.assembly;
    }
}

/// Timing breakdown for Phase 1 prescore.
#[derive(Debug, Default)]
pub struct PrescoreTimings {
    /// Time spent building extractions (chromatogram collection).
    pub extraction: Duration,
    /// Time spent in apex scoring (find_apex_location).
    pub scoring: Duration,
    /// Number of items that passed the fragmented_range filter.
    pub n_passed_filter: usize,
    /// Number of items where prescore returned Some (successful apex).
    pub n_scored: usize,
}

impl Serialize for PrescoreTimings {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        use serde::ser::SerializeStruct;
        let mut state = serializer.serialize_struct("PrescoreTimings", 4)?;
        state.serialize_field("extraction_thread_ms", &self.extraction.as_millis())?;
        state.serialize_field("scoring_thread_ms", &self.scoring.as_millis())?;
        state.serialize_field("n_passed_filter", &self.n_passed_filter)?;
        state.serialize_field("n_scored", &self.n_scored)?;
        state.end()
    }
}

impl std::ops::AddAssign for PrescoreTimings {
    fn add_assign(&mut self, rhs: Self) {
        self.extraction += rhs.extraction;
        self.scoring += rhs.scoring;
        self.n_passed_filter += rhs.n_passed_filter;
        self.n_scored += rhs.n_scored;
    }
}

/// Full pipeline report: per-phase timings and result-quality metrics.
/// All timing fields are in milliseconds.
#[derive(Debug, Default, Serialize)]
pub struct PipelineReport {
    // Per-file: index loading (ms)
    pub load_index_ms: u64,

    // Phase timings (all in ms)
    pub phase1_prescore_ms: u64,
    pub phase1_detail: PrescoreTimings,
    pub phase2_calibration_ms: u64,
    pub phase3_extraction_thread_ms: u64,
    pub phase3_scoring_thread_ms: u64,
    pub phase3_spectral_query_thread_ms: u64,
    pub phase3_assembly_thread_ms: u64,
    pub phase4_competition_ms: u64,
    pub phase5_rescore_ms: u64,
    pub phase6_output_ms: u64,

    // Result quality
    pub total_scored: usize,
    pub total_after_competition: usize,
    pub targets_at_1pct_qval: usize,
    pub targets_at_5pct_qval: usize,
    pub targets_at_10pct_qval: usize,
}

/// Top-level report for an entire CLI invocation.
/// Contains shared loading costs and per-file pipeline reports.
#[derive(Debug, Default, Serialize)]
pub struct RunReport {
    pub load_speclib_ms: u64,
    pub load_calib_lib_ms: u64,
    pub speclib_entries: usize,
    pub calib_lib_entries: usize,
    pub files: Vec<FileReport>,
}

/// Per-file report: file name + pipeline report.
#[derive(Debug, Serialize)]
pub struct FileReport {
    pub file_name: String,
    pub pipeline: PipelineReport,
}
