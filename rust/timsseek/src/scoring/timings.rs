//! Timing instrumentation for the scoring pipeline.
//!
//! This module provides timing measurement for each stage of the scoring process.
//! The timings are aggregated across parallel scoring operations to provide overall
//! performance metrics.

use serde::Serialize;
use std::time::Duration;

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

/// Full pipeline report: per-phase timings and result-quality metrics.
/// All timing fields are in milliseconds.
#[derive(Debug, Default, Serialize)]
pub struct PipelineReport {
    // Per-file: index loading (ms)
    pub load_index_ms: u64,

    // Phase timings (all in ms)
    pub phase1_prescore_ms: u64,
    pub phase2_calibration_ms: u64,
    pub phase3_extraction_ms: u64,
    pub phase3_scoring_ms: u64,
    pub phase3_spectral_query_ms: u64,
    pub phase3_assembly_ms: u64,
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
