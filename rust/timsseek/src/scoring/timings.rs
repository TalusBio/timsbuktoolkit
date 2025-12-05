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
/// println!("Prescore: {}ms", timings.prescore.as_millis());
/// println!("Localize: {}ms", timings.localize.as_millis());
/// ```
#[derive(Debug, Default)]
pub struct ScoreTimings {
    /// Time spent collecting chromatographic data (Stage 1: Prescore).
    pub prescore: Duration,

    /// Time spent finding peak apex (Stage 2: Localization).
    /// This is typically the bottleneck (~62% of total time).
    pub localize: Duration,

    /// Time spent refining search at detected apex (Stage 3: Secondary Query).
    pub secondary_query: Duration,

    /// Time spent assembling final results (Stage 4: Finalization).
    pub finalization: Duration,
}

impl Serialize for ScoreTimings {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        use serde::ser::SerializeStruct;
        let mut state = serializer.serialize_struct("ScoreTimings", 4)?;
        state.serialize_field("prescore_ms", &self.prescore.as_millis())?;
        state.serialize_field("localize_ms", &self.localize.as_millis())?;
        state.serialize_field("secondary_query_ms", &self.secondary_query.as_millis())?;
        state.serialize_field("finalization_ms", &self.finalization.as_millis())?;
        state.end()
    }
}

impl std::ops::AddAssign for ScoreTimings {
    fn add_assign(&mut self, rhs: Self) {
        self.prescore += rhs.prescore;
        self.localize += rhs.localize;
        self.secondary_query += rhs.secondary_query;
        self.finalization += rhs.finalization;
    }
}
