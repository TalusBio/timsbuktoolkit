//! First-class tracking of peptide-skipping reasons across scoring phases.
//!
//! Every scoring phase that can drop a peptide returns `Result<_, SkipReason>`;
//! accumulators fold these into a `SkipCounts` that is merged across Rayon
//! workers and surfaced in `PipelineReport`.

use crate::errors::DataProcessingError;
use serde::Serialize;
use timsquery::errors::DataProcessingError as TQDataProcessingError;

/// Why a peptide failed to produce a scored result in a given phase.
///
/// Variants are ordered roughly by when they are reachable in the pipeline
/// so the JSON field order in `SkipCounts` mirrors the control flow.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SkipReason {
    /// Precursor isolation window lies fully outside the instrument's
    /// fragmented mass range (pre-filter, no extraction).
    PrecursorOutOfFragmentedRange,
    /// Query RT (± tolerance) does not intersect the run's cycle range.
    RetentionTimeOutOfBounds,
    /// Speclib entry carries zero predicted fragments — nothing to score.
    NoExpectedFragments,
    /// Extraction ran but zero quad isolation windows overlapped the query's
    /// precursor range. Signals a library / instrument scan-schedule mismatch
    /// (all fragment m/z queries would miss every window) rather than peptide
    /// absence.
    FragmentsOutsideScanRange,
    /// Extraction ran, quad windows matched, but zero peaks landed in any
    /// chromatogram cell. Common for true-negative peptides. Lets us skip
    /// `filter_zero_intensity_ions` and `find_apex` work.
    NoObservedSignal,
    /// Extraction ran, but `find_apex` / `find_apex_location` saw no observed
    /// intensity at all (all traces empty).
    ApexEmptyData,
    /// Too few cycles in the RT window to compute traces reliably.
    ApexInsufficientData,
    /// NaN / Inf surfaced in scoring math (indicates upstream data pathology).
    ApexNonFiniteScore,
    /// `finalize_results` failed to assemble the ScoredCandidate.
    FinalizeError,
    /// Invariant violation (KeyNotFound / IndexOutOfBounds / VectorLength).
    /// Non-zero counts here indicate a bug, not data quality.
    InvariantViolation,
}

/// Map the timsquery-level error into a `SkipReason` sub-bucket.
fn classify_tq_error(err: &TQDataProcessingError) -> SkipReason {
    match err {
        TQDataProcessingError::ExpectedNonEmptyData => SkipReason::ApexEmptyData,
        TQDataProcessingError::InsufficientData { .. } => SkipReason::ApexInsufficientData,
        TQDataProcessingError::UnexpectedInfiniteError(_)
        | TQDataProcessingError::UnexpectedInfiniteErrors(_) => SkipReason::ApexNonFiniteScore,
        TQDataProcessingError::KeyNotFound
        | TQDataProcessingError::IndexOutOfBoundsError(_)
        | TQDataProcessingError::ExpectedVectorLength { .. }
        | TQDataProcessingError::ExpectedVectorSameLength => SkipReason::InvariantViolation,
    }
}

/// Map an apex-stage `DataProcessingError` into a `SkipReason`.
pub fn apex_error_to_skip(err: &DataProcessingError) -> SkipReason {
    match err {
        DataProcessingError::ExpectedNonEmptyData { .. } => SkipReason::ApexEmptyData,
        DataProcessingError::ExpectedFiniteNonNanData { .. } => SkipReason::ApexNonFiniteScore,
        DataProcessingError::TimsQueryDataProcessingError { error, .. } => classify_tq_error(error),
        DataProcessingError::ExpectedSlicesSameLength { .. }
        | DataProcessingError::ExpectedSetField { .. } => SkipReason::InvariantViolation,
    }
}

/// Same mapping, but scoped to the finalize phase so unmapped failures land
/// in the `FinalizeError` bucket rather than `ApexEmptyData`.
pub fn finalize_error_to_skip(err: &DataProcessingError) -> SkipReason {
    match err {
        DataProcessingError::ExpectedSlicesSameLength { .. }
        | DataProcessingError::ExpectedSetField { .. } => SkipReason::InvariantViolation,
        DataProcessingError::TimsQueryDataProcessingError { error, .. } => {
            match classify_tq_error(error) {
                SkipReason::InvariantViolation => SkipReason::InvariantViolation,
                _ => SkipReason::FinalizeError,
            }
        }
        _ => SkipReason::FinalizeError,
    }
}

/// Per-reason skip counts for a single phase. Mergeable across rayon workers
/// via `AddAssign`; serialized into `PipelineReport`.
#[derive(Debug, Default, Clone, Copy, Serialize)]
pub struct SkipCounts {
    pub precursor_out_of_fragmented_range: u32,
    pub retention_time_out_of_bounds: u32,
    pub no_expected_fragments: u32,
    pub fragments_outside_scan_range: u32,
    pub no_observed_signal: u32,
    pub apex_empty_data: u32,
    pub apex_insufficient_data: u32,
    pub apex_non_finite_score: u32,
    pub finalize_error: u32,
    pub invariant_violation: u32,
}

impl SkipCounts {
    pub fn bump(&mut self, reason: SkipReason) {
        let slot = match reason {
            SkipReason::PrecursorOutOfFragmentedRange => {
                &mut self.precursor_out_of_fragmented_range
            }
            SkipReason::RetentionTimeOutOfBounds => &mut self.retention_time_out_of_bounds,
            SkipReason::NoExpectedFragments => &mut self.no_expected_fragments,
            SkipReason::FragmentsOutsideScanRange => &mut self.fragments_outside_scan_range,
            SkipReason::NoObservedSignal => &mut self.no_observed_signal,
            SkipReason::ApexEmptyData => &mut self.apex_empty_data,
            SkipReason::ApexInsufficientData => &mut self.apex_insufficient_data,
            SkipReason::ApexNonFiniteScore => &mut self.apex_non_finite_score,
            SkipReason::FinalizeError => &mut self.finalize_error,
            SkipReason::InvariantViolation => &mut self.invariant_violation,
        };
        *slot = slot.saturating_add(1);
    }

    pub fn total(&self) -> u64 {
        self.named_fields().iter().map(|(n, _)| *n as u64).sum()
    }

    /// Single source of truth for field enumeration: exhaustive destructure
    /// so adding a `SkipReason` variant (which requires a new field) is a
    /// compile error here, propagating to `total` / `AddAssign` / `Display`.
    fn named_fields(&self) -> [(u32, &'static str); 10] {
        let Self {
            precursor_out_of_fragmented_range,
            retention_time_out_of_bounds,
            no_expected_fragments,
            fragments_outside_scan_range,
            no_observed_signal,
            apex_empty_data,
            apex_insufficient_data,
            apex_non_finite_score,
            finalize_error,
            invariant_violation,
        } = *self;
        [
            (precursor_out_of_fragmented_range, "precursor_oofr"),
            (retention_time_out_of_bounds, "rt_oob"),
            (no_expected_fragments, "no_expected_frags"),
            (fragments_outside_scan_range, "frags_outside_scan"),
            (no_observed_signal, "no_observed_signal"),
            (apex_empty_data, "apex_empty"),
            (apex_insufficient_data, "apex_insufficient"),
            (apex_non_finite_score, "apex_nonfinite"),
            (finalize_error, "finalize_err"),
            (invariant_violation, "invariant"),
        ]
    }
}

impl std::ops::AddAssign for SkipCounts {
    fn add_assign(&mut self, rhs: Self) {
        // Exhaustive destructure on `rhs` so adding a new field forces an
        // update here (the struct literal below would otherwise silently
        // drop the new field into `..`).
        let Self {
            precursor_out_of_fragmented_range,
            retention_time_out_of_bounds,
            no_expected_fragments,
            fragments_outside_scan_range,
            no_observed_signal,
            apex_empty_data,
            apex_insufficient_data,
            apex_non_finite_score,
            finalize_error,
            invariant_violation,
        } = rhs;
        self.precursor_out_of_fragmented_range += precursor_out_of_fragmented_range;
        self.retention_time_out_of_bounds += retention_time_out_of_bounds;
        self.no_expected_fragments += no_expected_fragments;
        self.fragments_outside_scan_range += fragments_outside_scan_range;
        self.no_observed_signal += no_observed_signal;
        self.apex_empty_data += apex_empty_data;
        self.apex_insufficient_data += apex_insufficient_data;
        self.apex_non_finite_score += apex_non_finite_score;
        self.finalize_error += finalize_error;
        self.invariant_violation += invariant_violation;
    }
}

impl std::fmt::Display for SkipCounts {
    /// One-line summary of non-zero buckets. Zero buckets omitted to keep
    /// stderr/log output readable on healthy runs.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut first = true;
        for (n, name) in self.named_fields() {
            if n == 0 {
                continue;
            }
            if !first {
                f.write_str(", ")?;
            }
            first = false;
            write!(f, "{name}={n}")?;
        }
        if first {
            f.write_str("none")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bump_and_total() {
        let mut c = SkipCounts::default();
        c.bump(SkipReason::RetentionTimeOutOfBounds);
        c.bump(SkipReason::RetentionTimeOutOfBounds);
        c.bump(SkipReason::NoExpectedFragments);
        c.bump(SkipReason::ApexNonFiniteScore);
        assert_eq!(c.retention_time_out_of_bounds, 2);
        assert_eq!(c.no_expected_fragments, 1);
        assert_eq!(c.apex_non_finite_score, 1);
        assert_eq!(c.total(), 4);
    }

    #[test]
    fn add_assign_merges() {
        let mut a = SkipCounts::default();
        a.bump(SkipReason::RetentionTimeOutOfBounds);
        let mut b = SkipCounts::default();
        b.bump(SkipReason::RetentionTimeOutOfBounds);
        b.bump(SkipReason::FinalizeError);
        a += b;
        assert_eq!(a.retention_time_out_of_bounds, 2);
        assert_eq!(a.finalize_error, 1);
        assert_eq!(a.total(), 3);
    }

    #[test]
    fn display_skips_zeros_and_joins_nonzero() {
        let mut c = SkipCounts::default();
        assert_eq!(format!("{}", c), "none");
        c.bump(SkipReason::RetentionTimeOutOfBounds);
        c.bump(SkipReason::ApexEmptyData);
        c.bump(SkipReason::ApexEmptyData);
        assert_eq!(format!("{}", c), "rt_oob=1, apex_empty=2");
    }

    #[test]
    fn classifies_timsquery_variants() {
        let tq = TQDataProcessingError::InsufficientData {
            real: 3,
            expected: 10,
        };
        assert_eq!(classify_tq_error(&tq), SkipReason::ApexInsufficientData);
        let tq = TQDataProcessingError::KeyNotFound;
        assert_eq!(classify_tq_error(&tq), SkipReason::InvariantViolation);
    }
}
