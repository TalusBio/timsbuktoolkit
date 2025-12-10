//! Chromatogram generation service

use std::sync::Arc;
use timscentroid::IndexedTimstofPeaks;
use timsquery::models::elution_group::TimsElutionGroup;
use timsquery::models::tolerance::Tolerance;
use timsquery::{
    KeyLike,
    QueriableData,
};
use tracing::instrument;

use crate::chromatogram_processor::{
    ChromatogramCollector,
    ChromatogramOutput,
    SmoothingMethod,
    apply_smoothing_chromatogram,
};
use crate::error::ViewerError;

/// Service for generating chromatograms from elution groups
pub struct ChromatogramService;

impl ChromatogramService {
    /// Generate a chromatogram for a single elution group
    ///
    /// # Arguments
    /// * `elution_group` - The elution group to generate a chromatogram for
    /// * `index` - Indexed timsTOF peaks data
    /// * `ms1_rts` - MS1 retention times in milliseconds
    /// * `tolerance` - Tolerance settings for querying
    /// * `smoothing` - Smoothing method to apply
    ///
    /// # Returns
    /// A `ChromatogramOutput` containing the generated chromatogram data
    #[instrument(skip_all, fields(eg_id = %elution_group.id()))]
    pub fn generate<T: KeyLike + std::fmt::Display>(
        elution_group: &TimsElutionGroup<T>,
        index: &IndexedTimstofPeaks,
        ms1_rts: Arc<[u32]>,
        tolerance: &Tolerance,
        smoothing: &SmoothingMethod,
    ) -> Result<ChromatogramOutput, ViewerError> {
        // Determine RT range for the query
        let rt_range_ms = match tolerance.rt_range_as_milis(elution_group.rt_seconds()) {
            timsquery::OptionallyRestricted::Unrestricted => {
                timsquery::TupleRange::try_new(*ms1_rts.first().unwrap(), *ms1_rts.last().unwrap())
                    .expect("Reference RTs should be sorted and valid")
            }
            timsquery::OptionallyRestricted::Restricted(r) => r,
        };

        // Create collector for this elution group
        let mut collector =
            ChromatogramCollector::new(elution_group.clone(), rt_range_ms, &ms1_rts).map_err(
                |e| ViewerError::General(format!("Failed to create collector: {:?}", e)),
            )?;

        // Query the index
        index.add_query(&mut collector, tolerance);

        // Convert to output format
        let mut output = ChromatogramOutput::try_new(collector, &ms1_rts).map_err(|e| {
            ViewerError::General(format!("Failed to generate chromatogram output: {:?}", e))
        })?;

        // Apply smoothing if configured
        apply_smoothing_chromatogram(&mut output, smoothing);

        tracing::info!(
            "Generated chromatogram for elution group {} with {} precursors and {} fragments",
            output.id,
            output.precursor_mzs.len(),
            output.fragment_mzs.len()
        );

        Ok(output)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Note: These tests require test data fixtures.
    // For now, we're documenting the test structure.
    // In a real implementation, you'd use test fixtures or mocks.

    #[test]
    fn test_chromatogram_service_structure() {
        // This test just verifies the service exists and is well-formed
        let _service = ChromatogramService;
        // Future tests would:
        // 1. Load test data
        // 2. Call ChromatogramService::generate
        // 3. Assert on the output structure
    }

    // Future test ideas:
    // #[test]
    // fn test_generate_chromatogram_with_smoothing() { ... }
    //
    // #[test]
    // fn test_generate_chromatogram_without_smoothing() { ... }
    //
    // #[test]
    // fn test_generate_chromatogram_invalid_rt_range() { ... }
}
