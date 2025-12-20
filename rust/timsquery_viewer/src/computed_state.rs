use std::collections::HashMap;
use std::sync::Arc;
use timscentroid::IndexedTimstofPeaks;
use timsquery::models::elution_group::TimsElutionGroup;
use timsquery::models::tolerance::Tolerance;
use timsquery::{
    KeyLike,
    QueriableData,
};
use tracing::instrument;

use crate::chromatogram_processor;
use crate::chromatogram_processor::{
    ChromatogramCollector,
    ChromatogramOutput,
    SmoothingMethod,
    apply_smoothing_chromatogram,
};
use crate::error::ViewerError;
use crate::file_loader::ElutionGroupData;
use crate::plot_renderer::{
    ChromatogramLines,
    MS2Spectrum,
};

/// Computed/cached state - derived from data and UI state
#[derive(Debug, Default)]
pub struct ComputedState {
    /// Computed chromatogram for the selected elution group (plot data)
    pub chromatogram: Option<ChromatogramLines>,
    /// X-axis bounds to apply on next plot render (min_rt, max_rt)
    pub chromatogram_x_bounds: Option<(f64, f64)>,
    /// Raw chromatogram output data (for MS2 extraction)
    pub chromatogram_output: Option<ChromatogramOutput>,
    /// Computed MS2 spectrum at selected RT
    pub ms2_spectrum: Option<MS2Spectrum>,
    /// Frame counter for auto-zooming the chromatogram
    pub auto_zoom_frame_counter: u8,
    /// Last clicked RT for MS2 spectrum extraction
    pub clicked_rt: Option<f64>,
    pub expected_intensities: Option<HashMap<String, f32>>,

    cache_key: Option<(u64, Tolerance, SmoothingMethod)>,
    last_requested_rt: Option<f64>,
}

impl ComputedState {
    pub fn update(
        &mut self,
        elution_groups: &ElutionGroupData,
        selected_idx: usize,
        index: &IndexedTimstofPeaks,
        ms1_rts: Arc<[u32]>,
        tolerance: &Tolerance,
        smoothing: &SmoothingMethod,
    ) {
        if !self.is_cache_valid(selected_idx, tolerance, smoothing) {
            self.reset_with_data(
                elution_groups,
                selected_idx,
                index,
                ms1_rts,
                tolerance,
                smoothing,
            );
        }

        if let Some(requested_rt) = self.clicked_rt {
            if let Some(last_rt) = self.last_requested_rt
                && (last_rt - requested_rt).abs() < f64::EPSILON
            {
                // Requested RT is the same as last generated, skip regeneration
                return;
            }
            self.generate_spectrum(requested_rt);
        }
    }

    fn is_cache_valid(
        &self,
        selected_idx: usize,
        tolerance: &Tolerance,
        smoothing: &SmoothingMethod,
    ) -> bool {
        // Check that the index is the same .... not urgent ...
        if let Some((cached_id, cached_tolerance, cached_smoothing)) = &self.cache_key {
            return *cached_id == selected_idx as u64
                && cached_tolerance == tolerance
                && cached_smoothing == smoothing;
        }
        false
    }

    /// Resets computed state when data or UI changes significantly
    fn clear(&mut self) {
        // TODO: Re-write as exhaustive pattern matching so I dont forget any field
        self.chromatogram = None;
        self.chromatogram_x_bounds = None;
        self.chromatogram_output = None;
        self.ms2_spectrum = None;
        self.auto_zoom_frame_counter = 0;
        self.clicked_rt = None;
        self.cache_key = None;
        self.last_requested_rt = None;
    }

    // TODO: replace form this being public to calling update.
    // Which internally calls this but handles cahcing...
    fn reset_with_data(
        &mut self,
        elution_groups: &ElutionGroupData,
        selected_idx: usize,
        index: &IndexedTimstofPeaks,
        ms1_rts: Arc<[u32]>,
        tolerance: &Tolerance,
        smoothing: &SmoothingMethod,
    ) {
        // TODO: add cache invalidation to avoid regenerating if nothing changed
        self.clear();
        let (elution_group, expected_intensities) = elution_groups.get_elem(selected_idx);
        let elution_group = elution_group.unwrap();
        self.cache_key = Some((selected_idx as u64, tolerance.clone(), *smoothing));
        self.expected_intensities = expected_intensities;
        match Self::generate_chromatogram(&elution_group, index, ms1_rts, tolerance, smoothing) {
            Ok(output) => {
                let chrom_lines = ChromatogramLines::from_chromatogram(&output);
                self.chromatogram_x_bounds = Some(chrom_lines.rt_seconds_range);
                self.chromatogram = Some(chrom_lines);
                self.chromatogram_output = Some(output);
                self.clicked_rt = Some(elution_group.rt_seconds() as f64);
                // Reset auto zoom frame counter to 5 to force bounds update for a few frames
                self.auto_zoom_frame_counter = 5;
            }
            Err(e) => {
                tracing::error!(
                    "Failed to generate chromatogram for elution group {}: {:?}",
                    elution_group.id(),
                    e
                );
            }
        }
    }

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
    fn generate_chromatogram<T: KeyLike + std::fmt::Display>(
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

    #[instrument(skip(self), fields(eg_id = %self.cache_key.as_ref().map(|(id, _, _)| *id).unwrap_or(0)))]
    fn generate_spectrum(&mut self, rt_seconds: f64) {
        self.last_requested_rt = Some(rt_seconds);
        if let Some(chrom_output) = &self.chromatogram_output {
            match chromatogram_processor::extract_ms2_spectrum_from_chromatogram(
                chrom_output,
                rt_seconds,
            ) {
                Ok(spectrum) => {
                    let num_peaks = spectrum.mz_values.len();
                    self.ms2_spectrum = Some(spectrum);
                    tracing::info!(
                        "Extracted MS2 spectrum at RT {:.2}s with {} peaks",
                        rt_seconds,
                        num_peaks
                    );
                }
                Err(e) => {
                    tracing::error!("Failed to extract MS2 spectrum: {:?}", e);
                    self.ms2_spectrum = None;
                }
            }
        } else {
            tracing::warn!("No chromatogram data available for MS2 extraction");
        }
    }
}
