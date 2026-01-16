use egui::Color32;
use std::collections::HashMap;
use timsquery::models::elution_group::TimsElutionGroup;
use timsquery::models::tolerance::Tolerance;
use timsquery::{
    QueriableData,
    TupleRange,
};
use timsseek::ExpectedIntensities;
use timsseek::errors::DataProcessingError;
use tracing::instrument;

use crate::chromatogram_processor;
use crate::chromatogram_processor::{
    ChromatogramCollector,
    ChromatogramOutput,
    SmoothingMethod,
    apply_smoothing_chromatogram,
};
use crate::error::ViewerError;
use crate::plot_renderer::{
    ChromatogramLines,
    MS2Spectrum,
    ScoreLines,
};
use timscentroid::rt_mapping::{
    MS1CycleIndex,
    RTIndex,
};
use timsquery::serde::IndexedPeaksHandle;
use timsseek::scoring::apex_finding::{
    ApexFinder,
    ApexScore,
    ScoringContext,
};

/// Result bundle from background chromatogram computation
pub(crate) struct ChromatogramComputationResult {
    pub selected_idx: u64,
    pub output: ChromatogramOutput,
    pub collector: ChromatogramCollector<String, f32>,
    pub expected_intensities: ExpectedIntensities<String>,
}

/// Computed/cached state - derived from data and UI state
#[derive(Debug, Default)]
pub struct ComputedState {
    /// Computed chromatogram for the selected elution group (plot data)
    pub chromatogram_lines: Option<ChromatogramLines>,
    pub score_lines: Option<ScoreLines>,
    pub chromatogram_x_bounds: Option<(f64, f64)>,
    pub ms2_spectrum: Option<MS2Spectrum>,
    pub auto_zoom_frame_counter: u8,
    pub clicked_rt: Option<f64>,
    pub expected_intensities: Option<ExpectedIntensities<String>>,
    pub apex_score: Option<ApexScore>,

    // Internal state (private)
    is_computing_chromatogram: bool,
    computing_index: Option<u64>,
    chromatogram_output: Option<ChromatogramOutput>,
    chromatogram_collector_buffer: Option<ChromatogramCollector<String, f32>>,
    apex_finder_buffer: Option<ApexFinder>,
    cache_key: Option<(u64, Tolerance, SmoothingMethod)>,
    last_requested_rt: Option<f64>,
    reference_lines: HashMap<String, (f64, Color32)>,
}

impl ComputedState {
    pub fn reference_lines(&self) -> &HashMap<String, (f64, Color32)> {
        &self.reference_lines
    }

    pub fn insert_reference_line(&mut self, name: String, rt: f64, color: Color32) {
        self.reference_lines.insert(name, (rt, color));
    }

    pub fn is_computing(&self) -> bool {
        self.is_computing_chromatogram
    }

    pub fn computing_index(&self) -> Option<u64> {
        self.computing_index
    }

    pub fn start_computing(&mut self, index: u64) {
        self.is_computing_chromatogram = true;
        self.computing_index = Some(index);
    }

    pub fn cancel_computing(&mut self) {
        self.is_computing_chromatogram = false;
        self.computing_index = None;
    }

    pub fn is_cache_valid(
        &self,
        selected_idx: usize,
        tolerance: &Tolerance,
        smoothing: &SmoothingMethod,
    ) -> bool {
        if let Some((cached_id, cached_tolerance, cached_smoothing)) = &self.cache_key {
            return *cached_id == selected_idx as u64
                && cached_tolerance == tolerance
                && cached_smoothing == smoothing;
        }
        false
    }

    /// Resets computed state when data or UI changes significantly
    pub fn clear(&mut self) {
        self.chromatogram_lines = None;
        self.chromatogram_x_bounds = None;
        self.chromatogram_output = None;
        self.ms2_spectrum = None;
        self.auto_zoom_frame_counter = 0;
        self.clicked_rt = None;
        self.cache_key = None;
        self.computing_index = None;
        self.last_requested_rt = None;
        self.apex_score = None;
        self.score_lines = None;
        self.expected_intensities = None;
        self.reference_lines.clear();
        self.is_computing_chromatogram = false;
    }

    pub(crate) fn build_collector(
        index: &IndexedPeaksHandle,
        elution_group: TimsElutionGroup<String>,
    ) -> Result<ChromatogramCollector<String, f32>, ViewerError> {
        let max_range = index.ms1_cycle_mapping().range_milis();
        // Create collector for this elution group
        let collector = ChromatogramCollector::new(
            elution_group,
            TupleRange::try_new(max_range.0, max_range.1)
                .expect("Reference RTs should be sorted and valid"),
            index.ms1_cycle_mapping(),
        )
        .map_err(|e| ViewerError::General(format!("Failed to create collector: {:?}", e)))?;
        Ok(collector)
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
    pub(crate) fn generate_chromatogram(
        collector: &mut ChromatogramCollector<String, f32>,
        elution_group: &TimsElutionGroup<String>,
        index: &IndexedPeaksHandle,
        tolerance: &Tolerance,
        smoothing: &SmoothingMethod,
    ) -> Result<ChromatogramOutput, ViewerError> {
        // Query the index
        index.add_query(collector, tolerance);

        // Convert to output format
        let mut output = ChromatogramOutput::try_new(collector, index.ms1_cycle_mapping())
            .map_err(|e| {
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

    #[instrument(skip_all, fields(eg_id = %context.query_values.eg.id()))]
    fn find_apex(
        apex_finder: &mut ApexFinder,
        context: &ScoringContext<String>,
        index: &IndexedPeaksHandle,
    ) -> Result<ApexScore, ViewerError> {
        apex_finder.find_apex(context, &|idx| {
            index
                .ms1_cycle_mapping()
                .rt_milis_for_index(&MS1CycleIndex::new(idx as u32))
                .unwrap()
        }).map_err(|x| {
            match x {
                DataProcessingError::ExpectedNonEmptyData { context: err_context } => {
                    tracing::warn!(
                        "{:#?}", context.query_values.eg, 
                    );
                    tracing::warn!(
                        "Apex finding failed for elution group {}: No valid data found in context {:?}",
                        context.query_values.eg.id(),
                        err_context
                    );
                }
                _ => {
                    tracing::error!(
                        "Apex finding failed for elution group {}: {:?}",
                        context.query_values.eg.id(),
                        x
                    );
                }
            };
            ViewerError::General("Apex finding error".into())

        })
    }

    /// Generate MS2 spectrum at the given retention time.
    /// Returns true if a new spectrum was generated, false if skipped (already generated for this RT).
    #[instrument(level = "trace", skip(self), fields(eg_id = %self.cache_key.as_ref().map(|(id, _, _)| *id).unwrap_or(0)))]
    pub fn generate_spectrum_at_rt(&mut self, rt_seconds: f64) -> bool {
        // Skip if we already generated spectrum for this RT
        // Use a tolerance of 1 microsecond (1e-6 seconds) instead of f64::EPSILON
        // which is far too strict for RT values and can fail on different FPU implementations
        const RT_TOLERANCE_SECONDS: f64 = 1e-6;
        if let Some(last_rt) = self.last_requested_rt
            && (last_rt - rt_seconds).abs() < RT_TOLERANCE_SECONDS
        {
            return false;
        }

        // Cache miss - generating new MS2 spectrum at this RT
        tracing::debug!(
            "MS2 spectrum cache MISS - generating spectrum at RT {:.2}s (previous: {:?})",
            rt_seconds,
            self.last_requested_rt
        );
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
                    true
                }
                Err(e) => {
                    tracing::error!("Failed to extract MS2 spectrum: {:?}", e);
                    self.ms2_spectrum = None;
                    false
                }
            }
        } else {
            tracing::warn!("No chromatogram data available for MS2 extraction");
            false
        }
    }

    /// Complete chromatogram computation with results from background thread
    pub(crate) fn complete_chromatogram_computation(
        &mut self,
        result: ChromatogramComputationResult,
        index: &IndexedPeaksHandle,
        tolerance: &Tolerance,
        smoothing: SmoothingMethod,
    ) {
        // Store chromatogram output and lines
        let chrom_lines = ChromatogramLines::from_chromatogram(&result.output);
        self.chromatogram_x_bounds = Some(chrom_lines.rt_seconds_range);
        self.chromatogram_lines = Some(chrom_lines);
        self.chromatogram_output = Some(result.output.clone());
        self.auto_zoom_frame_counter = 5;

        // Update cache key
        self.cache_key = Some((result.selected_idx, tolerance.clone(), smoothing));

        // Store for scoring
        self.expected_intensities = Some(result.expected_intensities.clone());
        self.chromatogram_collector_buffer = Some(result.collector.clone());

        // Compute scores on main thread
        self.compute_scores_from_buffers(
            index,
            &result.output,
            &result.expected_intensities,
            &result.collector,
        );

        // Clear computing state
        self.is_computing_chromatogram = false;
        self.computing_index = None;
    }

    /// Compute scores from buffers
    fn compute_scores_from_buffers(
        &mut self,
        index: &IndexedPeaksHandle,
        output: &ChromatogramOutput,
        expected_intensities: &ExpectedIntensities<String>,
        collector: &ChromatogramCollector<String, f32>,
    ) {
        let num_cycles = output.retention_time_results_seconds.len();

        // Prepare apex finder
        let apex_finder = if self.apex_finder_buffer.is_none() {
            self.apex_finder_buffer = Some(ApexFinder::new(num_cycles));
            self.apex_finder_buffer.as_mut().unwrap()
        } else {
            self.apex_finder_buffer.as_mut().unwrap()
        };

        // Build scoring context
        let scoring_ctx = ScoringContext {
            expected_intensities: expected_intensities.clone(),
            query_values: collector.clone(),
        };

        // Find apex
        let apex_score = match Self::find_apex(apex_finder, &scoring_ctx, index) {
            Ok(score) => score,
            Err(e) => {
                tracing::error!("Failed to compute apex score: {:?}", e);
                return;
            }
        };

        // Update RT reference (use insert to overwrite any previous apex RT)
        self.clicked_rt = Some(apex_score.retention_time_ms as f64 / 1000.0);
        self.reference_lines.insert(
            "Apex RT".into(),
            (apex_score.retention_time_ms as f64 / 1000.0, Color32::RED),
        );

        // Generate score lines
        self.score_lines = Some(ScoreLines::from_scores(
            apex_score,
            &apex_finder.traces,
            index.ms1_cycle_mapping(),
            collector.cycle_offset(),
        ));
        self.apex_score = Some(apex_score);
    }
}
