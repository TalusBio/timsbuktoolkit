use egui::Color32;
use std::collections::HashMap;
use timscentroid::IndexedTimstofPeaks;
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
use crate::file_loader::ElutionGroupData;
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

/// Computed/cached state - derived from data and UI state
#[derive(Debug, Default)]
pub struct ComputedState {
    /// Computed chromatogram for the selected elution group (plot data)
    pub chromatogram_lines: Option<ChromatogramLines>,
    pub score_lines: Option<ScoreLines>,
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
    pub expected_intensities: Option<ExpectedIntensities<String>>,
    pub apex_score: Option<ApexScore>,

    chromatogram_collector_buffer: Option<ChromatogramCollector<String, f32>>,
    apex_finder_buffer: Option<ApexFinder>,
    cache_key: Option<(u64, Tolerance, SmoothingMethod)>,
    last_requested_rt: Option<f64>,
    reference_lines: HashMap<String, (f64, Color32)>,
}

impl ComputedState {
    pub fn update(
        &mut self,
        elution_groups: &ElutionGroupData,
        selected_idx: usize,
        index: &IndexedPeaksHandle,
        tolerance: &Tolerance,
        smoothing: &SmoothingMethod,
    ) {
        if !self.is_cache_valid(selected_idx, tolerance, smoothing) {
            self.reset_with_data(elution_groups, selected_idx, index, tolerance, smoothing);
        }

        if let Some(requested_rt) = self.clicked_rt {
            if let Some(last_rt) = self.last_requested_rt
                && (last_rt - requested_rt).abs() < f64::EPSILON
            {
                // Requested RT is the same as last generated, skip regeneration
                return;
            }
            self.generate_spectrum(requested_rt);
            self.reference_lines
                .insert("Clicked RT".into(), (requested_rt, Color32::GREEN));
        }
    }

    pub fn reference_lines(&self) -> &HashMap<String, (f64, Color32)> {
        &self.reference_lines
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
        self.chromatogram_lines = None;
        self.chromatogram_x_bounds = None;
        self.chromatogram_output = None;
        self.ms2_spectrum = None;
        self.auto_zoom_frame_counter = 0;
        self.clicked_rt = None;
        self.cache_key = None;
        self.last_requested_rt = None;
        self.apex_score = None;
        self.reference_lines.clear();
    }

    fn reset_with_data(
        &mut self,
        elution_groups: &ElutionGroupData,
        selected_idx: usize,
        index: &IndexedPeaksHandle,
        tolerance: &Tolerance,
        smoothing: &SmoothingMethod,
    ) {
        // TODO: add cache invalidation to avoid regenerating if nothing changed
        self.clear();
        let (elution_group, expected_intensities) = elution_groups
            .get_elem(selected_idx)
            .expect("Valid elution group index");
        self.cache_key = Some((selected_idx as u64, tolerance.clone(), *smoothing));
        self.expected_intensities = Some(expected_intensities.clone());
        self.reference_lines
            .entry("Library RT".into())
            .or_insert((elution_group.rt_seconds() as f64, Color32::BLUE));

        let collector = match self.chromatogram_collector_buffer.as_mut() {
            // I can probably clean this statement...
            Some(collector) => {
                let max_range = index.ms1_cycle_mapping().range_milis();
                collector
                    .try_reset_with(
                        elution_group.clone(),
                        TupleRange::try_new(max_range.0, max_range.1)
                            .expect("Reference RTs should be sorted and valid"),
                        index.ms1_cycle_mapping(),
                    )
                    .unwrap();
                collector
            }
            None => {
                self.chromatogram_collector_buffer = Some(
                    Self::build_collector(index, elution_group.clone())
                        .expect("Failed to build chromatogram collector"),
                );
                self.chromatogram_collector_buffer.as_mut().unwrap()
            }
        };
        let num_cycles = collector.num_cycles();

        match Self::generate_chromatogram(collector, &elution_group, index, tolerance, smoothing) {
            Ok(output) => {
                let chrom_lines = ChromatogramLines::from_chromatogram(&output);
                self.chromatogram_x_bounds = Some(chrom_lines.rt_seconds_range);
                self.chromatogram_lines = Some(chrom_lines);
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
                return;
            }
        };

        // Prepare apex finder
        let apex_finder = if self.apex_finder_buffer.is_none() {
            self.apex_finder_buffer = Some(ApexFinder::new(num_cycles));
            self.apex_finder_buffer.as_mut().unwrap()
        } else {
            // We dont have to reset, calling the find_apex method will overwrite previous data
            self.apex_finder_buffer.as_mut().unwrap()
        };

        let scoring_ctx = ScoringContext {
            expected_intensities: expected_intensities.clone(),
            query_values: collector.clone(),
        };

        let apex_score = match Self::find_apex(apex_finder, &scoring_ctx, index) {
            Ok(score) => score,
            Err(e) => {
                tracing::error!(
                    "Failed to compute apex score for elution group {}: {:?}",
                    elution_group.id(),
                    e
                );
                return;
            }
        };
        self.clicked_rt = Some(apex_score.retention_time_ms as f64 / 1000.0);
        self.reference_lines
            .entry("Apex RT".into())
            .or_insert((apex_score.retention_time_ms as f64 / 1000.0, Color32::RED));
        self.score_lines = Some(ScoreLines::from_scores(
            apex_score,
            &apex_finder.traces,
            index.ms1_cycle_mapping(),
            collector.cycle_offset(),
        ));
        self.apex_score = Some(apex_score);
    }

    fn build_collector(
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
    fn generate_chromatogram(
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
