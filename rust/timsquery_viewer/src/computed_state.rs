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

#[derive(Debug)]
struct ScoringResult {
    apex_score: ApexScore,
    score_lines: ScoreLines,
    apex_rt_seconds: f64,
}

/// All state derived from a single chromatogram computation.
/// Grouped so that replacing one chromatogram with the next drops all derived state.
#[derive(Debug)]
struct ChromatogramResult {
    lines: ChromatogramLines,
    output: ChromatogramOutput,
    scoring: Option<ScoringResult>,
    expected_intensities: ExpectedIntensities<String>,
}

#[derive(Debug, Default)]
struct Ms2State {
    clicked_rt: Option<f64>,
    spectrum: Option<MS2Spectrum>,
    last_requested_rt: Option<f64>,
}

#[derive(Debug, Default)]
struct ScratchBuffers {
    apex_finder: Option<ApexFinder>,
}

/// State machine for chromatogram computation lifecycle.
///
/// Follows the same pattern as `IndexedDataState` and `ElutionGroupState`:
/// `None → Computing → Computed | Failed`
#[derive(Debug, Default)]
enum ChromatogramState {
    #[default]
    None,
    /// Background computation in progress
    Computing { index: u64 },
    /// Computation completed successfully; cache key stored for invalidation checks
    Computed {
        cache_key: (u64, Tolerance, SmoothingMethod),
    },
    /// Computation failed for this combination — prevents infinite retry
    Failed {
        cache_key: (u64, Tolerance, SmoothingMethod),
        error: String,
    },
}

/// Computed/cached state - derived from data and UI state
#[derive(Debug, Default)]
pub struct ComputedState {
    chromatogram_state: ChromatogramState,

    result: Option<ChromatogramResult>,
    pub auto_zoom_frame_counter: u8,
    reference_lines: HashMap<String, (f64, Color32)>,
    ms2: Ms2State,
    /// Reusable allocation — not reset between elution groups.
    scratch: ScratchBuffers,
}

impl ComputedState {
    pub fn ms2_spectrum(&self) -> Option<&MS2Spectrum> {
        self.ms2.spectrum.as_ref()
    }

    pub fn expected_intensities(&self) -> Option<&ExpectedIntensities<String>> {
        self.result.as_ref().map(|r| &r.expected_intensities)
    }

    pub fn clicked_rt(&self) -> Option<f64> {
        self.ms2.clicked_rt
    }

    pub fn set_clicked_rt(&mut self, rt: f64) {
        self.ms2.clicked_rt = Some(rt);
    }

    /// Borrows chromatogram data and the zoom counter simultaneously.
    /// Split out as a method so the borrow checker can see they are disjoint fields.
    pub fn chromatogram_render_context(
        &mut self,
    ) -> Option<(&ChromatogramLines, Option<&ApexScore>, &mut u8)> {
        let result = self.result.as_ref()?;
        let lines = &result.lines;
        let apex = result.scoring.as_ref().map(|s| &s.apex_score);
        let counter = &mut self.auto_zoom_frame_counter;
        Some((lines, apex, counter))
    }

    /// Same idea as [`Self::chromatogram_render_context`] for the score plot.
    pub fn score_render_context(&mut self) -> Option<(&ScoreLines, &mut u8)> {
        let score_lines = &self.result.as_ref()?.scoring.as_ref()?.score_lines;
        let counter = &mut self.auto_zoom_frame_counter;
        Some((score_lines, counter))
    }

    pub fn reference_lines(&self) -> &HashMap<String, (f64, Color32)> {
        &self.reference_lines
    }

    pub fn insert_reference_line(&mut self, name: String, rt: f64, color: Color32) {
        self.reference_lines.insert(name, (rt, color));
    }
}

impl ComputedState {
    pub fn is_computing(&self) -> bool {
        matches!(self.chromatogram_state, ChromatogramState::Computing { .. })
    }

    pub fn computing_index(&self) -> Option<u64> {
        match &self.chromatogram_state {
            ChromatogramState::Computing { index } => Some(*index),
            _ => None,
        }
    }

    pub fn start_computing(&mut self, index: u64) {
        self.chromatogram_state = ChromatogramState::Computing { index };
    }

    /// Record a failed computation so the same parameters won't be retried.
    pub fn fail_computing(
        &mut self,
        index: u64,
        tolerance: &Tolerance,
        smoothing: &SmoothingMethod,
        error: String,
    ) {
        self.chromatogram_state = ChromatogramState::Failed {
            cache_key: (index, tolerance.clone(), *smoothing),
            error,
        };
    }

    /// Returns the error message if in `Failed` state.
    pub fn failure_error(&self) -> Option<&str> {
        match &self.chromatogram_state {
            ChromatogramState::Failed { error, .. } => Some(error),
            _ => None,
        }
    }

    /// Returns the cached elution group ID (for tracing instrumentation).
    fn cached_eg_id(&self) -> u64 {
        match &self.chromatogram_state {
            ChromatogramState::Computed { cache_key } => cache_key.0,
            ChromatogramState::Failed { cache_key, .. } => cache_key.0,
            _ => 0,
        }
    }

    pub fn is_cache_valid(
        &self,
        selected_idx: usize,
        tolerance: &Tolerance,
        smoothing: &SmoothingMethod,
    ) -> bool {
        let cache_key = match &self.chromatogram_state {
            ChromatogramState::Computed { cache_key } => cache_key,
            ChromatogramState::Failed { cache_key, .. } => cache_key,
            _ => return false,
        };
        cache_key.0 == selected_idx as u64 && &cache_key.1 == tolerance && &cache_key.2 == smoothing
    }

    /// Resets computed state when data or UI changes significantly.
    /// Scratch buffers are preserved for reuse.
    pub fn clear(&mut self) {
        self.result = None;
        self.ms2 = Ms2State::default();
        self.auto_zoom_frame_counter = 0;
        self.reference_lines.clear();
        self.chromatogram_state = ChromatogramState::None;
    }

    pub(crate) fn build_collector(
        index: &IndexedPeaksHandle,
        elution_group: TimsElutionGroup<String>,
    ) -> Result<ChromatogramCollector<String, f32>, ViewerError> {
        let max_range = index.ms1_cycle_mapping().range_milis();
        let collector = ChromatogramCollector::new(
            elution_group,
            TupleRange::try_new(max_range.0, max_range.1)
                .expect("Reference RTs should be sorted and valid"),
            index.ms1_cycle_mapping(),
        )
        .map_err(|e| ViewerError::General(format!("Failed to create collector: {:?}", e)))?;
        Ok(collector)
    }

    #[instrument(skip_all, fields(eg_id = %elution_group.id()))]
    pub(crate) fn generate_chromatogram(
        collector: &mut ChromatogramCollector<String, f32>,
        elution_group: &TimsElutionGroup<String>,
        index: &IndexedPeaksHandle,
        tolerance: &Tolerance,
        smoothing: &SmoothingMethod,
    ) -> Result<ChromatogramOutput, ViewerError> {
        index.add_query(collector, tolerance);

        let mut output = ChromatogramOutput::try_new(collector, index.ms1_cycle_mapping())
            .map_err(|e| {
                let msg = match e {
                    timsquery::errors::DataProcessingError::ExpectedNonEmptyData => {
                        "No data found with the current tolerances. \
                         Try widening the mass or mobility tolerance, \
                         or removing the retention time restriction."
                            .to_string()
                    }
                    other => format!("Failed to generate chromatogram output: {:?}", other),
                };
                ViewerError::General(msg)
            })?;

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
            let user_msg = match &x {
                DataProcessingError::ExpectedNonEmptyData { context: err_context } => {
                    tracing::warn!(
                        "{:#?}", context.query_values.eg,
                    );
                    tracing::warn!(
                        "Apex finding failed for elution group {}: No valid data found in context {:?}",
                        context.query_values.eg.id(),
                        err_context
                    );
                    "No data found with the current tolerances. \
                     Try widening the mass or mobility tolerance, \
                     or removing the retention time restriction."
                        .to_string()
                }
                _ => {
                    tracing::error!(
                        "Apex finding failed for elution group {}: {:?}",
                        context.query_values.eg.id(),
                        x
                    );
                    format!("Apex finding failed: {:?}", x)
                }
            };
            ViewerError::General(user_msg)
        })
    }

    /// Generate MS2 spectrum at the given retention time.
    /// Returns true if a new spectrum was generated, false if skipped (already generated for this RT).
    #[instrument(level = "trace", skip(self), fields(eg_id = %self.cached_eg_id()))]
    pub fn generate_spectrum_at_rt(&mut self, rt_seconds: f64) -> bool {
        const RT_TOLERANCE_SECONDS: f64 = 1e-6;
        if let Some(last_rt) = self.ms2.last_requested_rt
            && (last_rt - rt_seconds).abs() < RT_TOLERANCE_SECONDS
        {
            return false;
        }

        tracing::debug!(
            "MS2 spectrum cache MISS - generating spectrum at RT {:.2}s (previous: {:?})",
            rt_seconds,
            self.ms2.last_requested_rt
        );
        self.ms2.last_requested_rt = Some(rt_seconds);

        let Some(result) = &self.result else {
            tracing::warn!("No chromatogram data available for MS2 extraction");
            return false;
        };

        match chromatogram_processor::extract_ms2_spectrum_from_chromatogram(
            &result.output,
            rt_seconds,
        ) {
            Ok(spectrum) => {
                let num_peaks = spectrum.mz_values.len();
                self.ms2.spectrum = Some(spectrum);
                tracing::info!(
                    "Extracted MS2 spectrum at RT {:.2}s with {} peaks",
                    rt_seconds,
                    num_peaks
                );
                true
            }
            Err(e) => {
                tracing::error!("Failed to extract MS2 spectrum: {:?}", e);
                self.ms2.spectrum = None;
                false
            }
        }
    }

    fn try_score(
        scratch: &mut ScratchBuffers,
        index: &IndexedPeaksHandle,
        output: &ChromatogramOutput,
        expected_intensities: &ExpectedIntensities<String>,
        collector: &ChromatogramCollector<String, f32>,
    ) -> Option<ScoringResult> {
        let num_cycles = output.retention_time_results_seconds.len();
        tracing::debug!(
            "try_score: num_cycles={} num_fragments={}",
            num_cycles,
            collector.fragments.mz_order.len(),
        );

        let apex_finder = scratch
            .apex_finder
            .get_or_insert_with(|| ApexFinder::new(num_cycles));

        let scoring_ctx = ScoringContext {
            expected_intensities: expected_intensities.clone(),
            query_values: collector.clone(),
        };

        let apex_score = match Self::find_apex(apex_finder, &scoring_ctx, index) {
            Ok(score) => {
                tracing::info!("Apex found at RT={:.2}ms", score.retention_time_ms);
                score
            }
            Err(e) => {
                tracing::warn!(
                    "Apex scoring skipped (chromatogram data still available): {:?}",
                    e
                );
                return None;
            }
        };

        let apex_rt_seconds = apex_score.retention_time_ms as f64 / 1000.0;

        let score_lines = ScoreLines::from_scores(
            apex_score,
            &apex_finder.traces,
            index.ms1_cycle_mapping(),
            collector.cycle_offset(),
        );

        Some(ScoringResult {
            apex_score,
            score_lines,
            apex_rt_seconds,
        })
    }

    /// Finalize a background chromatogram computation: build plot lines,
    /// attempt scoring, and replace `self.result` in one shot.
    pub(crate) fn complete_chromatogram_computation(
        &mut self,
        computation: ChromatogramComputationResult,
        index: &IndexedPeaksHandle,
        tolerance: &Tolerance,
        smoothing: SmoothingMethod,
    ) {
        let (n_rt, n_prec, n_frag) = (
            computation.output.retention_time_results_seconds.len(),
            computation.output.precursor_mzs.len(),
            computation.output.fragment_mzs.len(),
        );
        tracing::info!(
            "complete_chromatogram_computation: idx={} rt_points={} precursors={} fragments={} \
             est_output={:.1}KB est_collector={:.1}KB",
            computation.selected_idx,
            n_rt,
            n_prec,
            n_frag,
            (n_rt * (n_prec + n_frag) * std::mem::size_of::<f64>()) as f64 / 1024.0,
            (computation.collector.fragments.mz_order.len() * n_rt * std::mem::size_of::<f32>())
                as f64
                / 1024.0,
        );

        let lines = ChromatogramLines::from_chromatogram(&computation.output);
        let scoring = Self::try_score(
            &mut self.scratch,
            index,
            &computation.output,
            &computation.expected_intensities,
            &computation.collector,
        );

        let apex_rt_seconds = scoring.as_ref().map(|s| s.apex_rt_seconds);
        self.result = Some(ChromatogramResult {
            lines,
            output: computation.output.clone(),
            scoring,
            expected_intensities: computation.expected_intensities.clone(),
        });

        self.ms2 = Ms2State::default();
        self.reference_lines.clear();
        self.auto_zoom_frame_counter = 5;
        self.chromatogram_state = ChromatogramState::Computed {
            cache_key: (computation.selected_idx, tolerance.clone(), smoothing),
        };

        if let Some(apex_rt) = apex_rt_seconds {
            self.ms2.clicked_rt = Some(apex_rt);
            self.reference_lines
                .insert("Apex RT".into(), (apex_rt, Color32::RED));
        }

        tracing::info!(
            "complete_chromatogram_computation finished for idx={}",
            computation.selected_idx
        );
    }
}
