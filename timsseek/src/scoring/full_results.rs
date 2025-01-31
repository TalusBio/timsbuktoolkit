
use crate::fragment_mass::fragment_mass_builder::SafePosition;
use crate::scoring::calculate_scores::LongitudinalMainScoreElements;
use crate::scoring::search_results::IonSearchResults;
use timsquery::models::aggregators::raw_peak_agg::multi_chromatogram_agg::NaturalFinalizedMultiCMGArrays;
use serde::Serialize;


#[derive(Debug, Clone, Serialize)]
pub struct FullQueryResult {
    pub main_score_elements: LongitudinalMainScoreElements,
    pub longitudinal_main_score: Vec<f32>,
    pub extractions: NaturalFinalizedMultiCMGArrays<SafePosition>,
    pub search_results: IonSearchResults,
}
