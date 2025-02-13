use crate::fragment_mass::IonAnnot;
use crate::scoring::calculate_scores::LongitudinalMainScoreElements;
use crate::scoring::search_results::IonSearchResults;
use serde::Serialize;
use timsquery::models::aggregators::raw_peak_agg::multi_chromatogram_agg::NaturalFinalizedMultiCMGArrays;

#[derive(Debug, Clone, Serialize)]
pub struct FullQueryResult {
    pub main_score_elements: LongitudinalMainScoreElements,
    pub longitudinal_main_score: Vec<f32>,
    pub extractions: NaturalFinalizedMultiCMGArrays<IonAnnot>,
    pub search_results: IonSearchResults,
}
