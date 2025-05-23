use crate::IonAnnot;
use crate::scoring::calculate_scores::LongitudinalMainScoreElements;
use crate::scoring::search_results::IonSearchResults;
use serde::Serialize;
use timsquery::models::aggregators::ChromatogramCollector;

#[derive(Debug, Clone, Serialize)]
pub struct FullQueryResult {
    pub main_score_elements: LongitudinalMainScoreElements,
    pub longitudinal_main_score: Vec<f32>,
    pub extractions: ChromatogramCollector<IonAnnot, f32>,
    pub search_results: IonSearchResults,
}
