use crate::IonAnnot;
use crate::scoring::apex_finding::ElutionTraces;
use crate::scoring::results::ScoredCandidate;
use serde::Serialize;
use timsquery::models::aggregators::ChromatogramCollector;

#[derive(Debug, Clone, Serialize)]
pub struct FullQueryResult {
    pub main_score_elements: ElutionTraces,
    pub longitudinal_main_score: Vec<f32>,
    pub extractions: ChromatogramCollector<IonAnnot, f32>,
    pub search_results: ScoredCandidate,
}
