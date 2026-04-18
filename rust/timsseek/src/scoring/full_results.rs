use super::apex_finding::ElutionTraces;
use super::results::ScoredCandidate;
use crate::IonAnnot;
use serde::Serialize;
use timsquery::models::aggregators::ChromatogramCollector;

#[derive(Debug, Clone, Serialize)]
pub struct ViewerResult {
    pub traces: ElutionTraces,
    pub longitudinal_apex_profile: Vec<f32>,
    pub chromatograms: ChromatogramCollector<IonAnnot, f32>,
    pub scored: ScoredCandidate,
}
