use super::DigestSlice;
use micromzpaf::IonAnnot;
use serde::{
    Deserialize,
    Serialize,
};
use std::collections::HashMap;
use timsquery::ElutionGroup;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExpectedIntensities {
    // TODO: ... should I make this a Vec?
    // since i converted the elution group fragments from a
    // hashmap to a vec it would feel more consistent...
    pub fragment_intensities: HashMap<IonAnnot, f32>,
    pub precursor_intensities: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct QueryItemToScore {
    pub digest: DigestSlice,
    pub charge: u8,
    pub query: ElutionGroup<IonAnnot>,
    pub expected_intensity: ExpectedIntensities,
}

impl QueryItemToScore {
    pub fn sample() -> Self {
        let eg = ElutionGroup {
            id: 42,
            mobility: 0.75,
            rt_seconds: 123.4,
            precursors: vec![(-1, 450.0), (0, 450.5), (1, 451.0), (2, 451.5)].into(),
            fragments: vec![
                (IonAnnot::try_from("y1").unwrap(), 450.0),
                (IonAnnot::try_from("y2").unwrap(), 650.5),
                (IonAnnot::try_from("y3").unwrap(), 751.0),
                (IonAnnot::try_from("y4").unwrap(), 751.5),
            ]
            .into(),
        };
        let ei = ExpectedIntensities {
            fragment_intensities: HashMap::from_iter(
                [
                    (IonAnnot::try_from("y1").unwrap(), 1.0),
                    (IonAnnot::try_from("y2").unwrap(), 1.0),
                    (IonAnnot::try_from("y3").unwrap(), 1.0),
                    (IonAnnot::try_from("y4").unwrap(), 1.0),
                ]
                .iter()
                .cloned(),
            ),
            precursor_intensities: vec![0.0, 1.0, 0.4, 0.1],
        };
        let pepseq = "PEPTIDEPINKPEPTIDE".into();
        let digest = DigestSlice::from_string(pepseq, false, 1);
        let charge = 2;
        let query = eg;
        let expected_intensity = ei;
        QueryItemToScore {
            digest,
            charge,
            query,
            expected_intensity,
        }
    }
}
