use super::DigestSlice;
use micromzpaf::IonAnnot;
use serde::{
    Deserialize,
    Serialize,
};
use std::collections::HashMap;
use timsquery::tinyvec::tiny_vec;
use timsquery::{
    KeyLike,
    TimsElutionGroup,
};

// TODO: reimplement my own "keyed_vec" (essentially a vec that enforces
// unique keys on insertion) to avoid the HashMap overhead here.

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExpectedIntensities<T: KeyLike> {
    pub fragment_intensities: HashMap<T, f32>,
    pub precursor_intensities: HashMap<i8, f32>,
}

#[derive(Debug, Clone)]
pub struct QueryItemToScore {
    // Kinda hate this
    pub digest: DigestSlice,
    pub query: TimsElutionGroup<IonAnnot>,
    pub expected_intensity: ExpectedIntensities<IonAnnot>,
}

impl QueryItemToScore {
    pub fn sample() -> Self {
        let eg = TimsElutionGroup::builder()
            .id(42)
            .mobility_ook0(0.75)
            .rt_seconds(123.4)
            .fragment_labels(
                [
                    IonAnnot::try_from("y1").unwrap(),
                    IonAnnot::try_from("y2").unwrap(),
                    IonAnnot::try_from("y3").unwrap(),
                    IonAnnot::try_from("y4").unwrap(),
                ]
                .as_slice()
                .into(),
            )
            .fragment_mzs(vec![450.0, 650.5, 751.0, 751.5])
            .precursor_labels(tiny_vec!(-1, 0, 1, 2))
            .precursor(450.5, 2)
            .try_build()
            .unwrap();

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
            precursor_intensities: HashMap::from_iter(
                [(-1, 0.5), (0, 1.0), (1, 0.8), (2, 0.3)].iter().cloned(),
            ),
        };
        let pepseq = "PEPTIDEPINKPEPTIDE".into();
        let digest = DigestSlice::from_string(pepseq, false, 1);
        let query = eg;
        let expected_intensity = ei;
        QueryItemToScore {
            digest,
            query,
            expected_intensity,
        }
    }
}
