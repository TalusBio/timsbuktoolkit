use crate::tinyvec::{
    TinyVec,
    tiny_vec,
};
use crate::{
    KeyLike,
    TimsElutionGroup,
};

#[derive(Debug)]
pub enum ElutionGroupInputError {
    MismatchedFragmentLabelsLength {
        expected: usize,
        found: usize,
    },
    IonConversionError {
        inner: String,
    },
    /// The input shipped no `fragment_labels`. There is no way to synthesize
    /// them: a label names an ion series and ordinal, and a positional index
    /// carries no chemistry.
    MissingFragmentLabels,
}

/// User-friendly format for specifying elution groups in an input file
/// (copied from timsquery_cli for compatibility)
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ElutionGroupInput<T: KeyLike> {
    pub id: u64,
    pub mobility: f32,
    pub rt_seconds: f32,
    #[serde(alias = "precursor_mz")]
    #[serde(alias = "precursor_mono_mz")]
    pub precursor: f64,
    pub precursor_charge: u8,
    pub precursor_isotopes: Option<Vec<i8>>,
    #[serde(alias = "fragment_mzs")]
    #[serde(alias = "fragment_mono_mzs")]
    pub fragments: Vec<f64>,
    pub fragment_labels: Option<Vec<T>>,
}

impl<T: KeyLike, U: TryInto<T> + KeyLike> TryFrom<ElutionGroupInput<U>> for TimsElutionGroup<T> {
    type Error = ElutionGroupInputError;

    fn try_from(val: ElutionGroupInput<U>) -> Result<Self, Self::Error> {
        let builder = TimsElutionGroup::builder()
            .id(val.id)
            .mobility_ook0(val.mobility)
            .rt_seconds(val.rt_seconds)
            .precursor(val.precursor, val.precursor_charge)
            .precursor_labels(if let Some(isotopes) = val.precursor_isotopes {
                TinyVec::Heap(isotopes.into_iter().collect())
            } else {
                tiny_vec!(0i8)
            });

        let num_fragments = val.fragments.len();
        let builder = builder.fragment_mzs(val.fragments.clone());
        let builder = match val.fragment_labels {
            Some(ref labels) => {
                if labels.len() != num_fragments {
                    return Err(ElutionGroupInputError::MismatchedFragmentLabelsLength {
                        expected: num_fragments,
                        found: labels.len(),
                    });
                }
                let labels = labels
                    .iter()
                    .map(|x| {
                        x.clone().try_into().map_err(|_e| {
                            ElutionGroupInputError::IonConversionError {
                                inner: format!("Failed to convert fragment label: {:?}", x),
                            }
                        })
                    })
                    .collect::<Result<Vec<T>, ElutionGroupInputError>>()?;
                builder.fragment_labels(TinyVec::Heap(labels))
            }
            None => return Err(ElutionGroupInputError::MissingFragmentLabels),
        };

        Ok(builder.try_build().expect("I checked the sizes!"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ion::IonAnnot;

    fn input(fragment_labels: Option<Vec<IonAnnot>>) -> ElutionGroupInput<IonAnnot> {
        ElutionGroupInput {
            id: 0,
            mobility: 0.8,
            rt_seconds: 100.0,
            precursor: 500.0,
            precursor_charge: 2,
            precursor_isotopes: None,
            fragments: vec![100.0, 200.0],
            fragment_labels,
        }
    }

    /// Unlabelled input is rejected by name rather than by whatever fails
    /// first downstream. Nothing can stand in for a missing label.
    #[test]
    fn missing_fragment_labels_is_its_own_error() {
        let err = TimsElutionGroup::<IonAnnot>::try_from(input(None))
            .expect_err("no labels means no elution group");
        assert!(
            matches!(err, ElutionGroupInputError::MissingFragmentLabels),
            "expected MissingFragmentLabels, got {err:?}"
        );
    }

    #[test]
    fn label_count_must_match_fragment_count() {
        let one_label = vec![IonAnnot::try_from("y1").unwrap()];
        let err = TimsElutionGroup::<IonAnnot>::try_from(input(Some(one_label)))
            .expect_err("1 label for 2 fragments");
        assert!(
            matches!(
                err,
                ElutionGroupInputError::MismatchedFragmentLabelsLength {
                    expected: 2,
                    found: 1
                }
            ),
            "expected MismatchedFragmentLabelsLength, got {err:?}"
        );
    }
}
