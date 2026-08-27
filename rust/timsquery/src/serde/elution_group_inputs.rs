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
    AlreadyHasFragmentLabels,
    IonConversionError {
        inner: String,
    },
    MissingFragmentLabels,
    /// More fragments than a `u8` label space can index uniquely.
    TooManyFragmentsToLabel {
        count: usize,
    },
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

impl<T: KeyLike> ElutionGroupInput<T> {
    pub fn needs_fragment_labels(&self) -> bool {
        self.fragment_labels.is_none()
    }

    /// Synthesize positional `u8` labels for an input that shipped none.
    ///
    /// Labels must be unique within the elution group (see
    /// `ExpectedIntensities::try_from_pairs`), and `u8` can only index 256 of
    /// them. `i as u8` would wrap silently past that — emitting `0,1,..,255,0,1,..`
    /// and failing far downstream with a duplicate-key error — so the overflow
    /// is rejected here instead.
    pub fn try_fill_labels_u8(self) -> Result<ElutionGroupInput<u8>, ElutionGroupInputError> {
        let num_fragments = self.fragments.len();
        if self.fragment_labels.is_some() {
            return Err(ElutionGroupInputError::AlreadyHasFragmentLabels);
        }
        let fragment_labels: Vec<u8> = (0..num_fragments)
            .map(u8::try_from)
            .collect::<Result<_, _>>()
            .map_err(|_| ElutionGroupInputError::TooManyFragmentsToLabel {
                count: num_fragments,
            })?;

        Ok(ElutionGroupInput {
            id: self.id,
            mobility: self.mobility,
            rt_seconds: self.rt_seconds,
            precursor: self.precursor,
            precursor_charge: self.precursor_charge,
            precursor_isotopes: self.precursor_isotopes,
            fragments: self.fragments,
            fragment_labels: Some(fragment_labels),
        })
    }
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

    fn input_with_n_fragments(n: usize) -> ElutionGroupInput<u8> {
        ElutionGroupInput {
            id: 0,
            mobility: 0.8,
            rt_seconds: 100.0,
            precursor: 500.0,
            precursor_charge: 2,
            precursor_isotopes: None,
            fragments: vec![100.0; n],
            fragment_labels: None,
        }
    }

    /// Synthesized labels must be unique — that is the invariant
    /// `ExpectedIntensities::try_from_pairs` relies on. 256 fragments is the
    /// exact capacity of the `u8` label space.
    #[test]
    fn fill_labels_u8_is_unique_at_capacity() {
        let filled = input_with_n_fragments(256)
            .try_fill_labels_u8()
            .expect("256 fragments fit the u8 label space");
        let labels = filled.fragment_labels.unwrap();
        assert_eq!(labels.len(), 256);
        let unique: std::collections::HashSet<_> = labels.iter().copied().collect();
        assert_eq!(unique.len(), 256, "synthesized labels must all be distinct");
    }

    /// One past capacity must fail rather than wrap: a second `0` label would
    /// collide with the first and corrupt scoring downstream.
    #[test]
    fn fill_labels_u8_rejects_overflow_instead_of_wrapping() {
        let err = input_with_n_fragments(257)
            .try_fill_labels_u8()
            .expect_err("257 fragments cannot be labelled uniquely with a u8");
        assert!(
            matches!(
                err,
                ElutionGroupInputError::TooManyFragmentsToLabel { count: 257 }
            ),
            "expected TooManyFragmentsToLabel, got {err:?}"
        );
    }
}
