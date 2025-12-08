use crate::ion::IonAnnot;
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
    MismatchedFragmentLabelsLength { expected: usize, found: usize },
    AlreadyHasFragmentLabels,
    IonConversionError { inner: String },
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

impl<T: KeyLike> ElutionGroupInput<T> {
    pub fn needs_fragment_labels(&self) -> bool {
        self.fragment_labels.is_none()
    }

    pub fn try_fill_labels_u8(self) -> Result<ElutionGroupInput<u8>, ElutionGroupInputError> {
        let num_fragments = self.fragments.len();
        if let Some(_) = self.fragment_labels {
            return Err(ElutionGroupInputError::AlreadyHasFragmentLabels);
        }
        let fragment_labels: Vec<u8> = (0..num_fragments).map(|i| i as u8).collect();

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

    pub fn try_fill_labels_annot(
        self,
    ) -> Result<ElutionGroupInput<IonAnnot>, ElutionGroupInputError> {
        let tmp = self.try_fill_labels_u8()?;
        let new_frags = tmp
            .fragment_labels
            .unwrap()
            .into_iter()
            .map(|lbl| IonAnnot::try_new('?', Some(lbl), 1, 1).unwrap())
            .collect();
        Ok(ElutionGroupInput {
            id: tmp.id,
            mobility: tmp.mobility,
            rt_seconds: tmp.rt_seconds,
            precursor: tmp.precursor,
            precursor_charge: tmp.precursor_charge,
            precursor_isotopes: tmp.precursor_isotopes,
            fragments: tmp.fragments,
            fragment_labels: Some(new_frags),
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
