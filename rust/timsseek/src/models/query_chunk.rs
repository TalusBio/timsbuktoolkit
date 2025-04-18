use super::digest::DigestSlice;
use crate::data_sources::speclib::ExpectedIntensities;
use crate::fragment_mass::IonAnnot;
use rayon::iter::Zip as RayonZip;
use rayon::prelude::*;
use rayon::vec::IntoIter as RayonVecIntoIter;
use std::sync::Arc;
use timsquery::models::elution_group::ElutionGroup;

#[derive(Debug, Clone)]
pub struct NamedQueryChunk {
    pub digests: Vec<DigestSlice>,
    pub charges: Vec<u8>,
    pub queries: Vec<Arc<ElutionGroup<IonAnnot>>>,
    pub expected_intensities: Vec<ExpectedIntensities>,
}

impl NamedQueryChunk {
    pub fn new(
        digests: Vec<DigestSlice>,
        charges: Vec<u8>,
        queries: Vec<Arc<ElutionGroup<IonAnnot>>>,
        expected_intensities: Vec<ExpectedIntensities>,
    ) -> Self {
        assert_eq!(digests.len(), charges.len());
        assert_eq!(digests.len(), queries.len());
        assert_eq!(digests.len(), expected_intensities.len());
        Self {
            digests,
            charges,
            queries,
            expected_intensities,
        }
    }

    pub fn into_zip_par_iter(
        self,
    ) -> RayonZip<
        RayonVecIntoIter<ExpectedIntensities>,
        RayonZip<
            RayonVecIntoIter<Arc<ElutionGroup<IonAnnot>>>,
            RayonZip<RayonVecIntoIter<DigestSlice>, RayonVecIntoIter<u8>>,
        >,
    > {
        // IN THEORY I should implement IntoIter for this struct
        // but I failed at it (skill issues?) so this will do for now.
        // JSPP - 2024-11-21
        self.expected_intensities.into_par_iter().zip(
            self.queries.into_par_iter().zip(
                self.digests
                    .into_par_iter()
                    .zip(self.charges.into_par_iter()),
            ),
        )
    }

    pub fn len(&self) -> usize {
        self.queries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.queries.is_empty()
    }
}
