use crate::models::elution_group::tims_elution_group_builder::{
    SetPrecursorCharge,
    SetPrecursorMonoMz,
};
use crate::traits::KeyLike;
use crate::utils::constants::NEUTRON_MASS;
use bon::builder;
use serde::{
    Deserialize,
    Serialize,
};
use tinyvec::TinyVec;

/// A struct that represents an elution group.
///
/// The elution group is a single precursor ion that is framented.
/// The fragments m/z values are stored in a Vec where each element
/// contains the ion annotation and its m/z value.
#[derive(Debug, Serialize, Deserialize, Clone, bon::Builder)]
#[builder(finish_fn(vis = "", name = try_build_internal))]
pub struct TimsElutionGroup<T: KeyLike> {
    id: u64,
    mobility_ook0: f32,
    rt_seconds: f32,
    precursor_mono_mz: f64,
    precursor_charge: u8,

    // The baseline size of TinyVec<T> is 24 bits. (due to the stack size of a Vec<T>)
    // NOTE: This ca be different on different architectures. (would be 16 bytes on 32-bit
    // architectures)
    //
    // which means that when the array breaks even with the size of the array
    // (minus some space for the enum itself ...)
    // For IonAnnot which is 4 bytes, capacity of 3 == 24 bytes, but capacity of 4 == 32 bytes
    // and remains the same up to 7 (inclusive) (13 == 56 bytes)
    //
    // For an f64 which is 8 bytes, capacity of 1 == 24 bytes, but capacity of 2 == 32 bytes
    //
    // For an i8 [i8; 14] bunps is to 32
    //
    // Implementing all 4 fields as Vec<T> makes this struct 112 bytes on 64-bit architectures
    // regardless of the types (bc those are just fat pointers)
    //
    // Changing the labels to TinyVec with capaciry of 13 makes the struct 144 bytes for
    // the concrete type TimsElutionGroup<IonAnnot> but 408 bytes for TimsElutionGroup<String>
    //
    // In theory I can make this lighter if it was a genetic ...
    fragment_mzs: Vec<f64>,
    fragment_labels: TinyVec<[T; 13]>,
    precursor_labels: TinyVec<[i8; 13]>,
}

impl<T: KeyLike + Default, S: tims_elution_group_builder::IsComplete>
    TimsElutionGroupBuilder<T, S>
{
    pub fn try_build(self) -> Result<TimsElutionGroup<T>, crate::errors::DataProcessingError> {
        let candidate = self.try_build_internal();
        if candidate.fragment_labels.is_empty() | candidate.precursor_labels.is_empty() {
            return Err(crate::DataProcessingError::ExpectedNonEmptyData);
        }
        if candidate.fragment_mzs.len() != candidate.fragment_labels.len() {
            return Err(crate::DataProcessingError::ExpectedVectorSameLength);
        }

        Ok(candidate)
    }
}

impl<T: KeyLike + Default, S: tims_elution_group_builder::State> TimsElutionGroupBuilder<T, S> {
    pub fn precursor(
        self,
        mz_mono: f64,
        charge: u8,
    ) -> TimsElutionGroupBuilder<T, SetPrecursorCharge<SetPrecursorMonoMz<S>>>
    where
        S::PrecursorCharge: tims_elution_group_builder::IsUnset,
        S::PrecursorMonoMz: tims_elution_group_builder::IsUnset,
    {
        self.precursor_mono_mz(mz_mono).precursor_charge(charge)
    }
}

impl<T: KeyLike> TimsElutionGroup<T> {
    pub fn id(&self) -> u64 {
        self.id
    }

    pub fn precursor_count(&self) -> usize {
        self.precursor_labels.len()
    }

    pub fn fragment_count(&self) -> usize {
        self.fragment_mzs.len()
    }

    pub fn rt_seconds(&self) -> f32 {
        self.rt_seconds
    }

    pub fn precursor_charge(&self) -> u8 {
        self.precursor_charge
    }

    pub fn precursor_mz(&self) -> f64 {
        self.precursor_mono_mz
    }

    pub fn set_precursor_labels(&mut self, labels: impl Iterator<Item = i8>) {
        self.precursor_labels.clear();
        self.precursor_labels.extend(labels);
    }

    // NOTE: I am thinking about removing this and leave the rest as a trait
    pub fn set_rt_seconds(&mut self, rt_seconds: f32) {
        self.rt_seconds = rt_seconds;
    }

    pub fn mobility_ook0(&self) -> f32 {
        self.mobility_ook0
    }

    pub fn try_drop_fragment(&mut self, key: &T) -> Result<(), crate::errors::DataProcessingError> {
        let idx = self
            .fragment_labels
            .iter()
            .position(|label| label == key)
            .ok_or(crate::DataProcessingError::KeyNotFound)?;

        self.fragment_labels.remove(idx);
        self.fragment_mzs.remove(idx);
        Ok(())
    }

    pub fn try_drop_precursor(
        &mut self,
        key: i8,
    ) -> Result<(), crate::errors::DataProcessingError> {
        let idx = self
            .precursor_labels
            .iter()
            .position(|label| *label == key)
            .ok_or(crate::DataProcessingError::KeyNotFound)?;

        self.precursor_labels.remove(idx);
        Ok(())
    }

    fn precursor_mz_iter(&self) -> impl Iterator<Item = f64> + '_ {
        let offset = NEUTRON_MASS / self.precursor_charge as f64;
        self.precursor_labels.iter().map(move |isotope_index| {
            let idx = *isotope_index;
            self.precursor_mono_mz + (offset * idx as f64)
        })
    }

    pub fn get_precursor_mz_limits(&self) -> (f64, f64) {
        let mut min_precursor_mz = f64::MAX;
        let mut max_precursor_mz = f64::MIN;
        for (isotope_index, precursor_mz) in
            self.precursor_labels.iter().zip(self.precursor_mz_iter())
        {
            if *isotope_index < 0 {
                continue;
            }
            if precursor_mz < min_precursor_mz {
                min_precursor_mz = precursor_mz;
            }
            if precursor_mz > max_precursor_mz {
                max_precursor_mz = precursor_mz;
            }
        }
        (min_precursor_mz, max_precursor_mz)
    }

    pub fn mono_precursor_mz(&self) -> f64 {
        self.precursor_mono_mz
    }

    pub fn iter_precursors(&self) -> impl Iterator<Item = (i8, f64)> {
        self.precursor_labels
            .iter()
            .zip(self.precursor_mz_iter())
            .map(|(label, mz)| (*label, mz))
    }

    pub fn iter_fragments_refs(&self) -> impl Iterator<Item = (&T, &f64)> {
        self.fragment_labels.iter().zip(self.fragment_mzs.iter())
    }

    pub fn iter_fragments_refs_mut(&mut self) -> impl Iterator<Item = (&mut T, &mut f64)> {
        self.fragment_labels
            .iter_mut()
            .zip(self.fragment_mzs.iter_mut())
    }

    pub fn iter_fragments(&self) -> impl Iterator<Item = (&T, f64)> {
        self.fragment_labels
            .iter()
            .zip(self.fragment_mzs.iter())
            .map(|(label, mz)| (label, *mz))
    }

    pub fn with_rt_seconds(self, rt_seconds: f32) -> Self {
        Self { rt_seconds, ..self }
    }

    pub fn with_mobility(self, mobility: f32) -> Self {
        Self {
            mobility_ook0: mobility,
            ..self
        }
    }

    pub fn cast<U: KeyLike>(&self, f: impl Fn(&T) -> U) -> TimsElutionGroup<U> {
        let fragment_labels_converted: TinyVec<[U; 13]> =
            self.fragment_labels.iter().map(f).collect();

        TimsElutionGroup {
            id: self.id,
            mobility_ook0: self.mobility_ook0,
            rt_seconds: self.rt_seconds,
            precursor_mono_mz: self.precursor_mono_mz,
            precursor_charge: self.precursor_charge,
            fragment_mzs: self.fragment_mzs.clone(),
            fragment_labels: fragment_labels_converted,
            precursor_labels: self.precursor_labels.clone(),
        }
    }
}
