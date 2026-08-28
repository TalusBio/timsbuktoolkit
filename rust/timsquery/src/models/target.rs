use crate::models::target::target_builder::{
    SetPrecursorCharge,
    SetPrecursorMonoMz,
};
use crate::traits::KeyLike;
use crate::utils::constants::C13_C12_MASS_DIFF;
use serde::{
    Deserialize,
    Serialize,
};
use tinyvec::TinyVec;

/// One analyte's constraints in (m/z, RT, mobility) space, plus the fragment
/// m/z values to extract. The row form, as opposed to the columnar
/// [`crate::models::TargetColumns`].
#[derive(Debug, Serialize, Deserialize, Clone, bon::Builder)]
#[builder(finish_fn(vis = "", name = try_build_internal))]
pub struct Target<T: KeyLike> {
    #[builder(into)]
    id: crate::models::OwnedSourceId,
    #[serde(alias = "mobility")]
    mobility_ook0: f32,
    rt_seconds: f32,
    #[serde(alias = "precursor")]
    #[serde(alias = "precursor_mz")]
    precursor_mono_mz: f64,
    #[serde(alias = "charge")]
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
    // the concrete type Target<IonAnnot> but 408 bytes for Target<String>
    //
    // In theory I can make this lighter if it was a genetic ...
    #[serde(alias = "fragments")]
    fragment_mzs: Vec<f64>,
    fragment_labels: TinyVec<[T; 13]>,
    #[serde(alias = "precursor_isotopes")]
    precursor_labels: TinyVec<[i8; 13]>,
}

impl<T: KeyLike + Default, S: target_builder::IsComplete> TargetBuilder<T, S> {
    pub fn try_build(self) -> Result<Target<T>, crate::errors::DataProcessingError> {
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

impl<T: KeyLike + Default, S: target_builder::State> TargetBuilder<T, S> {
    pub fn precursor(
        self,
        mz_mono: f64,
        charge: u8,
    ) -> TargetBuilder<T, SetPrecursorCharge<SetPrecursorMonoMz<S>>>
    where
        S::PrecursorCharge: target_builder::IsUnset,
        S::PrecursorMonoMz: target_builder::IsUnset,
    {
        self.precursor_mono_mz(mz_mono).precursor_charge(charge)
    }
}

impl<T: KeyLike + Default> Target<T> {
    /// Build an empty scratch group DIRECTLY (bypassing the `bon` builder,
    /// which rejects empty fragment/precursor label sets). All scalars are
    /// zeroed and all vecs empty; the intended use is a reuse-in-place scratch
    /// buffer that gets refilled via [`Target::reset_from`] before
    /// every query, so the zeroed initial state is never observed.
    pub fn empty_like() -> Self {
        Self {
            id: Default::default(),
            mobility_ook0: 0.0,
            rt_seconds: 0.0,
            precursor_mono_mz: 0.0,
            precursor_charge: 0,
            fragment_mzs: Vec::new(),
            fragment_labels: TinyVec::new(),
            precursor_labels: TinyVec::new(),
        }
    }
}

impl<T: KeyLike> Target<T> {
    pub fn id(&self) -> crate::models::SourceId<'_> {
        self.id.as_ref()
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

    /// In-place copy reusing Vec/TinyVec capacity. The `clear()` + `push`
    /// pattern preserves the destination's heap buffer capacity across
    /// resets — after warm-up, zero alloc. Used by the isotope-offset
    /// scratch in timsseek. `G` is any `QueryGeom` (e.g. the columnar
    /// flyweight), not necessarily `Self` — the body reads `src` only
    /// through trait methods.
    pub fn reset_from<G: crate::traits::QueryGeom<Label = T>>(&mut self, src: &G) {
        self.id.set_from(src.output_id());
        self.mobility_ook0 = src.mobility_ook0();
        self.rt_seconds = src.rt_seconds();
        self.precursor_mono_mz = src.mono_precursor_mz();
        self.precursor_charge = src.precursor_charge();
        self.fragment_mzs.clear();
        self.fragment_labels.clear();
        for (lab, mz) in src.iter_fragments_refs() {
            // mz owned + shifted for decoy flyweights
            self.fragment_labels.push(lab.clone());
            self.fragment_mzs.push(mz);
        }
        self.precursor_labels.clear();
        self.precursor_labels
            .extend(src.iter_precursors().map(|(iso, _mz)| iso));
    }

    pub fn mobility_ook0(&self) -> f32 {
        self.mobility_ook0
    }

    fn precursor_mz_iter(&self) -> impl Iterator<Item = f64> + '_ {
        let offset = C13_C12_MASS_DIFF / self.precursor_charge as f64;
        self.precursor_labels.iter().map(move |isotope_index| {
            let idx = *isotope_index;
            self.precursor_mono_mz + (offset * idx as f64)
        })
    }

    pub fn precursor_mz_limits(&self) -> (f64, f64) {
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

    pub fn cast<U: KeyLike>(&self, f: impl Fn(&T) -> U) -> Target<U> {
        let fragment_labels_converted: TinyVec<[U; 13]> =
            self.fragment_labels.iter().map(f).collect();

        Target {
            id: self.id.clone(),
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
