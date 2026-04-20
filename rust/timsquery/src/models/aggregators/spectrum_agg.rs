use serde::Serialize;
use serde::ser::SerializeStruct;
use timscentroid::indexing::IndexedPeak;
use timscentroid::rt_mapping::RTIndex;
use tinyvec::TinyVec;

use crate::traits::queriable_data::{
    HasQueryData,
    PeakAddable,
};
use crate::utils::streaming_calculators::{
    RunningStatsCalculator,
    StreamingAggregatorError,
};
use crate::{
    KeyLike,
    TimsElutionGroup,
    ValueLike,
};
use std::ops::{
    Add,
    AddAssign,
};

/// Inline-capacity target matching `TimsElutionGroup`'s precursor/fragment
/// label TinyVecs — typical peptide ≤13 fragments / ≤3 precursors stays
/// stack-resident.
const SPEC_INLINE_CAP: usize = 13;

// TODO: rename to `SpectralAccumulator` — struct holds query scalars +
// label/mz lists alongside the accumulated intensities. "Collector" name
// predates the Query/Accumulator split.
#[derive(Debug, Clone, Serialize)]
pub struct SpectralCollector<T: KeyLike, V: Default + ValueLike> {
    // Query scalars carried from eg at construction / reset.
    pub id: u64,
    pub mobility_ook0: f32,
    pub rt_seconds: f32,
    pub precursor_mono_mz: f64,
    pub precursor_charge: u8,
    /// Cached from `TimsElutionGroup::get_precursor_mz_limits()` — skips
    /// negative-isotope labels, do NOT recompute from mono_mz + charge.
    pub precursor_mz_limits: (f64, f64),
    // Labels + mzs: arrays carry only intensities, so we need separate storage.
    pub precursor_labels: TinyVec<[i8; SPEC_INLINE_CAP]>,
    pub precursor_mzs: TinyVec<[f64; SPEC_INLINE_CAP]>,
    pub fragment_labels: TinyVec<[T; SPEC_INLINE_CAP]>,
    pub fragment_mzs: TinyVec<[f64; SPEC_INLINE_CAP]>,
    precursors: Vec<V>,
    fragments: Vec<V>,
}

impl<T: KeyLike, V: ValueLike + Default> SpectralCollector<T, V> {
    pub fn new(eg: &TimsElutionGroup<T>) -> Self {
        let mut out = Self {
            id: 0,
            mobility_ook0: 0.0,
            rt_seconds: 0.0,
            precursor_mono_mz: 0.0,
            precursor_charge: 0,
            precursor_mz_limits: (0.0, 0.0),
            precursor_labels: TinyVec::new(),
            precursor_mzs: TinyVec::new(),
            fragment_labels: TinyVec::new(),
            fragment_mzs: TinyVec::new(),
            precursors: Vec::new(),
            fragments: Vec::new(),
        };
        out.reset_with_overrides(eg, None, None);
        out
    }

    pub fn reset_with(&mut self, eg: &TimsElutionGroup<T>) {
        self.reset_with_overrides(eg, None, None);
    }

    /// Reset with optional rt/mobility overrides — replaces
    /// `item.query.clone().with_rt_seconds(r).with_mobility(m)` at callers.
    pub fn reset_with_overrides(
        &mut self,
        eg: &TimsElutionGroup<T>,
        rt_override: Option<f32>,
        mobility_override: Option<f32>,
    ) {
        self.id = eg.id();
        self.mobility_ook0 = mobility_override.unwrap_or_else(|| eg.mobility_ook0());
        self.rt_seconds = rt_override.unwrap_or_else(|| eg.rt_seconds());
        self.precursor_mono_mz = eg.mono_precursor_mz();
        self.precursor_charge = eg.precursor_charge();
        self.precursor_mz_limits = eg.get_precursor_mz_limits();

        self.precursor_labels.clear();
        self.precursor_mzs.clear();
        for (lbl, mz) in eg.iter_precursors() {
            self.precursor_labels.push(lbl);
            self.precursor_mzs.push(mz);
        }

        self.fragment_labels.clear();
        self.fragment_mzs.clear();
        for (lbl, mz) in eg.iter_fragments_refs() {
            self.fragment_labels.push(lbl.clone());
            self.fragment_mzs.push(*mz);
        }

        // Intensity buffers reuse capacity via resize.
        self.precursors.clear();
        self.precursors
            .resize(self.precursor_labels.len(), V::default());
        self.fragments.clear();
        self.fragments
            .resize(self.fragment_labels.len(), V::default());
    }

    pub fn iter_mut_precursors(&mut self) -> impl Iterator<Item = ((i8, f64), &mut V)> {
        self.precursor_labels
            .iter()
            .copied()
            .zip(self.precursor_mzs.iter().copied())
            .zip(self.precursors.iter_mut())
    }

    pub fn iter_mut_fragments(&mut self) -> impl Iterator<Item = ((&T, &f64), &mut V)> {
        self.fragment_labels
            .iter()
            .zip(self.fragment_mzs.iter())
            .zip(self.fragments.iter_mut())
    }

    pub fn iter_precursors(&self) -> impl Iterator<Item = ((i8, f64), &V)> {
        self.precursor_labels
            .iter()
            .copied()
            .zip(self.precursor_mzs.iter().copied())
            .zip(self.precursors.iter())
    }

    pub fn iter_fragments(&self) -> impl Iterator<Item = ((&T, &f64), &V)> {
        self.fragment_labels
            .iter()
            .zip(self.fragment_mzs.iter())
            .zip(self.fragments.iter())
    }
}

impl<T: KeyLike, V: Default + ValueLike> HasQueryData<T> for SpectralCollector<T, V> {
    fn id(&self) -> u64 {
        self.id
    }

    fn precursor_mz_limits(&self) -> (f64, f64) {
        self.precursor_mz_limits
    }

    fn mobility_ook0(&self) -> f32 {
        self.mobility_ook0
    }

    fn rt_seconds(&self) -> f32 {
        self.rt_seconds
    }

    fn iter_precursors(&self) -> impl Iterator<Item = (i8, f64)> + '_ {
        self.precursor_labels
            .iter()
            .copied()
            .zip(self.precursor_mzs.iter().copied())
    }

    fn iter_fragments<'a>(&'a self) -> impl Iterator<Item = (&'a T, f64)> + 'a
    where
        T: 'a,
    {
        self.fragment_labels
            .iter()
            .zip(self.fragment_mzs.iter().copied())
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct MzMobilityStatsCollector {
    // TODO: Make this guy smaller ... this is a surpirsingly
    // large struct ...
    pub mean_mz: Option<RunningStatsCalculator>,
    pub mean_mobility: Option<RunningStatsCalculator>,
}

impl Serialize for MzMobilityStatsCollector {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let mut state = serializer.serialize_struct("MzMobilityStatsCollector", 2)?;
        match (self.mean_mz(), self.mean_mobility()) {
            (Ok(mz), Ok(mob)) => {
                state.serialize_field("mean_mz", &mz)?;
                state.serialize_field("mean_mobility", &mob)?;
            }
            (Err(_), Err(_)) => {
                state.serialize_field("mean_mz", &f64::NAN)?;
                state.serialize_field("mean_mobility", &f64::NAN)?;
            }
            _ => unreachable!(),
        }
        state.end()
    }
}

impl MzMobilityStatsCollector {
    pub fn new(weight: f64, mean_mz: f64, mean_mobility: f64) -> Self {
        Self {
            mean_mz: Some(RunningStatsCalculator::new(weight, mean_mz)),
            mean_mobility: Some(RunningStatsCalculator::new(weight, mean_mobility)),
        }
    }

    pub fn add(&mut self, weight: f64, mz: f64, mobility: f64) {
        match (self.mean_mz.as_mut(), self.mean_mobility.as_mut()) {
            (Some(mmz), Some(mmob)) => {
                mmz.add(weight, mz);
                mmob.add(weight, mobility);
            }
            (None, None) => {
                self.mean_mz = Some(RunningStatsCalculator::new(weight, mz));
                self.mean_mobility = Some(RunningStatsCalculator::new(weight, mobility));
            }
            _ => unreachable!(),
        }
    }

    pub fn mean_mz(&self) -> Result<f64, StreamingAggregatorError> {
        self.mean_mz
            .map_or(Err(StreamingAggregatorError::NotEnoughData), |x| x.mean())
    }

    pub fn mean_mobility(&self) -> Result<f64, StreamingAggregatorError> {
        self.mean_mobility
            .map_or(Err(StreamingAggregatorError::NotEnoughData), |x| x.mean())
    }

    pub fn weight(&self) -> f64 {
        self.mean_mz.as_ref().map_or(0.0, |x| x.weight())
    }
}

impl<T: RTIndex> AddAssign<IndexedPeak<T>> for MzMobilityStatsCollector {
    fn add_assign(&mut self, other: IndexedPeak<T>) {
        self.add(
            other.intensity as f64,
            other.mz as f64,
            other.mobility_ook0.to_f64(),
        );
    }
}

impl<T: RTIndex> PeakAddable<T> for MzMobilityStatsCollector {}

impl Add for MzMobilityStatsCollector {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let mz = match (self.mean_mz, other.mean_mz) {
            (Some(mut mmz), Some(ommz)) => {
                mmz.add(
                    ommz.weight(),
                    ommz.mean().expect("mean should be initialized with value"),
                );
                Some(mmz)
            }
            (None, None) => None,
            (Some(x), None) => Some(x),
            (None, Some(x)) => Some(x),
        };

        let mobility = match (self.mean_mobility, other.mean_mobility) {
            (Some(mut mmob), Some(ommob)) => {
                mmob.add(
                    ommob.weight(),
                    ommob.mean().expect("mean should be initialized with value"),
                );
                Some(mmob)
            }
            (None, None) => None,
            (Some(x), None) => Some(x),
            (None, Some(x)) => Some(x),
        };

        Self {
            mean_mz: mz,
            mean_mobility: mobility,
        }
    }
}
