use serde::Serialize;
use serde::ser::SerializeStruct;

use crate::models::frames::peak_in_quad::ResolvedPeakInQuad;
use crate::utils::streaming_calculators::{
    RunningStatsCalculator,
    StreamingAggregatorError,
};
use crate::{
    ElutionGroup,
    KeyLike,
    ValueLike,
};
use std::ops::{Add, AddAssign};
use std::sync::Arc;

#[derive(Debug, Clone, Serialize)]
pub struct SpectralCollector<T: KeyLike, V: Default + ValueLike> {
    pub eg: Arc<ElutionGroup<T>>,
    precursors: Vec<V>,
    fragments: Vec<V>,
}

impl<T: KeyLike, V: ValueLike + Default> SpectralCollector<T, V> {
    pub fn new(eg: Arc<ElutionGroup<T>>) -> Self {
        let precursors = vec![Default::default(); eg.precursors.len()];
        let fragments = vec![Default::default(); eg.fragments.len()];
        Self {
            eg,
            precursors,
            fragments,
        }
    }

    pub fn iter_mut_precursors(&mut self) -> impl Iterator<Item = (&(i8, f64), &mut V)> {
        self.eg.precursors.iter().zip(self.precursors.iter_mut())
    }

    pub fn iter_mut_fragments(&mut self) -> impl Iterator<Item = (&(T, f64), &mut V)> {
        self.eg.fragments.iter().zip(self.fragments.iter_mut())
    }

    pub fn iter_precursors(&self) -> impl Iterator<Item = (&(i8, f64), &V)> {
        self.eg.precursors.iter().zip(self.precursors.iter())
    }

    pub fn iter_fragments(&self) -> impl Iterator<Item = (&(T, f64), &V)> {
        self.eg.fragments.iter().zip(self.fragments.iter())
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
        match (self.mean_mz, self.mean_mobility) {
            (Some(mut mmz), Some(mut mmob)) => {
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

impl Add for MzMobilityStatsCollector {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        let mz = match (self.mean_mz, other.mean_mz) {
            (Some(mut mmz), Some(ommz)) => {
                mmz.add(ommz.weight(), ommz.mean().expect("mean should be initialized with value"));
                Some(mmz)
            }
            (None, None) => None,
            (Some(x), None) => Some(x),
            (None, Some(x)) => Some(x),
        };

        let mobility = match (self.mean_mobility, other.mean_mobility) {
            (Some(mut mmob), Some(ommob)) => {
            mmob.add(ommob.weight(), ommob.mean().expect("mean should be initialized with value"));
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


impl AddAssign<ResolvedPeakInQuad> for MzMobilityStatsCollector {
    fn add_assign(&mut self, other: ResolvedPeakInQuad) {
        self.add(
            other.corrected_intensity as f64,
            other.mz as f64,
            other.mobility as f64,
        );
    }
}
