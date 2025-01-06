use super::super::chromatogram_agg::ScanTofStatsCalculatorPair;
use super::aggregator::ParitionedCMGAggregator;
use super::arrays::{
    PartitionedCMGArrayStats,
    PartitionedCMGArrays,
};
use crate::TimsqueryError;
use crate::errors::Result;
use crate::models::elution_group::ElutionGroup;
use crate::models::frames::raw_peak::RawPeak;
use crate::models::queries::MsLevelContext;
use crate::traits::aggregator::Aggregator;
use serde::Serialize;
use std::hash::Hash;

use timsrust::converters::{
    ConvertableDomain,
    Scan2ImConverter,
    Tof2MzConverter,
};

#[derive(Debug, Clone)]
pub struct MultiCMGStatsAgg<FH: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug> {
    pub converters: (Tof2MzConverter, Scan2ImConverter),
    pub ms1_stats: ParitionedCMGAggregator<usize>,
    pub ms2_stats: ParitionedCMGAggregator<FH>,
    pub id: u64,
    pub context: Option<MsLevelContext<usize, FH>>,
    pub buffer: Option<ScanTofStatsCalculatorPair>,
    // TODO: Make this a reference instead of a clone.... maybe.
    pub elution_group_ref: ElutionGroup<FH>,
}

#[derive(Debug, Clone)]
pub struct MultiCMGStatsFactory<FH: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug> {
    pub converters: (Tof2MzConverter, Scan2ImConverter),
    pub _phantom: std::marker::PhantomData<FH>,
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug> MultiCMGStatsFactory<FH> {
    pub fn build_with_elution_group(
        &self,
        elution_group: &ElutionGroup<FH>,
    ) -> MultiCMGStatsAgg<FH> {
        // TODO: RN this is a super hacky ... IDEALLY we would either keep a reference OR
        // preserve only the expected intensities.
        let elution_group_ref = elution_group.clone();

        let (fragment_keys, fragment_tof): (Vec<FH>, Vec<u32>) = elution_group
            .fragment_mzs
            .iter()
            .map(|(k, v)| {
                let tof = self.converters.0.invert(*v) as u32;
                (k.clone(), tof)
            })
            .unzip();

        let (precursor_keys, precursor_tof): (Vec<usize>, Vec<u32>) = elution_group
            .precursor_mzs
            .iter()
            .enumerate()
            .map(|(i, v)| (i, self.converters.0.invert(*v) as u32))
            .unzip();

        let expect_scan_index = self.converters.1.invert(elution_group.mobility) as usize;

        MultiCMGStatsAgg {
            converters: (self.converters.0, self.converters.1),
            ms1_stats: ParitionedCMGAggregator::new(
                precursor_keys,
                expect_scan_index,
                precursor_tof,
            ),
            ms2_stats: ParitionedCMGAggregator::new(fragment_keys, expect_scan_index, fragment_tof),
            id: elution_group.id,
            context: None,
            buffer: None,
            elution_group_ref,
        }
    }
}

// TODO: I hate this name ... I need to find something
// better - JSPP - 2024-12-20
#[derive(Debug, Clone, Serialize)]
pub struct NaturalFinalizedMultiCMGArrays<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub ms1_arrays: PartitionedCMGArrayStats<usize>,
    pub ms2_arrays: PartitionedCMGArrayStats<FH>,
    pub id: u64,
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug> Aggregator
    for MultiCMGStatsAgg<FH>
{
    type Context = MsLevelContext<usize, FH>;
    type Item = RawPeak;
    type Output = NaturalFinalizedMultiCMGArrays<FH>;

    fn add(&mut self, peak: impl Into<Self::Item>) {
        let peak = peak.into();
        let u64_intensity = peak.intensity as u64;
        let rt_miliseconds = (peak.retention_time * 1000.0) as u32;

        match self.get_context() {
            MsLevelContext::MS1(_i) => {
                self.ms1_stats.add(
                    rt_miliseconds,
                    peak.scan_index,
                    peak.tof_index,
                    u64_intensity,
                );
            }
            MsLevelContext::MS2(_i) => {
                self.ms2_stats.add(
                    rt_miliseconds,
                    peak.scan_index,
                    peak.tof_index,
                    u64_intensity,
                );
            }
        }
    }

    fn set_context(&mut self, context: Self::Context) {
        match &context {
            MsLevelContext::MS1(i) => {
                self.ms1_stats.set_context(*i).unwrap();
            }
            MsLevelContext::MS2(i) => {
                self.ms2_stats.set_context(i.clone()).unwrap();
            }
        }
        self.context = Some(context);
    }

    fn supports_context(&self) -> bool {
        true
    }

    fn get_context(&self) -> Self::Context {
        match &self.context {
            Some(context) => context.clone(),
            None => panic!("No context set"),
        }
    }

    fn finalize(self) -> NaturalFinalizedMultiCMGArrays<FH> {
        let mz_converter = &self.converters.0;
        let mobility_converter = &self.converters.1;

        let ms1_arrays = PartitionedCMGArrayStats::new(
            PartitionedCMGArrays::from(self.ms1_stats),
            mz_converter,
            mobility_converter,
        );
        let ms2_arrays = PartitionedCMGArrayStats::new(
            PartitionedCMGArrays::from(self.ms2_stats),
            mz_converter,
            mobility_converter,
        );

        NaturalFinalizedMultiCMGArrays {
            ms2_arrays,
            ms1_arrays,
            id: self.id,
        }
    }
}
