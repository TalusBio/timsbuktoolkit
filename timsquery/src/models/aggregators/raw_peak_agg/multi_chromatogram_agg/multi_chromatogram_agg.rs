use super::super::chromatogram_agg::ScanTofStatsCalculatorPair;
use super::aggregator::ParitionedCMGAggregator;
use super::arrays::{
    PartitionedCMGArrayStats,
    PartitionedCMGArrays,
};
use super::scored_arrays::{
    PartitionedCMGScoredStatsArrays,
    ScoresAtTime,
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

        let exp_fragment_intensities: Option<Vec<f32>> =
            if elution_group.expected_fragment_intensity.is_some() {
                let tmp_expect_inten = elution_group.expected_fragment_intensity.as_ref().unwrap();
                let out = fragment_keys
                    .iter()
                    .map(|x| *(tmp_expect_inten.get(x).unwrap_or(&0.0)))
                    .collect();
                Some(out)
            } else {
                None
            };

        MultiCMGStatsAgg {
            converters: (self.converters.0, self.converters.1),
            ms1_stats: ParitionedCMGAggregator::new(
                precursor_keys,
                elution_group.expected_precursor_intensity.to_owned(),
                expect_scan_index,
                precursor_tof,
            ),
            ms2_stats: ParitionedCMGAggregator::new(
                fragment_keys,
                exp_fragment_intensities,
                expect_scan_index,
                fragment_tof,
            ),
            id: elution_group.id,
            context: None,
            buffer: None,
            elution_group_ref,
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct MultiCMGArrayStats<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub ms1_stats: PartitionedCMGArrayStats<usize>,
    pub ms2_stats: PartitionedCMGArrayStats<FH>,
    pub id: u64,
}

#[derive(Debug, Clone, Serialize)]
pub struct NaturalFinalizedMultiCMGStatsArrays<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub ms1_stats: PartitionedCMGScoredStatsArrays<usize>,
    pub ms2_stats: PartitionedCMGScoredStatsArrays<FH>,
    pub cross_ms1: Vec<f64>,
    pub cross_ms2: Vec<f64>,
    pub id: u64,
}

#[derive(Debug, Clone, Serialize)]
pub struct ApexScores {
    pub ms1_scores: ScoresAtTime,
    pub ms2_scores: ScoresAtTime,
    pub main_score: f64,
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug>
    NaturalFinalizedMultiCMGStatsArrays<FH>
{
    pub fn finalized_score(&self) -> Result<ApexScores> {
        let main_score = self.cross_ms2.as_slice();
        let main_score_rts = self.ms1_stats.retention_time_miliseconds.as_slice();
        let apex_cross_ms2_index = main_score.iter().enumerate().max_by(|x, y| {
            let cmp = x.1.partial_cmp(y.1);
            if cmp.is_none() {
                panic!(
                    "Main score is not comparable {:?} vs {:?} ... [{:?}]",
                    x, y, main_score
                );
            }
            cmp.unwrap()
        });

        let apex_cross_ms2_index = match apex_cross_ms2_index {
            Some(x) => x,
            None => {
                return Err(TimsqueryError::Other(format!(
                    "Insufficient data for cross score, main_score: {:?}",
                    main_score
                )));
            }
        };

        let apex_rt = main_score_rts[apex_cross_ms2_index.0];

        let ms1_scores = self.ms1_stats.scores_at_rt(apex_rt);
        let ms2_scores = self.ms2_stats.scores_at_rt(apex_rt);

        let out = ApexScores {
            ms1_scores,
            ms2_scores,
            main_score: *apex_cross_ms2_index.1,
        };

        if out.main_score.is_infinite() {
            panic!("Main score is infinite, main_score: {:?}", main_score);
        }
        Ok(out)
    }
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug> Aggregator
    for MultiCMGStatsAgg<FH>
{
    type Context = MsLevelContext<usize, FH>;
    type Item = RawPeak;
    type Output = NaturalFinalizedMultiCMGStatsArrays<FH>;

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

    fn finalize(self) -> NaturalFinalizedMultiCMGStatsArrays<FH> {
        let mz_converter = &self.converters.0;
        let mobility_converter = &self.converters.1;

        // TODO clean this up ... 4 conversions in one line is feeling a bit weird.
        let ms1_stats = PartitionedCMGScoredStatsArrays::new(
            PartitionedCMGArrayStats::from(PartitionedCMGArrays::from(self.ms1_stats)),
            mz_converter,
            mobility_converter,
        )
        .unwrap();
        let ms2_stats = PartitionedCMGScoredStatsArrays::new(
            PartitionedCMGArrayStats::from(PartitionedCMGArrays::from(self.ms2_stats)),
            mz_converter,
            mobility_converter,
        )
        .unwrap();

        let cross_ms1 = ms2_stats.cross_scores(&ms1_stats).unwrap_or(vec![]);
        let cross_ms2 = ms1_stats.cross_scores(&ms2_stats).unwrap_or(vec![]);

        NaturalFinalizedMultiCMGStatsArrays {
            ms2_stats,
            ms1_stats,
            cross_ms1,
            cross_ms2,
            id: self.id,
        }
    }
}
