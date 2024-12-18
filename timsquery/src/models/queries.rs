use crate::traits::aggregator::{
    NoContext,
    ProvidesContext,
};
use crate::utils::tolerance_ranges::IncludedRange;
use std::collections::HashMap;
use std::hash::Hash;
use timsrust::converters::{
    ConvertableDomain,
    Scan2ImConverter,
    Tof2MzConverter,
};

#[derive(Debug, Clone)]
pub struct PrecursorIndexQuery {
    pub frame_index_range: IncludedRange<usize>,
    pub rt_range_seconds: IncludedRange<f32>,
    pub mz_index_ranges: Vec<IncludedRange<u32>>,
    pub mobility_index_range: IncludedRange<usize>,
    pub isolation_mz_range: IncludedRange<f32>,
}

/// Pretty Generic Definition of the context of a query.
///
/// The context is used to pass additional information to the aggregator.
/// And this implementation simply says  that there will be something passed.
/// if its an MS1 and maybe something else if its an MS2.
///
/// The aggregator should be able to handle this in its definition.
#[derive(Debug, Clone)]
pub enum MsLevelContext<T1, T2> {
    MS1(T1),
    MS2(T2),
}

#[derive(Debug, Clone)]
pub struct FragmentGroupIndexQuery<FH: Clone + Eq + Hash + Send + Sync> {
    pub mz_index_ranges: HashMap<FH, IncludedRange<u32>>,
    pub precursor_query: PrecursorIndexQuery,
}

impl<FH: Clone + Eq + Hash + Send + Sync> ProvidesContext for FragmentGroupIndexQuery<FH> {
    type Context = MsLevelContext<usize, FH>;
}

#[allow(clippy::from_over_into)]
impl Into<NoContext> for MsLevelContext<usize, u64> {
    fn into(self) -> NoContext {
        NoContext {}
    }
}

#[allow(clippy::from_over_into)]
impl Into<NoContext> for MsLevelContext<usize, usize> {
    fn into(self) -> NoContext {
        NoContext {}
    }
}

#[allow(clippy::from_over_into)]
impl Into<NoContext> for MsLevelContext<usize, String> {
    fn into(self) -> NoContext {
        NoContext {}
    }
}

impl<FH: Clone + Eq + Hash + Send + Sync> FragmentGroupIndexQuery<FH> {
    // TODO find if there is a good way to specify in the type that the
    // Only a specific context is returned from each function.

    pub fn iter_ms1_mzs(
        &self,
    ) -> impl Iterator<Item = (MsLevelContext<usize, FH>, IncludedRange<u32>)> + '_ {
        let out = self
            .precursor_query
            .mz_index_ranges
            .iter()
            .enumerate()
            .map(|(i, mz_range)| (MsLevelContext::MS1(i), *mz_range));
        out
    }

    pub fn iter_ms2_mzs(
        &self,
    ) -> impl Iterator<Item = (MsLevelContext<usize, FH>, IncludedRange<u32>)> + '_ {
        let out = self
            .mz_index_ranges
            .iter()
            .map(|(i, mz_range)| (MsLevelContext::MS2(i.clone()), *mz_range));
        out
    }
}

#[derive(Debug, Clone)]
pub struct NaturalPrecursorQuery {
    pub rt_range: IncludedRange<f32>,
    pub mobility_range: IncludedRange<f32>,
    pub mz_ranges: Vec<IncludedRange<f64>>,
    pub isolation_mz_range: IncludedRange<f32>,
}

impl NaturalPrecursorQuery {
    pub fn as_precursor_query(
        &self,
        mz_converter: &Tof2MzConverter,
        im_converter: &Scan2ImConverter,
    ) -> PrecursorIndexQuery {
        PrecursorIndexQuery {
            frame_index_range: (
                im_converter.invert(self.mobility_range.start()).round() as usize,
                im_converter.invert(self.mobility_range.end()).round() as usize,
            )
                .into(),
            rt_range_seconds: self.rt_range,
            mz_index_ranges: self
                .mz_ranges
                .iter()
                .map(|mz_range| {
                    (
                        mz_converter.invert(mz_range.start()).round() as u32,
                        mz_converter.invert(mz_range.end()).round() as u32,
                    )
                        .into()
                })
                .collect(),
            mobility_index_range: (
                im_converter.invert(self.mobility_range.start()).round() as usize,
                im_converter.invert(self.mobility_range.end()).round() as usize,
            )
                .into(),
            isolation_mz_range: self.isolation_mz_range,
        }
    }
}
