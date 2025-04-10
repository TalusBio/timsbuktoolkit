use super::quad_index::{
    TransposedQuadIndex,
    TransposedQuadIndexBuilder,
};
use crate::ToleranceAdapter;
use crate::errors::Result;
use crate::models::adapters::FragmentIndexAdapter;
use crate::models::elution_group::ElutionGroup;
use crate::models::frames::expanded_frame::{
    ExpandedFrameSlice,
    FrameProcessingConfig,
    SortingStateTrait,
    par_read_and_expand_frames,
};
use crate::models::frames::peak_in_quad::PeakInQuad;
use crate::models::frames::raw_peak::RawPeak;
use crate::models::frames::single_quad_settings::{
    SingleQuadrupoleSetting,
    SingleQuadrupoleSettingIndex,
    get_matching_quad_settings,
    matches_quad_settings,
};
use crate::models::queries::{
    FragmentGroupIndexQuery,
    MsLevelContext,
};
use crate::traits::aggregator::Aggregator;
use crate::traits::key_like::KeyLike;
use crate::traits::queriable_data::QueriableData;
use crate::utils::display::{
    GlimpseConfig,
    glimpse_vec,
};
use crate::utils::tolerance_ranges::IncludedRange;
use rayon::prelude::*;
use serde::Serialize;
use std::collections::HashMap;
use std::fmt::{
    Debug,
    Display,
};
use std::hash::Hash;
use std::time::Instant;
use timsrust::Metadata;
use timsrust::converters::{
    Scan2ImConverter,
    Tof2MzConverter,
};
use timsrust::readers::{
    FrameReader,
    MetadataReader,
};
use tracing::{
    debug,
    info,
    instrument,
};

// TODO break this module apart ... its getting too big for my taste
// - JSP: 2024-11-19

pub struct QuadSplittedTransposedIndex {
    precursor_index: TransposedQuadIndex,
    fragment_indices: HashMap<SingleQuadrupoleSettingIndex, TransposedQuadIndex>,
    flat_quad_settings: Vec<SingleQuadrupoleSetting>,
    pub mz_converter: Tof2MzConverter,
    pub im_converter: Scan2ImConverter,
    adapter: FragmentIndexAdapter,
}

impl Debug for QuadSplittedTransposedIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "QuadSplittedTransposedIndex(num_quads: {})",
            self.flat_quad_settings.len()
        )
    }
}

impl QuadSplittedTransposedIndex {
    pub fn query_ms2_peaks<F>(
        &self,
        tof_range: IncludedRange<u32>,
        precursor_mz_range: IncludedRange<f64>,
        scan_range: Option<IncludedRange<usize>>,
        rt_range_seconds: Option<IncludedRange<f32>>,
        f: &mut F,
    ) where
        F: FnMut(PeakInQuad),
    {
        let matching_quads: Vec<SingleQuadrupoleSettingIndex> = self
            .get_matching_quad_settings(precursor_mz_range, scan_range)
            .collect();
        self.query_peaks_in_precursors(&matching_quads, tof_range, scan_range, rt_range_seconds, f);
    }

    fn query_ms1_peaks<F>(
        &self,
        tof_range: IncludedRange<u32>,
        scan_range: Option<IncludedRange<usize>>,
        rt_range_seconds: Option<IncludedRange<f32>>,
        f: &mut F,
    ) where
        F: FnMut(PeakInQuad),
    {
        self.precursor_index
            .query_peaks(tof_range, scan_range, rt_range_seconds)
            .for_each(f);
    }

    fn query_peaks_in_precursors<F>(
        &self,
        matching_quads: &[SingleQuadrupoleSettingIndex],
        tof_range: IncludedRange<u32>,
        scan_range: Option<IncludedRange<usize>>,
        rt_range_seconds: Option<IncludedRange<f32>>,
        f: &mut F,
    ) where
        F: FnMut(PeakInQuad),
    {
        for quad in matching_quads {
            let tqi = self
                .fragment_indices
                .get(quad)
                .expect("Only existing quads should be queried.");
            tqi.query_peaks(tof_range, scan_range, rt_range_seconds)
                .for_each(&mut *f);
        }
    }

    fn get_matching_quad_settings(
        &self,
        precursor_mz_range: IncludedRange<f64>,
        scan_range: Option<IncludedRange<usize>>,
    ) -> impl Iterator<Item = SingleQuadrupoleSettingIndex> + '_ {
        get_matching_quad_settings(&self.flat_quad_settings, precursor_mz_range, scan_range)
    }
}

impl Display for QuadSplittedTransposedIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut disp_str = String::new();
        disp_str.push_str("QuadSplittedTransposedIndex\n");

        disp_str.push_str(&format!("mz_converter: {:?}\n", self.mz_converter));
        disp_str.push_str(&format!("im_converter: {:?}\n", self.im_converter));
        disp_str.push_str("flat_quad_settings: \n");
        disp_str.push_str(&glimpse_vec(
            &self.flat_quad_settings,
            Some(GlimpseConfig {
                max_items: 10,
                padding: 4,
                new_line: true,
            }),
        ));

        disp_str.push_str("precursor_index: \n");
        disp_str.push_str(&format!(" -- {}\n", self.precursor_index));
        let mut num_shown = 0;
        for (qs, tqi) in self.fragment_indices.iter() {
            disp_str.push_str(&format!(" - {:?}: \n", qs));
            disp_str.push_str(&format!(" -- {}\n", tqi));
            num_shown += 1;
            if num_shown > 5 {
                disp_str.push_str(&format!(
                    " ........ len = {}\n",
                    self.fragment_indices.len()
                ));
                break;
            }
        }

        write!(f, "{}", disp_str)
    }
}

impl QuadSplittedTransposedIndex {
    #[instrument(name = "QuadSplittedTransposedIndex::from_path", level = "debug")]
    pub fn from_path(path: &str) -> Result<Self> {
        let st = Instant::now();
        info!("Building transposed quad index from path {}", path);
        let tmp = QuadSplittedTransposedIndexBuilder::from_path(path)?;
        let out = tmp.build();
        let elapsed = st.elapsed();
        info!("Transposed quad index built in {:#?}", elapsed);
        debug!("{}", out);
        Ok(out)
    }

    #[instrument(
        name = "QuadSplittedTransposedIndex::from_path_centroided",
        level = "debug"
    )]
    pub fn from_path_centroided(path: &str) -> Result<Self> {
        let st = Instant::now();
        info!(
            "Building CENTROIDED transposed quad index from path {}",
            path
        );
        let tmp = QuadSplittedTransposedIndexBuilder::from_path_centroided(path)?;
        let out = tmp.build();
        let elapsed = st.elapsed();
        info!("Transposed CENTROIDED quad index built in {:#?}", elapsed);
        debug!("{}", out);
        Ok(out)
    }
}

#[derive(Debug, Clone, Default)]
pub struct QuadSplittedTransposedIndexBuilder {
    indices: HashMap<Option<SingleQuadrupoleSetting>, TransposedQuadIndexBuilder>,
    mz_converter: Option<Tof2MzConverter>,
    im_converter: Option<Scan2ImConverter>,
    metadata: Option<Metadata>,
    // TODO use during build to make sure we
    // have a the right number of peaks in the end.
    // ... this means I need to implement len for TransposedQuadIndexBuilder
    added_peaks: u64,
}

impl QuadSplittedTransposedIndexBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    fn add_frame_slice<S: SortingStateTrait>(&mut self, frame_slice: ExpandedFrameSlice<S>) {
        // Add key if it doesnt exist ...
        self.indices
            .entry(frame_slice.quadrupole_settings)
            .or_insert(TransposedQuadIndexBuilder::new(
                frame_slice.quadrupole_settings,
            ));

        self.added_peaks += frame_slice.len() as u64;
        self.indices
            .get_mut(&frame_slice.quadrupole_settings)
            .unwrap()
            .add_frame_slice(frame_slice);
    }

    #[instrument(
        name = "QuadSplittedTransposedIndexBuilder::from_path",
        level = "debug"
    )]
    fn from_path(path: &str) -> Result<Self> {
        Self::from_path_base(path, FrameProcessingConfig::NotCentroided)
    }

    #[instrument(
        name = "QuadSplittedTransposedIndexBuilder::from_path_centroided",
        level = "debug"
    )]
    fn from_path_centroided(path: &str) -> Result<Self> {
        let config = FrameProcessingConfig::default_centroided();
        Self::from_path_base(path, config)
    }

    // TODO: I think i should split this into two functions, one that starts the builder
    // and one that adds the frameslices, maybe even have a config struct that dispatches
    // the right preprocessing steps.
    #[instrument(
        name = "QuadSplittedTransposedIndexBuilder::from_path_base",
        level = "debug"
    )]
    fn from_path_base(path: &str, centroid_config: FrameProcessingConfig) -> Result<Self> {
        let file_reader = FrameReader::new(path)?;

        let sql_path = std::path::Path::new(path).join("analysis.tdf");
        let meta_converters = MetadataReader::new(&sql_path)?;

        let out_meta_converters = meta_converters.clone();
        let mut final_out = Self {
            indices: HashMap::new(),
            mz_converter: Some(meta_converters.mz_converter),
            im_converter: Some(meta_converters.im_converter),
            metadata: Some(out_meta_converters),
            added_peaks: 0,
        };

        let centroid_config = centroid_config
            .with_converters(meta_converters.im_converter, meta_converters.mz_converter);

        let split_frames = par_read_and_expand_frames(&file_reader, centroid_config)?;

        // TODO use the rayon contructor to fold
        let out2: Result<Vec<Self>> = split_frames
            .into_par_iter()
            .map(|(_q, frameslices)| {
                // TODO:Refactor so the internal index is built first and then added.
                // This should save a couple of thousand un-necessary hashmap lookups.
                let mut out = Self::new();
                for frameslice in frameslices {
                    out.add_frame_slice(frameslice);
                }
                Ok(out)
            })
            .collect();

        let out2 = out2?.into_iter().fold(Self::new(), |mut x, y| {
            x.fold(y);
            x
        });
        final_out.fold(out2);

        Ok(final_out)
    }

    pub fn fold(&mut self, other: Self) {
        for (qs, builder) in other.indices.into_iter() {
            self.indices
                .entry(qs)
                .and_modify(|bl| bl.fold(builder.clone()))
                .or_insert(builder);
        }
        self.added_peaks += other.added_peaks;
    }

    #[instrument(skip(self), level = "debug")]
    pub fn build(self) -> QuadSplittedTransposedIndex {
        let mut indices = HashMap::new();
        let mut flat_quad_settings = Vec::new();
        let built: Vec<(TransposedQuadIndex, Option<SingleQuadrupoleSetting>)> = self
            .indices
            .into_par_iter()
            .map(|(qs, builder)| (builder.build(), qs))
            .collect();

        let mut precursor_index: Option<TransposedQuadIndex> = None;
        for (qi, qs) in built.into_iter() {
            if let Some(qs) = qs {
                indices.insert(qs.index, qi);
                flat_quad_settings.push(qs);
            } else {
                precursor_index = Some(qi);
            }
        }

        flat_quad_settings.sort_by(|a, b| {
            a.ranges
                .isolation_mz
                .partial_cmp(&b.ranges.isolation_mz)
                .unwrap()
        });

        QuadSplittedTransposedIndex {
            precursor_index: precursor_index.expect("Precursor peaks should be present"),
            fragment_indices: indices,
            flat_quad_settings,
            mz_converter: self.mz_converter.unwrap(),
            im_converter: self.im_converter.unwrap(),
            adapter: FragmentIndexAdapter::from(self.metadata.unwrap()),
        }
    }
}

impl<FH: KeyLike> QueriableData<FragmentGroupIndexQuery<FH>, RawPeak, MsLevelContext<usize, FH>>
    for QuadSplittedTransposedIndex
{
    fn query(&self, fragment_query: &FragmentGroupIndexQuery<FH>) -> Vec<RawPeak> {
        let precursor_mz_range = IncludedRange::new(
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
        );
        let scan_range = Some(fragment_query.precursor_query.mobility_index_range);

        fragment_query
            .iter_ms2_mzs()
            .flat_map(|(_fh, tof_range)| {
                let mut local_vec: Vec<RawPeak> = vec![];
                self.query_ms2_peaks(
                    tof_range,
                    precursor_mz_range,
                    scan_range,
                    Some(fragment_query.precursor_query.rt_range_seconds),
                    &mut |x| local_vec.push(RawPeak::from(x)),
                );

                local_vec
            })
            .collect()
    }

    fn add_query<A, O, AG, C2>(
        &self,
        fragment_query: &FragmentGroupIndexQuery<FH>,
        aggregator: &mut AG,
    ) where
        A: From<RawPeak> + Send + Sync + Clone + Copy,
        AG: Aggregator<Item = A, Output = O, Context = C2>,
        MsLevelContext<usize, FH>: Into<C2>,
    {
        let precursor_mz_range = IncludedRange::new(
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
        );
        let scan_range = Some(fragment_query.precursor_query.mobility_index_range);

        fragment_query.iter_ms2_mzs().for_each(|(fh, tof_range)| {
            if aggregator.supports_context() {
                aggregator.set_context(fh.into());
            }

            self.query_ms2_peaks(
                tof_range,
                precursor_mz_range,
                scan_range,
                Some(fragment_query.precursor_query.rt_range_seconds),
                &mut |peak| aggregator.add(RawPeak::from(peak)),
            );
        })
    }

    fn add_query_multi_group<A, O, AG, C2>(
        &self,
        fragment_queries: &[FragmentGroupIndexQuery<FH>],
        aggregator: &mut [AG],
    ) where
        A: From<RawPeak> + Send + Sync + Clone + Copy,
        AG: Aggregator<Item = A, Output = O, Context = C2>,
        MsLevelContext<usize, FH>: Into<C2>,
    {
        let prec_mz_ranges = fragment_queries
            .iter()
            .map(|x| {
                (
                    x.precursor_query.isolation_mz_range.start() as f64,
                    x.precursor_query.isolation_mz_range.end() as f64,
                )
                    .into()
            })
            .collect::<Vec<_>>();

        let scan_ranges = fragment_queries
            .iter()
            .map(|x| Some(x.precursor_query.mobility_index_range))
            .collect::<Vec<_>>();

        let rt_ranges = fragment_queries
            .iter()
            .map(|x| x.precursor_query.rt_range_seconds)
            .collect::<Vec<_>>();

        // Query the ms1 mzs first.
        aggregator.par_iter_mut().enumerate().for_each(|(i, agg)| {
            fragment_queries[i]
                .iter_ms1_mzs()
                .for_each(|(fh, mz_range)| {
                    if agg.supports_context() {
                        agg.set_context(fh.into());
                    }
                    self.query_ms1_peaks(mz_range, scan_ranges[i], Some(rt_ranges[i]), &mut |x| {
                        agg.add(RawPeak::from(x))
                    });
                });
        });

        for quad_setting in self.flat_quad_settings.iter() {
            let local_index = quad_setting.index;

            let tqi = self
                .fragment_indices
                .get(&local_index)
                .expect("Only existing quads should be queried.");

            aggregator.par_iter_mut().enumerate().for_each(|(i, agg)| {
                if !matches_quad_settings(quad_setting, prec_mz_ranges[i], scan_ranges[i]) {
                    return;
                }

                for (fh, tof_range) in fragment_queries[i].iter_ms2_mzs() {
                    if agg.supports_context() {
                        agg.set_context(fh.into());
                    }
                    tqi.query_peaks(tof_range, scan_ranges[i], Some(rt_ranges[i]))
                        .for_each(|x| agg.add(RawPeak::from(x)));
                }
            });
        }
    }
}

impl<FH: Clone + Serialize + Eq + Hash + Send + Sync + std::fmt::Debug>
    ToleranceAdapter<FragmentGroupIndexQuery<FH>, ElutionGroup<FH>>
    for QuadSplittedTransposedIndex
{
    fn query_from_elution_group(
        &self,
        tol: &dyn crate::traits::tolerance::Tolerance,
        elution_group: &ElutionGroup<FH>,
    ) -> FragmentGroupIndexQuery<FH> {
        self.adapter.query_from_elution_group(tol, elution_group)
    }
}
