use crate::errors::Result;
use crate::models::frames::expanded_frame::{
    ExpandedFrameSlice,
    FrameProcessingConfig,
    SortedState,
    par_read_and_expand_frames,
};
use crate::models::frames::peak_in_quad::PeakInQuad;
use crate::models::frames::single_quad_settings::{
    SingleQuadrupoleSetting,
    SingleQuadrupoleSettingIndex,
    get_matching_quad_settings,
};
use crate::utils::tolerance_ranges::IncludedRange;
use std::collections::HashMap;
use std::time::Instant;
use timsrust::converters::{
    Scan2ImConverter,
    Tof2MzConverter,
};
use timsrust::readers::{
    FrameReader,
    MetadataReader,
};
use tracing::{
    error,
    info,
    instrument,
};
use std::sync::Arc;

#[derive(Debug)]
pub struct ExpandedRawFrameIndex {
    bundled_ms1_frames: ExpandedSliceBundle,
    bundled_frames: HashMap<SingleQuadrupoleSettingIndex, ExpandedSliceBundle>,
    flat_quad_settings: Vec<SingleQuadrupoleSetting>,
    pub cycle_rt_ms: Arc<[u32]>,
    pub mz_converter: Tof2MzConverter,
    pub im_converter: Scan2ImConverter,
}

#[derive(Debug, Clone)]
pub struct ExpandedSliceBundle {
    slices: Vec<ExpandedFrameSlice<SortedState>>,
    frame_indices: Vec<usize>,
}

impl ExpandedSliceBundle {
    pub fn new(mut slices: Vec<ExpandedFrameSlice<SortedState>>) -> Self {
        slices.sort_unstable_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap());
        let frame_indices = slices.iter().map(|x| x.frame_index).collect();
        Self {
            slices,
            frame_indices,
        }
    }

    pub fn query_peaks<F>(
        &self,
        tof_range: IncludedRange<u32>,
        scan_range: Option<IncludedRange<usize>>,
        frame_index_range: IncludedRange<usize>,
        f: &mut F,
    ) where
        F: FnMut(PeakInQuad),
    {
        // Binary search the rt if needed.
        let frame_indices = self.frame_indices.as_slice();
        let low = frame_indices.partition_point(|x| *x < frame_index_range.start());
        let high = frame_indices.partition_point(|x| *x <= frame_index_range.end());

        for i in low..high {
            let slice = &self.slices[i];
            slice.query_peaks(tof_range, scan_range, f);
        }
    }
}

impl ExpandedRawFrameIndex {
    pub fn query_peaks<F>(
        &self,
        tof_range: IncludedRange<u32>,
        precursor_mz_range: IncludedRange<f64>,
        scan_range: Option<IncludedRange<usize>>,
        frame_index_range: IncludedRange<usize>,
        f: &mut F,
    ) where
        F: FnMut(PeakInQuad),
    {
        let matching_quads: Vec<SingleQuadrupoleSettingIndex> =
            get_matching_quad_settings(&self.flat_quad_settings, precursor_mz_range, scan_range)
                .collect();
        self.query_precursor_peaks(&matching_quads, tof_range, scan_range, frame_index_range, f);
    }

    fn query_ms1_peaks<F>(
        &self,
        tof_range: IncludedRange<u32>,
        scan_range: Option<IncludedRange<usize>>,
        frame_index_range: IncludedRange<usize>,
        f: &mut F,
    ) where
        F: FnMut(PeakInQuad),
    {
        self.bundled_ms1_frames
            .query_peaks(tof_range, scan_range, frame_index_range, f);
    }

    fn query_precursor_peaks<F>(
        &self,
        matching_quads: &[SingleQuadrupoleSettingIndex],
        tof_range: IncludedRange<u32>,
        scan_range: Option<IncludedRange<usize>>,
        frame_index_range: IncludedRange<usize>,
        f: &mut F,
    ) where
        F: FnMut(PeakInQuad),
    {
        for quad in matching_quads {
            let tqi = self
                .bundled_frames
                .get(quad)
                .expect("Only existing quads should be queried.");

            tqi.query_peaks(tof_range, scan_range, frame_index_range, &mut *f)
        }
    }

    #[instrument(name = "ExpandedRawFrameIndex::from_path_centroided")]
    pub fn from_path_centroided(path: &str) -> Result<Self> {
        let config = FrameProcessingConfig::default_centroided();
        Self::from_path_base(path, config)
    }

    #[instrument(name = "ExpandedRawFrameIndex::from_path")]
    pub fn from_path(path: &str) -> Result<Self> {
        Self::from_path_base(path, FrameProcessingConfig::NotCentroided)
    }

    #[instrument(name = "ExpandedRawFrameIndex::from_path_base")]
    pub fn from_path_base(path: &str, centroid_config: FrameProcessingConfig) -> Result<Self> {
        info!(
            "Building ExpandedRawFrameIndex from path {} config {:?}",
            path, centroid_config,
        );
        let file_reader = match FrameReader::new(path) {
            Ok(x) => x,
            Err(e) => {
                error!("Failed to open file reader for path {}. Error: {}", path, e);
                return Err(e.into());
            }
        };

        let sql_path = std::path::Path::new(path).join("analysis.tdf");
        let meta_converters = match MetadataReader::new(&sql_path) {
            Ok(x) => x,
            Err(e) => {
                error!(
                    "Failed to open metadata reader for path {}. Error: {}",
                    sql_path.display(),
                    e
                );
                return Err(e.into());
            }
        };
        let centroid_config = centroid_config
            .with_converters(meta_converters.im_converter, meta_converters.mz_converter);

        let st = Instant::now();
        let centroided_split_frames = par_read_and_expand_frames(&file_reader, centroid_config)?;
        let centroided_elap = st.elapsed();
        info!("Reading + Centroiding took {:#?}", centroided_elap);

        let mut out_ms2_frames = HashMap::new();
        let mut out_ms1_frames: Option<ExpandedSliceBundle> = None;

        let mut flat_quad_settings = Vec::new();
        centroided_split_frames
            .into_iter()
            .for_each(|(q, frameslices)| match q {
                None => {
                    out_ms1_frames = Some(ExpandedSliceBundle::new(frameslices));
                }
                Some(q) => {
                    flat_quad_settings.push(q);
                    out_ms2_frames.insert(q.index, ExpandedSliceBundle::new(frameslices));
                }
            });

        let mut cycle_rts_ms: Vec<u32> = out_ms1_frames.as_ref().unwrap().slices.iter().map(|x| (1000.0 * x.rt) as u32).collect();
        cycle_rts_ms.sort_unstable();

        let out = Self {
            bundled_ms1_frames: out_ms1_frames.expect("At least one ms1 frame should be present"),
            bundled_frames: out_ms2_frames,
            flat_quad_settings,
            mz_converter: meta_converters.mz_converter,
            im_converter: meta_converters.im_converter,
            cycle_rt_ms: cycle_rts_ms.into(),
        };

        Ok(out)
    }
}

