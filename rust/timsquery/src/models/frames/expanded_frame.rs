use super::peak_in_quad::PeakInQuad;
use super::single_quad_settings::{
    ExpandedFrameQuadSettings,
    SingleQuadrupoleSetting,
    expand_quad_settings,
};
use crate::errors::{
    Result,
    UnsupportedDataError,
};
use crate::sort_vecs_by_first;
use crate::utils::compress_explode::explode_vec;
use crate::utils::frame_processing::{
    PeakArrayRefs,
    lazy_centroid_weighted_frame,
};
use crate::utils::sorting::top_n;
use crate::utils::tolerance_ranges::{
    IncludedRange,
    scan_tol_range,
    tof_tol_range,
};
use rayon::prelude::*;
use std::borrow::Cow;
use std::collections::HashMap;
use std::marker::PhantomData;
use std::sync::Arc;
use timsrust::converters::{
    Scan2ImConverter,
    Tof2MzConverter,
};
use timsrust::readers::{
    FrameReader,
    FrameReaderError,
};
use timsrust::{
    AcquisitionType,
    Frame,
    MSLevel,
    QuadrupoleSettings,
};
use tracing::{
    info,
    instrument,
    trace,
    warn,
};

/// A frame after expanding the mobility data and re-sorting it by tof.
#[derive(Debug, Clone)]
pub struct ExpandedFrame {
    pub tof_indices: Vec<u32>,
    pub scan_numbers: Vec<usize>,
    pub intensities: Vec<u32>,
    pub frame_index: usize,
    pub rt: f64,
    pub acquisition_type: AcquisitionType,
    pub ms_level: MSLevel,
    pub quadrupole_settings: Arc<QuadrupoleSettings>,
    pub intensity_correction_factor: f64,
    pub window_group: u8,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum UnsortedState {}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SortedState {}

pub trait SortingStateTrait {}
impl SortingStateTrait for UnsortedState {}
impl SortingStateTrait for SortedState {}

#[derive(Debug, Clone)]
pub struct ExpandedFrameSlice<S: SortingStateTrait> {
    pub tof_indices: Vec<u32>, // I could Arc<[u32]> if I didnt want to sort by it ...
    pub scan_numbers: Vec<usize>,
    pub intensities: Vec<u32>,
    pub frame_index: usize,
    pub rt: f64,
    pub acquisition_type: AcquisitionType,
    pub ms_level: MSLevel,
    pub quadrupole_settings: Option<SingleQuadrupoleSetting>,
    pub intensity_correction_factor: f64,
    pub window_group: u8,
    pub window_subindex: u8,
    _sorting_state: PhantomData<S>,
}

// Example of how to use it with your specific types
fn sort_by_tof_macro(
    tof_indices: Vec<u32>,
    scan_numbers: Vec<usize>,
    intensities: Vec<u32>,
) -> (Vec<u32>, Vec<usize>, Vec<u32>) {
    sort_vecs_by_first!(tof_indices, scan_numbers, intensities)
}

impl<T: SortingStateTrait> ExpandedFrameSlice<T> {
    pub fn sort_by_tof(self) -> ExpandedFrameSlice<SortedState> {
        // let mut indices = argsort_by(&self.tof_indices, |x| *x);
        // sort_by_indices_multi!(
        //     &mut indices,
        //     &mut self.tof_indices,
        //     &mut self.scan_numbers,
        //     &mut self.intensities
        // );

        let (tof_indices, scan_numbers, intensities) =
            sort_by_tof_macro(self.tof_indices, self.scan_numbers, self.intensities);

        ExpandedFrameSlice {
            tof_indices,
            scan_numbers,
            intensities,
            frame_index: self.frame_index,
            rt: self.rt,
            acquisition_type: self.acquisition_type,
            ms_level: self.ms_level,
            quadrupole_settings: self.quadrupole_settings,
            intensity_correction_factor: self.intensity_correction_factor,
            window_group: self.window_group,
            window_subindex: self.window_subindex,
            _sorting_state: PhantomData,
        }
    }

    pub fn len(&self) -> usize {
        self.tof_indices.len()
    }

    pub fn is_empty(&self) -> bool {
        self.tof_indices.is_empty()
    }
}

impl ExpandedFrameSlice<SortedState> {
    pub fn query_peaks<F>(
        &self,
        tof_range: IncludedRange<u32>,
        scan_range: Option<IncludedRange<usize>>,
        f: &mut F,
    ) where
        F: FnMut(PeakInQuad),
    {
        // TODO abstract this operation... Range -> Range
        let peak_range = {
            let peak_ind_start = self.tof_indices.partition_point(|x| *x < tof_range.start());
            let peak_ind_end = self.tof_indices.partition_point(|x| *x <= tof_range.end());
            peak_ind_start..peak_ind_end
        };

        for peak_ind in peak_range {
            let scan_index = self.scan_numbers[peak_ind];
            if let Some(scan_range) = &scan_range {
                if !scan_range.contains(scan_index) {
                    continue;
                }
            }

            let intensity = self.intensities[peak_ind];
            let tof_index = self.tof_indices[peak_ind];
            let retention_time = self.rt;
            f(PeakInQuad {
                scan_index,
                intensity,
                tof_index,
                retention_time: retention_time as f32,
            });
        }
    }
}

fn trim_scan_edges(scan_start: usize, scan_end: usize) -> (usize, usize) {
    // Poor man's fix to different quad windows bleeding onto each other ...
    // I will trim the smallest of 20 scans or 10% of the range size.
    //
    // old dumber implementation
    // let peak_start = frame.scan_offsets[qs.ranges.scan_start + 10];
    // let peak_end = frame.scan_offsets[qs.ranges.scan_end - 10];
    // (peak_start, peak_end)
    //
    let scan_range = scan_end - scan_start;
    assert!(scan_range > 0, "Expected scan range to be positive...");
    let offset_use = (scan_range / 10).min(20);

    let new_start = scan_start + offset_use;
    let new_end = scan_end - offset_use;

    assert!(new_start < new_end, "Expected new_start < new_end");

    (new_start, new_end)
}

fn expand_unfragmented_frame(frame: Frame) -> ExpandedFrameSlice<SortedState> {
    let scan_numbers = explode_vec(&frame.scan_offsets);
    let intensities = frame.intensities;
    let tof_indices = frame.tof_indices;
    let curr_slice = ExpandedFrameSlice {
        tof_indices,
        scan_numbers,
        intensities,
        frame_index: frame.index,
        rt: frame.rt_in_seconds,
        acquisition_type: frame.acquisition_type,
        ms_level: frame.ms_level,
        quadrupole_settings: None,
        intensity_correction_factor: frame.intensity_correction_factor,
        window_group: frame.window_group,
        window_subindex: 0u8,
        _sorting_state: PhantomData::<UnsortedState>,
    };

    curr_slice.sort_by_tof()
}

fn expand_fragmented_frame(
    frame: Frame,
    quads: Vec<SingleQuadrupoleSetting>,
) -> Vec<ExpandedFrameSlice<SortedState>> {
    let mut out = Vec::with_capacity(quads.len());
    let exploded_scans = explode_vec(&frame.scan_offsets);
    for (i, qs) in quads.into_iter().enumerate() {
        // This whole block is kind of ugly but I think its good enough for now ...
        let (slice_scan_start, slice_scan_end) =
            trim_scan_edges(qs.ranges.scan_start, qs.ranges.scan_end);

        let tof_index_index_slice_start = frame.scan_offsets[slice_scan_start];
        let tof_index_index_slice_end = frame.scan_offsets[slice_scan_end];

        let slice_tof_indices =
            frame.tof_indices[tof_index_index_slice_start..tof_index_index_slice_end].to_vec();
        let slice_scan_numbers =
            exploded_scans[tof_index_index_slice_start..tof_index_index_slice_end].to_vec();
        let slice_intensities =
            frame.intensities[tof_index_index_slice_start..tof_index_index_slice_end].to_vec();

        let curr_slice = ExpandedFrameSlice {
            tof_indices: slice_tof_indices,
            scan_numbers: slice_scan_numbers,
            intensities: slice_intensities,
            frame_index: frame.index,
            rt: frame.rt_in_seconds,
            acquisition_type: frame.acquisition_type,
            ms_level: frame.ms_level,
            quadrupole_settings: Some(qs),
            intensity_correction_factor: frame.intensity_correction_factor,
            window_group: frame.window_group,
            window_subindex: i as u8,
            _sorting_state: PhantomData::<UnsortedState>,
        };
        out.push(curr_slice.sort_by_tof());
    }
    out
}

pub fn expand_and_split_frame(frame: Frame) -> Vec<ExpandedFrameSlice<SortedState>> {
    let quad_settings = &frame.quadrupole_settings;
    // Q: Can I save a vec allocation if I make the expanded quad_settings
    // To return a vec of builders? or Enum(settings|frame_slice)?
    let expanded_quad_settings = expand_quad_settings(quad_settings);

    match expanded_quad_settings {
        ExpandedFrameQuadSettings::Unfragmented => vec![expand_unfragmented_frame(frame)],
        ExpandedFrameQuadSettings::Fragmented(quads) => expand_fragmented_frame(frame, quads),
    }
}

impl ExpandedFrame {
    pub fn from_frame(frame: Frame) -> Self {
        let scan_numbers = explode_vec(&frame.scan_offsets);

        let sorted = sort_vecs_by_first!(frame.tof_indices, scan_numbers, frame.intensities);
        let tof_indices = sorted.0;
        let scan_numbers = sorted.1;
        let intensities = sorted.2;

        ExpandedFrame {
            tof_indices,
            scan_numbers,
            intensities,
            frame_index: frame.index,
            rt: frame.rt_in_seconds,
            acquisition_type: frame.acquisition_type,
            ms_level: frame.ms_level,
            quadrupole_settings: frame.quadrupole_settings,
            intensity_correction_factor: frame.intensity_correction_factor,
            window_group: frame.window_group,
        }
    }
}

type QuadBundledFrameslices =
    HashMap<Option<SingleQuadrupoleSetting>, Vec<ExpandedFrameSlice<SortedState>>>;

#[instrument(skip(frame_iter))]
pub fn par_expand_and_arrange_frames(
    frame_iter: impl ParallelIterator<Item = Frame>,
) -> QuadBundledFrameslices {
    let split = frame_iter.flat_map(expand_and_split_frame);

    // Implementation so the folding occurs in parrallel.
    // Inspired by/Copied from: https://stackoverflow.com/a/70097253/4295016
    let mut out: QuadBundledFrameslices = split
        .fold(
            HashMap::new,
            |mut acc: QuadBundledFrameslices, x: ExpandedFrameSlice<SortedState>| {
                acc.entry(x.quadrupole_settings).or_default().push(x);
                acc
            },
        )
        .reduce_with(|mut acc1, acc2| {
            for (qs, frameslices) in acc2.into_iter() {
                acc1.entry(qs).or_default().extend(frameslices);
            }
            acc1
        })
        .expect("At least one frame is iterated over");

    // Finally sort each of them internally by retention time.
    for (_, es) in out.iter_mut() {
        es.par_sort_unstable_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap());
    }
    out
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CentroidingSettings {
    pub ims_tol_pct: f64,
    pub mz_tol_ppm: f64,
    pub window_width: usize,
}

impl Default for CentroidingSettings {
    fn default() -> Self {
        CentroidingSettings {
            ims_tol_pct: 1.5,
            mz_tol_ppm: 15.0,
            window_width: 3,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FrameProcessingConfig {
    Centroided {
        settings: CentroidingSettings,
        ims_converter: Option<Scan2ImConverter>,
        mz_converter: Option<Tof2MzConverter>,
    },
    NotCentroided,
}

impl FrameProcessingConfig {
    pub fn with_converters(
        self,
        ims_converter: Scan2ImConverter,
        mz_converter: Tof2MzConverter,
    ) -> Self {
        match self {
            FrameProcessingConfig::Centroided { settings, .. } => {
                FrameProcessingConfig::Centroided {
                    settings,
                    ims_converter: Some(ims_converter),
                    mz_converter: Some(mz_converter),
                }
            }
            FrameProcessingConfig::NotCentroided => FrameProcessingConfig::NotCentroided,
        }
    }

    pub fn default_centroided() -> Self {
        FrameProcessingConfig::Centroided {
            settings: Default::default(),
            ims_converter: Default::default(),
            mz_converter: Default::default(),
        }
    }

    pub fn default_not_centroided() -> Self {
        FrameProcessingConfig::NotCentroided
    }
}

pub fn warn_and_skip_badframes(
    frame_iter: impl ParallelIterator<Item = std::result::Result<Frame, FrameReaderError>>,
) -> impl ParallelIterator<Item = Frame> {
    frame_iter.filter_map(|x| {
        // Log the info of the frame that broke ...
        match x {
            Ok(frame) => Some(frame),
            Err(e) => {
                warn!("Failed to read frame {:?}", e);
                None
            }
        }
    })
}

#[instrument(
    skip(frame_reader),
    fields(
        num_frames = %frame_reader.len(),
        // get_path was removed in timsrust 0.4.2
        // path = frame_reader.get_path().to_str(),
    )
)]
pub fn par_read_and_expand_frames(
    frame_reader: &FrameReader,
    centroiding_config: FrameProcessingConfig,
) -> Result<HashMap<Option<SingleQuadrupoleSetting>, Vec<ExpandedFrameSlice<SortedState>>>> {
    let dia_windows = match frame_reader.get_dia_windows() {
        Some(dia_windows) => dia_windows,
        None => {
            return Err(UnsupportedDataError::NoMS2DataError.into());
        }
    };

    let mut all_expanded_frames = HashMap::new();
    for dia_window in dia_windows.into_iter() {
        info!("Processing dia window: {:?}", dia_window);
        let curr_iter = frame_reader.parallel_filter(|x| x.quadrupole_settings == dia_window);
        let curr_iter = warn_and_skip_badframes(curr_iter);

        let expanded_frames = match centroiding_config {
            FrameProcessingConfig::Centroided {
                settings,
                ims_converter,
                mz_converter,
            } => par_expand_and_centroid_frames(
                curr_iter,
                settings.ims_tol_pct,
                settings.mz_tol_ppm,
                settings.window_width,
                &ims_converter.unwrap(),
                &mz_converter.unwrap(),
            ),
            FrameProcessingConfig::NotCentroided => par_expand_and_arrange_frames(curr_iter),
        };

        all_expanded_frames.extend(expanded_frames);
    }

    let slice_infos: Vec<ExpandedQuadSliceInfo> = all_expanded_frames
        .values()
        .map(|x| ExpandedQuadSliceInfo::new(x))
        .collect();

    let tbl = Table::new(slice_infos);
    println!("Slice info: \n{}", tbl);

    info!("Processing MS1 frames");
    let ms1_iter = frame_reader.parallel_filter(|x| x.ms_level == MSLevel::MS1);
    let ms1_iter = warn_and_skip_badframes(ms1_iter);
    let expanded_ms1_frames = match centroiding_config {
        FrameProcessingConfig::Centroided {
            settings,
            ims_converter,
            mz_converter,
        } => par_expand_and_centroid_frames(
            ms1_iter,
            settings.ims_tol_pct,
            settings.mz_tol_ppm,
            settings.window_width,
            &ims_converter.unwrap(),
            &mz_converter.unwrap(),
        ),
        FrameProcessingConfig::NotCentroided => par_expand_and_arrange_frames(ms1_iter),
    };
    all_expanded_frames.extend(expanded_ms1_frames);
    info!("Done reading and expanding frames");

    Ok(all_expanded_frames)
}

use tabled::{
    Table,
    Tabled,
};

#[derive(Debug, Clone, Copy)]
pub struct ExpandedQuadSliceInfo {
    pub quad_settings: Option<SingleQuadrupoleSetting>,
    pub cycle_time_seconds: f64,
    pub peak_width_seconds: Option<f64>,
}

impl Tabled for ExpandedQuadSliceInfo {
    const LENGTH: usize = 7;

    // Required methods
    fn fields(&self) -> Vec<Cow<'_, str>> {
        match self.quad_settings {
            Some(qs) => {
                let qr = qs.ranges;
                let qi = qs.index;

                let scan_ranges = format!("{:?}-{:?}", qr.scan_start, qr.scan_end);
                let quad_ranges = format!("{:?}-{:?}", qr.isolation_low, qr.isolation_high);
                vec![
                    Cow::Owned(qi.major_index.to_string()),
                    Cow::Owned(qi.sub_index.to_string()),
                    Cow::Owned(quad_ranges),
                    Cow::Owned(scan_ranges),
                    Cow::Owned(qr.collision_energy.to_string()),
                    Cow::Owned(self.cycle_time_seconds.to_string()),
                    Cow::Owned(
                        self.peak_width_seconds
                            .as_ref()
                            .map(|x| x.to_string())
                            .unwrap_or("None".to_string()),
                    ),
                ]
            }
            None => vec![
                Cow::Borrowed("None"),
                Cow::Borrowed("None"),
                Cow::Owned("None".to_string()),
                Cow::Owned("None".to_string()),
                Cow::Owned("None".to_string()),
                Cow::Owned("None".to_string()),
            ],
        }
    }

    fn headers() -> Vec<Cow<'static, str>> {
        vec![
            Cow::Borrowed("Major index"),
            Cow::Borrowed("Minor index"),
            Cow::Borrowed("Quad ranges"),
            Cow::Borrowed("Scan ranges"),
            Cow::Borrowed("CE"),
            Cow::Borrowed("Cycle time"),
            Cow::Borrowed("Peak width"),
        ]
    }
}

impl ExpandedQuadSliceInfo {
    #[instrument(skip(frameslices), ret, level = "debug")]
    pub fn new(frameslices: &[ExpandedFrameSlice<SortedState>]) -> Self {
        let avg_cycle_time = frameslices
            .windows(2)
            .map(|x| {
                let a = &x[0];
                let b = &x[1];
                let diff = b.rt - a.rt;
                assert!(diff > 0.0);
                diff
            })
            .sum::<f64>()
            / frameslices.len() as f64;

        let peak_width = Self::estimate_peak_width(frameslices);

        Self {
            quad_settings: frameslices[0].quadrupole_settings,
            cycle_time_seconds: avg_cycle_time,
            peak_width_seconds: peak_width,
        }
    }

    /// Estimate the peak width of elutions within the quad.
    ///
    /// This is a pretty dumb algorithm ... basically...
    /// 1. picks 100 equally spaced points between the 20% and 80% of the frames.
    /// 2. At each of those points finds the top 10 peaks.
    /// 3. For each of them iteratively goes forward and backwards in time until
    ///    the intensity of the peak drops under 1% of its intensity.
    /// 4. The time difference between the forward and backward point is the peak width.
    /// 5. Finds the median of the peak widths.
    fn estimate_peak_width(frameslices: &[ExpandedFrameSlice<SortedState>]) -> Option<f64> {
        let start_idx = frameslices.len() / 5;
        let end_idx = frameslices.len() * 4 / 5;
        let step = (end_idx - start_idx) / 100;

        if step == 0 {
            warn!("Estimating peak width failed, step is 0");
            return None;
        }

        #[derive(Debug, Clone)]
        struct Peak {
            tof: u32,
            patience: u32,
            max_rt: f64,
            min_rt: f64,
            intensity: u32,
            fwdone: bool,
            bwdone: bool,
            any_update: bool,
        }

        let mut peaks: Vec<Peak> = Vec::with_capacity(500);

        for i in 0..49 {
            let local_idx = start_idx + i * step;

            let local_frame = &frameslices[local_idx];
            let local_rt = local_frame.rt;
            let (top_intens, top_indices) = top_n(&local_frame.intensities, 10);
            let tofs: Vec<_> = top_indices
                .iter()
                .zip(top_intens.iter())
                .filter_map(|(tof, inten)| if *inten > 100u32 { Some(tof) } else { None })
                .map(|x| local_frame.tof_indices[*x])
                .collect();
            let tofs: Vec<u32> = tofs
                .as_slice()
                .windows(2)
                .filter_map(|x| {
                    let a = x[0];
                    let b = x[1];
                    if a.abs_diff(b) > 5 { Some(a) } else { None }
                })
                .collect();

            let mut local_peaks: Vec<Peak> = tofs
                .iter()
                .map(|tof| {
                    let mut local_inten = 0u32;
                    local_frame.query_peaks(IncludedRange::new(tof - 1, tof + 1), None, &mut |x| {
                        local_inten += x.intensity
                    });
                    assert!(local_inten > 0);
                    Peak {
                        tof: *tof,
                        patience: 1,
                        max_rt: local_rt,
                        min_rt: local_rt,
                        intensity: local_inten,
                        fwdone: false,
                        bwdone: false,
                        any_update: false,
                    }
                })
                .collect();

            // Forward lookup
            for lookup_frame in frameslices.iter().skip(local_idx + 1) {
                let local_rt = lookup_frame.rt;
                let mut any = false;
                for peak in local_peaks.iter_mut() {
                    if peak.fwdone {
                        continue;
                    }
                    let mut local_int = 0u32;
                    lookup_frame.query_peaks(
                        IncludedRange::new(peak.tof - 1, peak.tof + 1),
                        None,
                        &mut |x| local_int += x.intensity,
                    );
                    if local_int > (peak.intensity / 100) {
                        peak.max_rt = local_rt;
                        peak.patience = 1;
                        peak.any_update = true;
                    } else if peak.patience > 0 {
                        peak.patience -= 1;
                    } else {
                        peak.fwdone = true;
                    }

                    any = true;
                }
                if !any {
                    break;
                }
            }

            // reset patience
            for peak in local_peaks.iter_mut() {
                peak.patience = 1;
            }

            // Backward lookup
            let mut j = local_idx - 1;
            while j > 1 {
                // for j in (0..local_idx).rev() {
                let lookup_frame = &frameslices[j];
                let local_rt = lookup_frame.rt;
                let mut any = false;
                for peak in local_peaks.iter_mut() {
                    if peak.bwdone {
                        continue;
                    }
                    let mut local_int = 0;
                    lookup_frame.query_peaks(
                        IncludedRange::new(peak.tof - 1, peak.tof + 1),
                        None,
                        &mut |x| local_int += x.intensity,
                    );

                    if local_int > (peak.intensity / 100) {
                        peak.min_rt = local_rt;
                        peak.patience = 1;
                        peak.any_update = true;
                    } else if peak.patience > 0 {
                        peak.patience -= 1;
                    } else {
                        peak.bwdone = true;
                    }

                    any = true;
                }
                if !any {
                    break;
                }
                j -= 1;
            }

            for peak in local_peaks.into_iter() {
                if peak.fwdone && peak.bwdone && peak.any_update {
                    peaks.push(peak);
                }
            }
        }

        if peaks.is_empty() {
            return None;
        }

        // get median
        peaks.sort_by(|x, y| {
            (x.max_rt - x.min_rt)
                .partial_cmp(&(y.max_rt - y.min_rt))
                .unwrap()
        });
        let median = peaks[peaks.len() / 2].max_rt - peaks[peaks.len() / 2].min_rt;
        Some(median)
    }
}

#[instrument(skip(frames))]
pub fn par_expand_and_centroid_frames(
    frames: impl ParallelIterator<Item = Frame>,
    ims_tol_pct: f64,
    mz_tol_ppm: f64,
    window_width: usize,
    ims_converter: &Scan2ImConverter,
    mz_converter: &Tof2MzConverter,
) -> HashMap<Option<SingleQuadrupoleSetting>, Vec<ExpandedFrameSlice<SortedState>>> {
    let split_frames = par_expand_and_arrange_frames(frames);

    let out: HashMap<Option<SingleQuadrupoleSetting>, Vec<ExpandedFrameSlice<SortedState>>> =
        split_frames
            .into_iter()
            .map(|(qs, frameslices)| {
                // NOTE: Since the the centroiding runs in paralel over the windows, its ok if this
                // outer loop is done in series.
                let start_peaks: usize = frameslices.iter().map(|x| x.len()).sum();
                let centroided = par_lazy_centroid_frameslices(
                    &frameslices,
                    window_width,
                    ims_tol_pct,
                    mz_tol_ppm,
                    ims_converter,
                    mz_converter,
                );
                let end_peaks: usize = centroided.iter().map(|x| x.len()).sum();
                trace!(
                    "Peak counts for quad {:?}: raw={}/centroid={}",
                    qs, start_peaks, end_peaks
                );
                (qs, centroided)
            })
            .collect();

    out
}

fn centroid_frameslice_window(
    frameslices: &[ExpandedFrameSlice<SortedState>],
    ims_tol_pct: f64,
    mz_tol_ppm: f64,
    ims_converter: &Scan2ImConverter,
    mz_converter: &Tof2MzConverter,
) -> ExpandedFrameSlice<SortedState> {
    assert!(frameslices.len() > 1, "Expected at least 2 frameslices");
    let reference_index = frameslices.len() / 2;

    let peak_refs: Vec<PeakArrayRefs> = frameslices
        .iter()
        .map(|x| PeakArrayRefs::new(&x.tof_indices, &x.scan_numbers, &x.intensities))
        .collect();

    let ((tof_array, intensity_array), ims_array) = lazy_centroid_weighted_frame(
        &peak_refs,
        reference_index,
        |tof| tof_tol_range(tof, mz_tol_ppm, mz_converter),
        |scan| scan_tol_range(scan, ims_tol_pct, ims_converter),
    );

    ExpandedFrameSlice {
        tof_indices: tof_array,
        scan_numbers: ims_array,
        intensities: intensity_array,
        frame_index: frameslices[reference_index].frame_index,
        rt: frameslices[reference_index].rt,
        acquisition_type: frameslices[reference_index].acquisition_type,
        ms_level: frameslices[reference_index].ms_level,
        quadrupole_settings: frameslices[reference_index].quadrupole_settings,
        intensity_correction_factor: frameslices[reference_index].intensity_correction_factor,
        window_group: frameslices[reference_index].window_group,
        window_subindex: frameslices[reference_index].window_subindex,
        _sorting_state: PhantomData::<SortedState>,
    }
}

#[instrument(skip(frameslices))]
pub fn par_lazy_centroid_frameslices(
    frameslices: &[ExpandedFrameSlice<SortedState>],
    window_width: usize,
    ims_tol_pct: f64,
    mz_tol_ppm: f64,
    ims_converter: &Scan2ImConverter,
    mz_converter: &Tof2MzConverter,
) -> Vec<ExpandedFrameSlice<SortedState>> {
    assert!(
        frameslices.windows(2).all(|window| {
            let a = &window[0];
            let b = &window[1];
            a.rt < b.rt && a.quadrupole_settings == b.quadrupole_settings
        }),
        "All frames should be sorted by rt and have the same quad settings"
    );

    assert!(frameslices.len() > window_width);

    let local_lambda = |fss: &[ExpandedFrameSlice<SortedState>]| {
        centroid_frameslice_window(fss, ims_tol_pct, mz_tol_ppm, ims_converter, mz_converter)
    };

    frameslices
        .par_windows(window_width)
        .map(local_lambda)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_expanded_frame_from_frame() {
        let frame = Frame {
            tof_indices: vec![1, 2, 3, 4, 1, 2],
            scan_offsets: vec![0, 2, 4, 6],
            intensities: vec![10, 20, 30, 40, 50, 60],
            index: 0,
            rt_in_seconds: 0.0,
            acquisition_type: AcquisitionType::DIAPASEF,
            ms_level: MSLevel::MS1,
            quadrupole_settings: Arc::new(QuadrupoleSettings::default()),
            intensity_correction_factor: 1.0,
            window_group: 0,
        };

        let expanded_frame = ExpandedFrame::from_frame(frame);
        assert_eq!(expanded_frame.tof_indices, vec![1, 1, 2, 2, 3, 4]);
        assert_eq!(expanded_frame.scan_numbers, vec![0, 2, 0, 2, 1, 1]);
        assert_eq!(expanded_frame.intensities, vec![10, 50, 20, 60, 30, 40]);
        assert_eq!(expanded_frame.frame_index, 0);
        assert_eq!(expanded_frame.rt, 0.0);
        assert_eq!(expanded_frame.acquisition_type, AcquisitionType::DIAPASEF);
        assert_eq!(expanded_frame.ms_level, MSLevel::MS1);
        assert_eq!(
            *expanded_frame.quadrupole_settings,
            QuadrupoleSettings::default()
        );
        assert_eq!(expanded_frame.intensity_correction_factor, 1.0);
        assert_eq!(expanded_frame.window_group, 0);
    }
}
