use std::iter::repeat;

use super::expanded_frame::{
    ExpandedFrameSlice,
    SortingStateTrait,
};
use super::single_quad_settings::SingleQuadrupoleSetting;
use timsrust::{
    AcquisitionType,
    MSLevel,
};

pub struct ExpandedWindowGroup {
    pub tof_indices: Vec<u32>,
    pub scan_numbers: Vec<usize>,
    pub intensities: Vec<u32>,
    pub retention_times: Vec<f64>,
    // pub frame_indices: Vec<usize>,
    pub acquisition_type: AcquisitionType,
    pub ms_level: MSLevel,
    pub quadrupole_settings: Option<SingleQuadrupoleSetting>,
    // pub intensity_correction_factor: f64,
    pub window_group: u8,
}

impl ExpandedWindowGroup {
    pub fn from_expanded_frame_slices<S: SortingStateTrait>(
        expanded_frame_slices: Vec<ExpandedFrameSlice<S>>,
    ) -> Self {
        let num_peaks = expanded_frame_slices
            .iter()
            .map(|x| x.tof_indices.len())
            .sum();
        let mut tof_indices = Vec::with_capacity(num_peaks);
        let mut scan_numbers = Vec::with_capacity(num_peaks);
        let mut intensities = Vec::with_capacity(num_peaks);
        let mut retention_times = Vec::with_capacity(num_peaks);

        let acq_type = expanded_frame_slices[0].acquisition_type;
        let ms_level = expanded_frame_slices[0].ms_level;
        let quad_settings = expanded_frame_slices[0].quadrupole_settings;
        let window_group = expanded_frame_slices[0].window_group;

        for slice in expanded_frame_slices {
            let local_len = slice.tof_indices.len();
            tof_indices.extend(slice.tof_indices);
            scan_numbers.extend(slice.scan_numbers);
            intensities.extend(slice.intensities);
            retention_times.extend(repeat(slice.rt).take(local_len));
        }

        ExpandedWindowGroup {
            tof_indices,
            scan_numbers,
            intensities,
            retention_times,
            acquisition_type: acq_type,
            ms_level,
            quadrupole_settings: quad_settings,
            // intensity_correction_factor: expanded_frame_slices[0].intensity_correction_factor,
            window_group,
        }
    }
}
