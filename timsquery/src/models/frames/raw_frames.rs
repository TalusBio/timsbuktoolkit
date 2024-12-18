use timsrust::{
    Frame,
    QuadrupoleSettings,
};

use super::raw_peak::RawPeak;
use crate::utils::tolerance_ranges::IncludedRange;
use tracing::trace;

pub fn scans_matching_quad(
    quad_settings: &QuadrupoleSettings,
    quad_range: IncludedRange<f64>,
) -> Option<IncludedRange<usize>> {
    let mut min_start = usize::MAX;
    let mut max_end = 0;

    for i in 0..quad_settings.scan_ends.len() {
        let half_width = quad_settings.isolation_width[i] / 2.0;
        let quad_start = quad_settings.isolation_mz[i] - half_width;
        let quad_end = quad_settings.isolation_mz[i] + half_width;

        if quad_start <= quad_range.end() && quad_end >= quad_range.start() {
            let start = quad_settings.scan_starts[i];
            let end = quad_settings.scan_ends[i];
            min_start = min_start.min(start);
            max_end = max_end.max(end);
        }
    }

    if min_start == usize::MAX {
        None
    } else {
        Some((min_start, max_end).into())
    }
}

pub fn frame_elems_matching(
    frame: &Frame,
    tof_range: IncludedRange<u32>,
    scan_range: IncludedRange<usize>,
    quad_range: Option<IncludedRange<f64>>,
) -> impl Iterator<Item = RawPeak> + '_ {
    trace!(
        "frame_elems_matching tof_range: {:?}, scan_range: {:?}, quad_range: {:?}",
        tof_range, scan_range, quad_range
    );
    let quad_scan_range = quad_range
        .and_then(|quad_range| scans_matching_quad(&frame.quadrupole_settings, quad_range));

    let scan_range_use = if quad_scan_range.is_none() {
        0..0
    } else {
        let quad_scan_range = quad_scan_range.unwrap();

        // Only checkinghere bc its common for them to get flipped
        // bc skill issues. (and the highest scan is actually the lowest 1/k0)
        let min_scan = scan_range.start();
        let max_scan = scan_range.end();

        let min_quad_scan = quad_scan_range.start();
        let max_quad_scan = quad_scan_range.end();

        let (scan_ind_start, scan_ind_end): (usize, usize) =
            (min_quad_scan.max(min_scan), max_quad_scan.min(max_scan));

        if scan_ind_end < scan_ind_start {
            // This can happen if the quad range matches but not the scan range.
            // TODO refactor this logic...
            0..0
        } else {
            scan_ind_start..scan_ind_end
        }
    };

    let retention_time = frame.rt as f32;

    scan_range_use.flat_map(move |scan_index| {
        let scan_is = frame.scan_offsets[scan_index];
        let scan_ie = frame.scan_offsets[scan_index + 1];

        let outs: Vec<RawPeak> = (scan_is..scan_ie)
            .map(move |peak_index| {
                let tof_ind = frame.tof_indices[peak_index];
                let intensity = frame.intensities[peak_index];

                (tof_ind, intensity, scan_index)
            })
            .filter(|(tof_ind, _, _)| *tof_ind >= tof_range.start() && *tof_ind < tof_range.end())
            .map(|(tof_ind, intensity, scan_index)| RawPeak {
                scan_index,
                tof_index: tof_ind,
                intensity,
                retention_time,
            })
            .collect();

        outs
    })
}
