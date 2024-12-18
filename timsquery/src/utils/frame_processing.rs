use crate::sort_vecs_by_first;
use tracing::error;

use super::tolerance_ranges::IncludedRange;

pub type TofIntensityVecs = (Vec<u32>, Vec<u32>);
pub type CentroidedVecs = (TofIntensityVecs, Vec<usize>);

fn sort_n_check(
    agg_intensity: Vec<u32>,
    agg_tof: Vec<u32>,
    agg_ims: Vec<usize>,
) -> ((Vec<u32>, Vec<u32>), Vec<usize>) {
    let (tof_array, ims_array, intensity_array) =
        sort_vecs_by_first!(agg_tof, agg_ims, agg_intensity);

    if let Some(x) = intensity_array.last() {
        assert!(
            *x > 0,
            "Expected all intensities to be positive non-zero, got {:?}",
            intensity_array
        );
        let max_tof = tof_array.iter().max().expect("At least one element");
        assert!(*max_tof < (u32::MAX - 1));
    }
    ((tof_array, intensity_array), ims_array)
}

pub struct PeakArrayRefs<'a> {
    pub tof_array: &'a [u32],
    pub ims_array: &'a [usize],
    pub intensity_array: &'a [u32],
}

impl<'a> PeakArrayRefs<'a> {
    pub fn new(
        tof_array: &'a [u32],
        ims_array: &'a [usize],
        intensity_array: &'a [u32],
    ) -> PeakArrayRefs<'a> {
        // TODO make the asserts optional at compile time.
        assert!(tof_array.len() == intensity_array.len());
        assert!(tof_array.len() == ims_array.len());
        assert!(
            tof_array.windows(2).all(|x| x[0] <= x[1]),
            "Expected tof array to be sorted"
        );
        Self {
            tof_array,
            ims_array,
            intensity_array,
        }
    }

    fn len(&self) -> usize {
        self.tof_array.len()
    }
}

struct Peak {
    center_tof: u32,
    center_ims: usize,
    // max_intensity: u32,
    agg_intensity: u32,
    // point_count: usize,
}
struct OrderItem {
    intensity: u32,
    weight: u32,
    tof: u32,
    ims: usize,
}

fn find_gaussian_peaks(
    items: &[OrderItem],
    min_points: usize,
    min_agg_intensity: u32,
    tof_tol_range_fn: impl Fn(u32) -> IncludedRange<u32>,
    ims_tol_range_fn: impl Fn(usize) -> IncludedRange<usize>,
) -> Vec<Peak> {
    let mut taken = vec![false; items.len()];
    let mut taken_queue = Vec::new();
    let mut peaks = Vec::new();
    let mut i = 0;

    while i < items.len() {
        let current = &items[i];
        if taken[i] {
            i += 1;
            continue;
        }

        // Fast pre-filtering: check if this point could be a peak center
        // I could make a smarter version of this ... like to get the cumulative
        // sum of intensities and then get the intensity of each window very quickly.
        // if current.weight < min_apex_intensity {
        //     i += 1;
        //     continue;
        // }

        // Get tolerance ranges for current point
        let tof_range = tof_tol_range_fn(current.tof);
        let ims_range = ims_tol_range_fn(current.ims);

        // Find all points within tolerance ranges
        let st = tof_range.start();
        let ep = tof_range.end();
        let window_start = items[..i].partition_point(|item| item.tof < st);
        let window_end = items[i..].partition_point(|item| item.tof < ep) + i;

        // Count points and check intensity distribution
        let mut point_count = 0;
        let mut is_peak = true;
        let max_weight = current.weight;
        let mut agg_intensity = 0;

        for j in window_start..window_end {
            if taken[j] {
                continue;
            };
            let point = &items[j];

            if ims_range.contains(point.ims) {
                taken_queue.push(j);
                point_count += 1;
                agg_intensity += point.intensity;

                // Track maximum intensity
                if point.weight > max_weight {
                    is_peak = false; // Current point isn't the peak if we found higher intensity
                    break;
                }
            }
        }

        // If this is a valid peak, record it
        if is_peak && agg_intensity >= min_agg_intensity && point_count >= min_points {
            peaks.push(Peak {
                center_tof: current.tof,
                center_ims: current.ims,
                agg_intensity,
            });
            for j in taken_queue.iter() {
                taken[*j] = true;
            }
            taken_queue.clear();
        } else {
            // If this wasn't a peak, jump to the maximum intensity point we found
            // This optimization helps us avoid checking every point in dense regions
            // NOTE: This optimization might make sense but in our data might skip
            // peaks that are very close in tof but not in mobility.
            //
            // i = if max_intensity_idx > i {
            //     max_intensity_idx
            // } else {
            //     i + 1
            // };
            i += 1;
            taken_queue.clear();
        }
    }

    peaks
}

// TODO: Instrument to add/return aggregation metrics/stats
pub fn lazy_centroid_weighted_frame<'a>(
    peak_refs: &'a [PeakArrayRefs<'a>],
    reference_index: usize,
    tof_tol_range_fn: impl Fn(u32) -> IncludedRange<u32>,
    ims_tol_range_fn: impl Fn(usize) -> IncludedRange<usize>,
) -> CentroidedVecs {
    let slice_sizes: Vec<usize> = peak_refs.iter().map(|x| x.len()).collect();
    let tot_size: usize = slice_sizes.iter().sum();
    let ref_arrays = &peak_refs[reference_index];
    let num_arrays = peak_refs.len();
    assert!(num_arrays > 1);
    let arr_len = ref_arrays.len();
    if arr_len == 0 {
        error!("No peaks in reference array when centroiding");
        return ((Vec::new(), Vec::new()), Vec::new());
    }

    let mut order: Vec<OrderItem> = Vec::with_capacity(tot_size);
    for (major_idx, peak_ref) in peak_refs.iter().enumerate() {
        for (minor_idx, &intensity) in peak_ref.intensity_array.iter().enumerate() {
            let tof = peak_ref.tof_array[minor_idx];
            let ims = peak_ref.ims_array[minor_idx];
            let weight = if major_idx == reference_index {
                peak_ref.intensity_array[minor_idx]
            } else {
                0
            };
            order.push(OrderItem {
                tof,
                intensity: weight,
                ims,
                weight: intensity,
            });
        }
    }

    // Sort by intensity
    order.sort_unstable_by(|a, b| a.tof.cmp(&b.tof));
    assert!(order[tot_size - 1].tof > order[0].tof);

    let outpeaks = find_gaussian_peaks(&order, 1, 100, &tof_tol_range_fn, ims_tol_range_fn);

    let ((agg_intensity, agg_tof), agg_ims) = outpeaks
        .into_iter()
        .map(|x| ((x.agg_intensity, x.center_tof), x.center_ims))
        .unzip();

    sort_n_check(agg_intensity, agg_tof, agg_ims)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper function to create a simple tolerance range for testing
    fn test_tof_tolerance(tof: u32) -> IncludedRange<u32> {
        let tolerance = 2;
        (tof.saturating_sub(tolerance), tof.saturating_add(tolerance)).into()
    }

    fn test_ims_tolerance(ims: usize) -> IncludedRange<usize> {
        let tolerance = 1;
        (ims.saturating_sub(tolerance), ims.saturating_add(tolerance)).into()
    }

    // Helper function to create PeakArrayRefs
    fn create_peak_refs<'a>(
        tof: &'a [u32],
        ims: &'a [usize],
        intensity: &'a [u32],
    ) -> PeakArrayRefs<'a> {
        PeakArrayRefs::new(tof, ims, intensity)
    }

    #[test]
    fn test_basic_centroid() {
        // Test case with two simple peaks that should be merged
        let peak_refs = vec![
            create_peak_refs(
                &[100, 101, 200, 201],   // tof
                &[1, 1, 2, 2],           // ims
                &[1000, 2000, 500, 600], // intensity
            ),
            create_peak_refs(
                &[101, 201], // tof
                &[1, 2],     // ims
                &[900, 600], // intensity
            ),
        ];

        let result = lazy_centroid_weighted_frame(
            &peak_refs,
            0, // reference_index
            test_tof_tolerance,
            test_ims_tolerance,
        );

        let ((tof_array, intensity_array), ims_array) = result;

        assert_eq!(
            tof_array.len(),
            2,
            "Expected 2 centroids, got {:?}, {:?}, {:?}",
            tof_array,
            intensity_array,
            ims_array
        );
        assert_eq!(intensity_array.len(), 2);
        assert_eq!(ims_array.len(), 2);

        // Check that the peaks were properly merged and weighted
        assert!(tof_array[0] >= 100 && tof_array[0] <= 101);
        assert_eq!(ims_array[0], 1);
        assert_eq!(intensity_array[0], 3000); // Should be sum of intensities
        assert_eq!(intensity_array[1], 1100); // Should be sum of intensities
    }

    // // NOTE: Current (experimental) implementation doe not have limit
    //
    // #[test]
    // fn test_max_peaks_limit() {
    //     // Test that the function respects the max_peaks parameter
    //     let peak_refs = vec![
    //         create_peak_refs(&[100, 200, 300], &[1, 2, 3], &[1000, 900, 800]),
    //         create_peak_refs(&[101, 201, 301], &[1, 2, 3], &[950, 850, 750]),
    //     ];

    //     let result = lazy_centroid_weighted_frame(
    //         &peak_refs,
    //         0,
    //         2, // max_peaks set to 2
    //         test_tof_tolerance,
    //         test_ims_tolerance,
    //     );

    //     let ((tof_array, _), _) = result;
    //     assert_eq!(
    //         tof_array.len(),
    //         2,
    //         "Should only return max_peaks number of peaks"
    //     );
    // }

    #[test]
    fn test_empty_input() {
        // Test handling of empty arrays
        let peak_refs = vec![
            create_peak_refs(&[], &[], &[]),
            create_peak_refs(&[], &[], &[]),
        ];

        let result =
            lazy_centroid_weighted_frame(&peak_refs, 0, test_tof_tolerance, test_ims_tolerance);

        let ((tof_array, intensity_array), ims_array) = result;
        assert_eq!(tof_array.len(), 0);
        assert_eq!(intensity_array.len(), 0);
        assert_eq!(ims_array.len(), 0);
    }

    #[test]
    fn test_intensity_weighted_centroid() {
        // Test that centroids are properly weighted by intensity
        let peak_refs = vec![
            create_peak_refs(
                &[100, 200],
                &[1, 2],
                &[1000, 100], // High intensity for first peak
            ),
            create_peak_refs(
                &[102, 202],
                &[1, 2],
                &[100, 1000], // High intensity for second peak
            ),
        ];

        let result =
            lazy_centroid_weighted_frame(&peak_refs, 0, test_tof_tolerance, test_ims_tolerance);

        let ((tof_array, _), ims_array) = result;

        // First centroid should be closer to 100 due to higher intensity
        assert!(tof_array[0] - 100 < 102 - tof_array[0]);
        assert_eq!(ims_array[0], 1);
    }

    #[test]
    #[should_panic]
    fn test_invalid_reference_index() {
        let peak_refs = vec![
            create_peak_refs(&[100], &[1], &[1000]),
            create_peak_refs(&[101], &[1], &[900]),
        ];

        // This should panic due to invalid reference index
        lazy_centroid_weighted_frame(
            &peak_refs,
            2, // Invalid index
            test_tof_tolerance,
            test_ims_tolerance,
        );
    }

    // Pretty decent test suggested by claude but im not sure if I really
    // need to suuport it ... since scan indices are usually less than 2_000
    // #[test]
    // fn test_large_input_values() {
    //     // Test handling of large values close to u32::MAX
    //     let peak_refs = vec![
    //         create_peak_refs(
    //             &[u32::MAX - 10, u32::MAX - 5],
    //             &[usize::MAX - 10, usize::MAX - 5],
    //             &[1000, 900],
    //         ),
    //         create_peak_refs(
    //             &[u32::MAX - 9, u32::MAX - 4],
    //             &[usize::MAX - 9, usize::MAX - 4],
    //             &[950, 850],
    //         ),
    //     ];

    //     let result =
    //         lazy_centroid_weighted_frame(&peak_refs, 0, 10, test_tof_tolerance, test_ims_tolerance);

    //     let ((tof_array, _), _) = result;
    //     assert!(!tof_array.is_empty(), "Should handle large values properly");
    // }
}
