use half::f16;
use timscentroid::indexing::{
    IndexedPeak,
    IndexedPeakGroup,
    IndexedTimstofPeaks,
};
use timscentroid::lazy::LazyIndexedTimstofPeaks;
use timscentroid::utils::OptionallyRestricted::*;
use timscentroid::utils::TupleRange;

// For now, we'll skip MS2 geometry tests since QuadrupoleIsolationScheme
// doesn't have a public constructor. We can test MS1 queries instead.

#[test]
fn test_small_dataset_eager_vs_lazy_ms1() {
    // Create a small, controlled dataset
    let mut peaks = Vec::new();

    // Add 10 peaks across different m/z ranges
    for i in 0..10 {
        peaks.push(IndexedPeak {
            mz: 400.0 + (i as f32) * 10.0, // 400, 410, 420, ..., 490
            intensity: 100.0 * (i as f32 + 1.0),
            mobility_ook0: f16::from_f32(1.0 + (i as f32) * 0.01),
            cycle_index: i % 3, // Cycles 0, 1, 2
        });
    }

    let cycle_to_rt_ms = vec![0, 100, 200]; // 3 cycles
    let bucket_size = 4; // Small buckets: 4 peaks each

    // Build MS1 group using canonical constructor
    let (ms1_group, stats) = IndexedPeakGroup::testing_new(peaks, cycle_to_rt_ms, bucket_size);

    println!(
        "Created MS1 group: {} peaks, {} buckets",
        stats.num_peaks, stats.num_buckets
    );

    // Create index with empty MS2 (test MS1 only)
    let index = IndexedTimstofPeaks::from_parts(ms1_group, vec![]);

    // Serialize to temp directory
    let temp_dir = std::env::temp_dir().join("timscentroid_test_lazy_ms1");
    if temp_dir.exists() {
        std::fs::remove_dir_all(&temp_dir).unwrap();
    }

    let mut config = timscentroid::serialization::SerializationConfig::default();
    // Just some weird number ...
    config.row_group_size = 3;

    index
        .save_to_directory_with_config(&temp_dir, config)
        .unwrap();

    // Load eagerly
    let eager_index = IndexedTimstofPeaks::load_from_directory(&temp_dir).unwrap();

    // Load lazily
    let lazy_index = LazyIndexedTimstofPeaks::load_from_directory(&temp_dir).unwrap();

    // Query both with same parameters
    let mz_range = TupleRange::try_new(420.0, 460.0).unwrap(); // Should match peaks at 420, 430, 440, 450, 460

    // Original index results
    let orig_index_results = index
        .query_peaks_ms1(mz_range, Unrestricted, Unrestricted)
        .collect::<Vec<_>>();

    assert!(orig_index_results.iter().any(|p| p.mz == 420.0));
    assert!(orig_index_results.iter().any(|p| p.mz == 430.0));
    assert!(orig_index_results.iter().any(|p| p.mz == 440.0));
    assert!(orig_index_results.iter().any(|p| p.mz == 450.0));
    assert!(orig_index_results.iter().any(|p| p.mz == 460.0));
    assert_eq!(
        orig_index_results.len(),
        5,
        "Should find 5 peaks in original index"
    );

    // Query MS1 with eager
    let eager_results: Vec<_> = eager_index
        .query_peaks_ms1(mz_range, Unrestricted, Unrestricted)
        .collect();

    // Query MS1 with lazy
    let lazy_results: Vec<_> = lazy_index
        .query_peaks_ms1(mz_range, Unrestricted, Unrestricted)
        .collect();

    println!("Original index found {} peaks", orig_index_results.len());
    println!("Eager found {} peaks", eager_results.len());
    println!("Lazy found {} peaks", lazy_results.len());

    // Print details for debugging
    println!("\nEager peaks:");
    for peak in &eager_results {
        println!(
            "  mz={:.1}, intensity={:.1}, mobility={:.3}, cycle={}",
            peak.mz, peak.intensity, peak.mobility_ook0, peak.cycle_index
        );
    }

    println!("\nLazy peaks:");
    for peak in &lazy_results {
        println!(
            "  mz={:.1}, intensity={:.1}, mobility={:.3}, cycle={}",
            peak.mz, peak.intensity, peak.mobility_ook0, peak.cycle_index
        );
    }

    // Assert they match
    assert_eq!(
        eager_results.len(),
        lazy_results.len(),
        "Eager and lazy should return same number of peaks!"
    );
    assert_eq!(
        orig_index_results.len(),
        lazy_results.len(),
        "Lazy results should match original index!"
    );

    // Cleanup
    std::fs::remove_dir_all(&temp_dir).unwrap();
}

// I could mayube go nuclear and fuzz test it ... just generate random peak sets and ensure
// that the bucket calculations are consistent with each other AND with a naive greedy implementation
// with varying bucket sizes and row group sizes ...

#[test]
fn test_bucket_boundary_case() {
    // This test checks the bucket boundary issue
    // If we have 10 peaks and bucket_size=4:
    // Bucket 0: peaks 0-3 (mz 400-430)
    // Bucket 1: peaks 4-7 (mz 440-470)
    // Bucket 2: peaks 8-9 (mz 480-490)
    //
    // Row group with row_group_size=2M will contain all buckets
    // But the bucket calculation might have off-by-one errors

    let peaks: Vec<_> = (0..10)
        .map(|i| IndexedPeak {
            mz: 400.0 + (i as f32) * 10.0,
            intensity: 100.0,
            mobility_ook0: f16::from_f32(1.0),
            cycle_index: 0,
        })
        .collect();

    let cycle_to_rt_ms = vec![0];
    let bucket_size = 4;

    let (group, stats) = IndexedPeakGroup::testing_new(peaks, cycle_to_rt_ms, bucket_size);

    // Should have 3 buckets: [0-3], [4-7], [8-9]
    assert_eq!(stats.num_buckets, 3, "Should have 3 buckets");
    assert_eq!(stats.num_peaks, 10, "Should have 10 peaks");

    println!(
        "Test passed: {} peaks organized into {} buckets",
        stats.num_peaks, stats.num_buckets
    );
}
