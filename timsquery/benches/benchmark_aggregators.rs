
use timsquery::models::aggregators::{MultiCMGStatsAgg, RawPeakIntensityAggregator, RawPeakVectorAggregator};
use timsquery::models::aggregators::raw_peak_agg::multi_chromatogram_agg::aggregator::ParitionedCMGAggregator;
use timsquery::models::frames::raw_peak::RawPeak;
use std::time::Instant;
use timsquery::Aggregator;

fn setup_peaks() -> Vec<RawPeak> {
    (1..=10_000).map(|i| {
        RawPeak {
            scan_index: i,
            tof_index: i as u32,
            intensity: i as u32,
            retention_time: i as f32,
        }}).collect()
}

fn bench_paritioned_cmg_agg() {
    // RN this is glaciation-slow
    // For sure its is worth optimizing this guy
    let peaks = setup_peaks();

    for _ in 0..50 {
        let now = Instant::now();
        let mut agg = ParitionedCMGAggregator::new(vec![1usize], 400, vec![400]);
        for peak in peaks.iter() {
            agg.add(peak.retention_time as u32, peak.scan_index, peak.tof_index, peak.intensity as u64);
        }
        let elapsed = now.elapsed();
        println!("Elapsed Aggregating ParitionedCMGAggregator: {:.2?}", elapsed);
    }

}

fn bench_raw_peak_intensity_agg() {
    let peaks = setup_peaks();


    for _ in 0..50 {
        let now = Instant::now();
        let mut agg = RawPeakIntensityAggregator::new(42);
        for peak in peaks.iter() {
            agg.add(*peak);
        }
        let elapsed = now.elapsed();
        println!("Elapsed Aggregating RawPeakIntensityAggregator: {:.2?}", elapsed);
    }

}


fn bench_raw_peak_vector_agg() {
    let peaks = setup_peaks();

    for _ in 0..50 {
        let now = Instant::now();
        let mut agg = RawPeakVectorAggregator::new(42);
        for peak in peaks.iter() {
            agg.add(*peak);
        }
        let elapsed = now.elapsed();
        println!("Elapsed Aggregating RawPeakVectorAggregator: {:.2?}", elapsed);
    }

}

fn main() {
    bench_paritioned_cmg_agg();
    bench_raw_peak_intensity_agg();
    bench_raw_peak_vector_agg();
}