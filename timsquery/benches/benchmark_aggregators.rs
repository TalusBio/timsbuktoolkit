use std::time::{
    Duration,
    Instant,
};
use timsquery::Aggregator;
use timsquery::models::aggregators::raw_peak_agg::multi_chromatogram_agg::aggregator::ParitionedCMGAggregator;
use timsquery::models::aggregators::{
    MultiCMGStatsAgg,
    RawPeakIntensityAggregator,
    RawPeakVectorAggregator,
};
use timsquery::models::frames::raw_peak::RawPeak;

fn setup_peaks() -> Vec<RawPeak> {
    (1..=10_000)
        .map(|i| RawPeak {
            scan_index: i,
            tof_index: i as u32,
            intensity: i as u32,
            retention_time: i as f32,
        })
        .collect()
}

fn bench_paritioned_cmg_agg() -> Duration {
    // RN this is glaciation-slow
    // For sure its is worth optimizing this guy
    let peaks = setup_peaks();
    let mut tot = Duration::from_millis(0);

    for _ in 0..50 {
        let now = Instant::now();
        let mut agg = ParitionedCMGAggregator::new(vec![1usize], 400, vec![400]);
        for peak in peaks.iter() {
            agg.add(
                peak.retention_time as u32,
                peak.scan_index,
                peak.tof_index,
                peak.intensity as u64,
            );
        }
        let elapsed = now.elapsed();
        tot += elapsed;
        println!(
            "Elapsed Aggregating ParitionedCMGAggregator: {:.2?}",
            elapsed
        );
    }

    tot
}

fn bench_raw_peak_intensity_agg() -> Duration {
    let peaks = setup_peaks();
    let mut tot = Duration::from_millis(0);

    for _ in 0..50 {
        let now = Instant::now();
        let mut agg = RawPeakIntensityAggregator::new(42);
        for peak in peaks.iter() {
            agg.add(*peak);
        }
        let elapsed = now.elapsed();
        tot += elapsed;
        println!(
            "Elapsed Aggregating RawPeakIntensityAggregator: {:.2?}",
            elapsed
        );
    }

    tot
}

fn bench_raw_peak_vector_agg() -> Duration {
    let peaks = setup_peaks();
    let mut tot = Duration::from_millis(0);

    for _ in 0..50 {
        let now = Instant::now();
        let mut agg = RawPeakVectorAggregator::new(42);
        for peak in peaks.iter() {
            agg.add(*peak);
        }
        let elapsed = now.elapsed();
        tot += elapsed;
        println!(
            "Elapsed Aggregating RawPeakVectorAggregator: {:.2?}",
            elapsed
        );
    }

    tot
}

fn main() {
    let tot_bpc = bench_paritioned_cmg_agg();
    let tot_rpi = bench_raw_peak_intensity_agg();
    let tot_va = bench_raw_peak_vector_agg();

    println!(
        "Total: {:.2?} ParitionedCMGAggregator: {:.2?} RawPeakIntensityAggregator: {:.2?} RawPeakVectorAggregator: {:.2?}",
        tot_bpc + tot_rpi + tot_va,
        tot_bpc,
        tot_rpi,
        tot_va
    );
}
