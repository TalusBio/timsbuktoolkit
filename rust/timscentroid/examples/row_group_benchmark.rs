/// Benchmark different row group sizes for Parquet serialization
///
/// Compares disk usage, load times, and query efficiency across different row group sizes
/// to find the optimal balance for lazy loading from cloud storage.
use half::f16;
use rand::SeedableRng;
use timscentroid::lazy::LazyIndexedTimstofPeaks;
use timscentroid::serialization::SerializationConfig;
use timscentroid::utils::{
    OptionallyRestricted,
    TupleRange,
};
use timscentroid::{
    CentroidingConfig,
    IndexedTimstofPeaks,
    StorageLocation,
};
use timsrust::TimsTofPath;

use OptionallyRestricted::{
    Restricted,
    Unrestricted,
};

const NUM_QUERIES: usize = 500;

fn main() {
    const DATA_FILE: &str =
        "/Users/sebastianpaez/data/decompressed_timstof/250225_Desnaux_200ng_Hela_ICC_on_DIA.d/";
    const RANDOM_SEED: u64 = 42;

    let row_group_sizes = vec![
        ("4K", 4_096),
        ("10K", 10_000),
        ("50K", 50_000),
        ("100K", 100_000),
        ("500K", 500_000),
        ("1M", 1_000_000),
    ];

    println!(
        "Row Group Size Benchmark - {} queries per config",
        NUM_QUERIES
    );
    println!("Data: {}\n", DATA_FILE);

    // Load and index data
    let file = TimsTofPath::new(DATA_FILE).expect("Failed to open TimsTOF file");
    let config = CentroidingConfig {
        max_peaks: 20_000,
        mz_ppm_tol: 5.0,
        im_pct_tol: 3.0,
        early_stop_iterations: 200,
    };

    let start = std::time::Instant::now();
    let (index, stats) = IndexedTimstofPeaks::from_timstof_file(&file, config);
    println!("Indexed in {:?}: {}\n", start.elapsed(), stats);

    // Generate deterministic queries
    let queries = {
        let mut rng = rand::rngs::StdRng::seed_from_u64(RANDOM_SEED);
        build_queries(&mut rng, NUM_QUERIES)
    };

    // Test each row group size
    let mut results = Vec::new();
    for (label, row_group_size) in &row_group_sizes {
        println!("Testing {} ({} peaks/group)", label, row_group_size);
        let result = test_row_group_size(&index, *row_group_size, label, &queries);
        results.push((label.to_string(), *row_group_size, result));
    }

    print_results(&results);
}

#[derive(Debug, Clone)]
struct BenchmarkResult {
    disk_mb: f64,
    eager_load_ms: f64,
    lazy_init_ms: f64,
    query_ms: f64,
    bytes_per_peak: f64,
    num_gets: usize,
}

fn test_row_group_size(
    index: &IndexedTimstofPeaks,
    row_group_size: usize,
    label: &str,
    queries: &[((f32, f32), (f32, f32), (f16, f16))],
) -> BenchmarkResult {
    let output_dir = std::path::PathBuf::from(format!("./benchmark_rg_{}", label));
    if output_dir.exists() {
        std::fs::remove_dir_all(&output_dir).unwrap();
    }

    // Serialize
    let config = SerializationConfig {
        compression: parquet::basic::Compression::SNAPPY,
        row_group_size,
        write_batch_size: 8192,
    };
    index
        .save_to_directory_with_config(&output_dir, config)
        .unwrap();

    let disk_mb = calculate_directory_size(&output_dir).unwrap() as f64 / 1_048_576.0;

    // Eager load
    let start = std::time::Instant::now();
    let location = StorageLocation::from_path(&output_dir);
    let _eager = IndexedTimstofPeaks::load_from_storage(location.clone()).unwrap();
    let eager_load_ms = start.elapsed().as_secs_f64() * 1000.0;

    // Lazy load + queries
    let start = std::time::Instant::now();
    let lazy_index = LazyIndexedTimstofPeaks::load_from_storage(location)
        .unwrap()
        .with_instrumentation(format!("RG_{}", label));
    let lazy_init_ms = start.elapsed().as_secs_f64() * 1000.0;

    let query_start = std::time::Instant::now();
    let mut total_peaks = 0;
    for (prec, frag, im) in queries {
        let prec = TupleRange::try_new(prec.0, prec.1).unwrap();
        let frag = TupleRange::try_new(frag.0, frag.1).unwrap();
        let im = TupleRange::try_new(im.0, im.1).unwrap();
        let results = lazy_index.query_peaks_ms2(prec, frag, Unrestricted, Restricted(im));
        for (_wg, peaks_vec) in results {
            total_peaks += peaks_vec.len();
        }
    }
    let query_ms = query_start.elapsed().as_secs_f64() * 1000.0 / queries.len() as f64;

    let metrics = lazy_index.metrics().unwrap().snapshot();
    let bytes_per_peak = metrics.bytes_read as f64 / total_peaks.max(1) as f64;

    println!(
        "  {:.1} MB, eager: {:.0}ms, lazy: {:.0}ms, query: {:.2}ms, {:.0} B/peak, {} GETs",
        disk_mb, eager_load_ms, lazy_init_ms, query_ms, bytes_per_peak, metrics.get_count
    );

    std::fs::remove_dir_all(&output_dir).unwrap();

    BenchmarkResult {
        disk_mb,
        eager_load_ms,
        lazy_init_ms,
        query_ms,
        bytes_per_peak,
        num_gets: metrics.get_count,
    }
}

fn print_results(results: &[(String, usize, BenchmarkResult)]) {
    println!(
        "\n{:<8} {:<10} {:<10} {:<10} {:<10} {:<12} {:<6}",
        "Size", "Disk(MB)", "Eager(ms)", "Lazy(ms)", "Query(ms)", "Bytes/Peak", "GETs"
    );
    println!("{}", "-".repeat(72));

    for (label, _size, r) in results {
        println!(
            "{:<8} {:<10.1} {:<10.0} {:<10.0} {:<10.2} {:<12.0} {:<6}",
            label,
            r.disk_mb,
            r.eager_load_ms,
            r.lazy_init_ms,
            r.query_ms,
            r.bytes_per_peak,
            r.num_gets
        );
    }

    // Find best configurations
    let best_disk = results
        .iter()
        .min_by(|a, b| a.2.disk_mb.partial_cmp(&b.2.disk_mb).unwrap())
        .unwrap();
    let best_query = results
        .iter()
        .min_by(|a, b| a.2.bytes_per_peak.partial_cmp(&b.2.bytes_per_peak).unwrap())
        .unwrap();

    println!(
        "\nBest disk usage: {} ({:.1} MB)",
        best_disk.0, best_disk.2.disk_mb
    );
    println!(
        "Best query efficiency: {} ({:.0} bytes/peak, {} GETs)",
        best_query.0, best_query.2.bytes_per_peak, best_query.2.num_gets
    );

    println!("\nS3 estimated (150ms/GET):");
    for (label, _size, r) in results {
        let s3_time = (r.num_gets as f64 * 150.0) + r.query_ms;
        println!("  {:<8} -> ~{:.0} ms/query", label, s3_time);
    }
}

fn build_queries(
    rng: &mut impl rand::Rng,
    num_queries: usize,
) -> Vec<((f32, f32), (f32, f32), (f16, f16))> {
    (0..num_queries)
        .map(|_| {
            let prec_start: f32 = rng.random_range(600.0..800.0);
            let frag_start: f32 = rng.random_range(600.0..800.0);
            let im_start = f16::from_f32(rng.random_range(0.7..1.1));
            (
                (prec_start, prec_start + 0.04),
                (frag_start, frag_start + 0.04),
                (im_start, im_start + f16::from_f32(0.1)),
            )
        })
        .collect()
}

fn calculate_directory_size(path: &std::path::Path) -> std::io::Result<u64> {
    let mut total_size = 0u64;
    for entry in std::fs::read_dir(path)? {
        let entry = entry?;
        let metadata = entry.metadata()?;
        if metadata.is_file() {
            total_size += metadata.len();
        } else if metadata.is_dir() {
            total_size += calculate_directory_size(&entry.path())?;
        }
    }
    Ok(total_size)
}
