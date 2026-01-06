/// Compare local vs S3 lazy loading performance
use half::f16;
use timscentroid::StorageLocation;
use timscentroid::lazy::LazyIndexedTimstofPeaks;
use timscentroid::utils::{
    OptionallyRestricted,
    TupleRange,
};

use OptionallyRestricted::{
    Restricted,
    Unrestricted,
};

fn main() {
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .init();

    let serialized_dir = std::path::Path::new("./serialized_peaks");
    let s3_url = "s3://terraform-workstations-bucket/jspaezp/serialized_peaks_test";

    if !serialized_dir.exists() {
        eprintln!("Error: Run `cargo run --example serialization --release` first");
        std::process::exit(1);
    }

    println!("Lazy Loading Benchmark: Local vs S3\n");
    const NUM_ITERATIONS: usize = 5;

    // Benchmark init times
    println!("Local init:");
    let local_times = benchmark_init(serialized_dir, "Local", NUM_ITERATIONS);
    let local_avg = avg(&local_times);

    println!("\nS3 init ({}):", s3_url);
    let s3_times = benchmark_init_s3(s3_url, "S3", NUM_ITERATIONS);
    let s3_avg = avg(&s3_times);

    println!(
        "\nInit time - Local: {:.2?}, S3: {:.2?} (Local is {:.1}x faster)",
        local_avg,
        s3_avg,
        s3_avg.as_secs_f64() / local_avg.as_secs_f64()
    );

    // Load indices for query tests
    let index_local =
        LazyIndexedTimstofPeaks::load_from_storage(StorageLocation::from_path(serialized_dir))
            .unwrap()
            .with_instrumentation("Local");
    let index_s3 =
        LazyIndexedTimstofPeaks::load_from_storage(StorageLocation::from_url(s3_url).unwrap())
            .unwrap()
            .with_instrumentation("S3");

    println!("\nQuery tests:");
    compare_results(&index_local, &index_s3);
}

fn benchmark_init(dir: &std::path::Path, label: &str, n: usize) -> Vec<std::time::Duration> {
    (0..n)
        .map(|i| {
            let start = std::time::Instant::now();
            let _index =
                LazyIndexedTimstofPeaks::load_from_storage(StorageLocation::from_path(dir))
                    .unwrap()
                    .with_instrumentation(label);
            let elapsed = start.elapsed();
            println!("  {}: {:.2?}", i + 1, elapsed);
            elapsed
        })
        .collect()
}

fn benchmark_init_s3(url: &str, label: &str, n: usize) -> Vec<std::time::Duration> {
    (0..n)
        .map(|i| {
            let start = std::time::Instant::now();
            let _index =
                LazyIndexedTimstofPeaks::load_from_storage(StorageLocation::from_url(url).unwrap())
                    .unwrap()
                    .with_instrumentation(label);
            let elapsed = start.elapsed();
            println!("  {}: {:.2?}", i + 1, elapsed);
            elapsed
        })
        .collect()
}

fn avg(times: &[std::time::Duration]) -> std::time::Duration {
    times.iter().sum::<std::time::Duration>() / times.len() as u32
}

fn compare_results(index_local: &LazyIndexedTimstofPeaks, index_s3: &LazyIndexedTimstofPeaks) {
    let queries = build_queries(&mut rand::rng());

    // Sync queries
    println!("\n  Sync API:");
    let local_sync = test_querying(index_local, &queries);
    println!(
        "    Local: {} peaks, {:.2?}/query",
        local_sync.npeaks, local_sync.query_time
    );
    let s3_sync = test_querying(index_s3, &queries);
    println!(
        "    S3:    {} peaks, {:.2?}/query",
        s3_sync.npeaks, s3_sync.query_time
    );

    // Async queries
    println!("\n  Async API:");
    let runtime = tokio::runtime::Runtime::new().unwrap();
    let local_async = runtime.block_on(test_querying_async(index_local, &queries));
    println!(
        "    Local: {} peaks, {:.2?}/query",
        local_async.npeaks, local_async.query_time
    );
    let s3_async = runtime.block_on(test_querying_async(index_s3, &queries));
    println!(
        "    S3:    {} peaks, {:.2?}/query",
        s3_async.npeaks, s3_async.query_time
    );

    println!("\n  Speedup (async vs sync):");
    println!(
        "    Local: {:.1}x",
        local_sync.query_time.as_secs_f64() / local_async.query_time.as_secs_f64()
    );
    println!(
        "    S3:    {:.1}x",
        s3_sync.query_time.as_secs_f64() / s3_async.query_time.as_secs_f64()
    );

    index_local.print_metrics("Local");
    index_s3.print_metrics("S3");
}

#[derive(Debug)]
struct QueryResult {
    npeaks: usize,
    query_time: std::time::Duration,
}

fn test_querying(
    index: &LazyIndexedTimstofPeaks,
    queries: &[((f32, f32), (f32, f32), (f16, f16))],
) -> QueryResult {
    let start = std::time::Instant::now();
    let mut npeaks = 0;

    for (prec, frag, im) in queries {
        let prec = TupleRange::try_new(prec.0, prec.1).unwrap();
        let frag = TupleRange::try_new(frag.0, frag.1).unwrap();
        let im = TupleRange::try_new(im.0, im.1).unwrap();
        for (_wg, peaks_vec) in index.query_peaks_ms2(prec, frag, Unrestricted, Restricted(im)) {
            npeaks += peaks_vec.len();
        }
    }

    QueryResult {
        npeaks,
        query_time: start.elapsed() / queries.len() as u32,
    }
}

async fn test_querying_async(
    index: &LazyIndexedTimstofPeaks,
    queries: &[((f32, f32), (f32, f32), (f16, f16))],
) -> QueryResult {
    let start = std::time::Instant::now();

    let tasks: Vec<_> = queries
        .iter()
        .map(|(prec, frag, im)| {
            let prec = TupleRange::try_new(prec.0, prec.1).unwrap();
            let frag = TupleRange::try_new(frag.0, frag.1).unwrap();
            let im = TupleRange::try_new(im.0, im.1).unwrap();
            index.query_peaks_ms2_async(prec, frag, Unrestricted, Restricted(im))
        })
        .collect();

    let all_results = futures::future::join_all(tasks).await;
    let npeaks = all_results
        .iter()
        .filter_map(|r| r.as_ref().ok())
        .flat_map(|results| results.iter())
        .map(|(_wg, peaks)| peaks.len())
        .sum();

    QueryResult {
        npeaks,
        query_time: start.elapsed() / queries.len() as u32,
    }
}

fn build_queries(rng: &mut impl rand::Rng) -> Vec<((f32, f32), (f32, f32), (f16, f16))> {
    (0..10)
        .map(|_| {
            let prec_start: f32 = rng.random_range(600.0..800.0);
            let frag_start: f32 = rng.random_range(600.0..800.0);
            let im_start = f16::from_f32(rng.random_range(0.7..1.1));
            (
                (prec_start, prec_start + 0.05),
                (frag_start, frag_start + 0.05),
                (im_start, im_start + f16::from_f32(0.1)),
            )
        })
        .collect()
}
