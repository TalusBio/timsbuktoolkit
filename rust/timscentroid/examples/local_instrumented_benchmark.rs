use OptionallyRestricted::{
    Restricted,
    Unrestricted,
};
use half::f16;
use rand::SeedableRng;
use timscentroid::StorageLocation;
use timscentroid::lazy::LazyIndexedTimstofPeaks;
use timscentroid::utils::{
    OptionallyRestricted,
    TupleRange,
};

fn main() {
    let serialized_dir = std::path::Path::new("./serialized_peaks");

    // Check if serialized data exists locally
    if !serialized_dir.exists() {
        eprintln!("Error: Serialized data not found at {:?}", serialized_dir);
        eprintln!("Please run `cargo run --example serialization --release` first");
        std::process::exit(1);
    }

    println!("=== Local Instrumented Benchmark ===\n");

    // Configuration
    const NUM_QUERIES: usize = 100; // Increased from 10 now that we're fast!
    const RANDOM_SEED: u64 = 42; // Deterministic seed for reproducibility

    // Load index with instrumentation
    println!("Loading lazy index...");
    let location = StorageLocation::from_path(serialized_dir);
    let index = LazyIndexedTimstofPeaks::load_from_storage(location)
        .unwrap()
        .with_instrumentation("Local Benchmark");

    // Generate deterministic test queries
    let queries = {
        let mut rng = rand::rngs::StdRng::seed_from_u64(RANDOM_SEED);
        build_queries(&mut rng, NUM_QUERIES)
    };

    println!("Running {} queries (Sync API)...\n", queries.len());

    // Run sync queries
    let start = std::time::Instant::now();
    let mut total_peaks = 0;
    let mut total_intensity = 0.0;

    for (i, (prec, frag, im)) in queries.iter().enumerate() {
        let prec = TupleRange::try_new(prec.0, prec.1).unwrap();
        let frag = TupleRange::try_new(frag.0, frag.1).unwrap();
        let im = TupleRange::try_new(im.0, im.1).unwrap();

        let results = index.query_peaks_ms2(prec, frag, Unrestricted, Restricted(im));

        for (_wg, peaks_vec) in results {
            for peak in peaks_vec {
                total_intensity += peak.intensity as f64;
                total_peaks += 1;
            }
        }

        if i == 0 {
            // Print metrics after first query to see initialization overhead
            println!("=== After First Query ===");
            index.print_metrics("First Query");
            println!();
        }
    }

    let elapsed = start.elapsed();
    println!("\n=== Final Results ===");
    println!("Total time: {:?}", elapsed);
    println!("Time per query: {:?}", elapsed / queries.len() as u32);
    println!("Total peaks: {}", total_peaks);
    println!(
        "Peaks per query: {:.2}",
        total_peaks as f64 / queries.len() as f64
    );
    println!("Total intensity: {:.2e}", total_intensity);

    println!("\n=== Final Metrics ===");
    index.print_metrics("All Queries");

    // Now test async API
    println!("\n\n=== Testing Async API ===\n");

    let location = StorageLocation::from_path(serialized_dir);
    let index_async = LazyIndexedTimstofPeaks::load_from_storage(location)
        .unwrap()
        .with_instrumentation("Local Async Benchmark");

    // Regenerate same queries for async test
    let queries_async = {
        let mut rng = rand::rngs::StdRng::seed_from_u64(RANDOM_SEED);
        build_queries(&mut rng, NUM_QUERIES)
    };

    let runtime = tokio::runtime::Runtime::new().unwrap();

    println!(
        "Running {} queries concurrently (Async API)...\n",
        queries_async.len()
    );

    let start = std::time::Instant::now();
    let (total_peaks_async, total_intensity_async) = runtime.block_on(async {
        let mut tasks = Vec::new();

        for (prec, frag, im) in &queries_async {
            let prec = TupleRange::try_new(prec.0, prec.1).unwrap();
            let frag = TupleRange::try_new(frag.0, frag.1).unwrap();
            let im = TupleRange::try_new(im.0, im.1).unwrap();

            tasks.push(index_async.query_peaks_ms2_async(prec, frag, Unrestricted, Restricted(im)));
        }

        let all_results = futures::future::join_all(tasks).await;

        let mut total_peaks = 0;
        let mut total_intensity = 0.0;

        for result in all_results {
            match result {
                Ok(results) => {
                    for (_wg, peaks_vec) in results {
                        for peak in peaks_vec {
                            total_intensity += peak.intensity as f64;
                            total_peaks += 1;
                        }
                    }
                }
                Err(e) => {
                    eprintln!("Query error: {}", e);
                }
            }
        }

        (total_peaks, total_intensity)
    });

    let elapsed_async = start.elapsed();

    println!("\n=== Async API Results ===");
    println!("Total time: {:?}", elapsed_async);
    println!("Time per query: {:?}", elapsed_async / queries.len() as u32);
    println!("Total peaks: {}", total_peaks_async);
    println!(
        "Peaks per query: {:.2}",
        total_peaks_async as f64 / queries.len() as f64
    );
    println!("Total intensity: {:.2e}", total_intensity_async);

    println!("\n=== Async API Metrics ===");
    index_async.print_metrics("All Async Queries");

    // Compare
    println!("\n=== Performance Comparison ===");
    let speedup = elapsed.as_secs_f64() / elapsed_async.as_secs_f64();
    println!("Async is {:.2}x faster than sync", speedup);

    // Verify correctness
    if total_peaks == total_peaks_async && (total_intensity - total_intensity_async).abs() < 0.001 {
        println!("✓ Sync and async results match!");
    } else {
        println!("✗ Results differ!");
        println!(
            "  Sync:  {} peaks, {:.2e} intensity",
            total_peaks, total_intensity
        );
        println!(
            "  Async: {} peaks, {:.2e} intensity",
            total_peaks_async, total_intensity_async
        );
    }
}

fn build_queries(
    rng: &mut impl rand::Rng,
    num_queries: usize,
) -> Vec<((f32, f32), (f32, f32), (f16, f16))> {
    let mut out = Vec::with_capacity(num_queries);
    for _ in 0..num_queries {
        let prec_start: f32 = rng.random_range(600.0..800.0);
        let prec_end = prec_start + 0.05;
        let prec = (prec_start, prec_end).try_into().unwrap();

        let frag_start: f32 = rng.random_range(600.0..800.0);
        let frag_end = frag_start + 0.05;
        let frag = (frag_start, frag_end).try_into().unwrap();

        let im_start_i: f32 = rng.random_range(0.7..1.1);
        let im_start: f16 = f16::from_f32(im_start_i);
        let im_end = im_start + (f16::from_f32(0.1f32));
        let im = (im_start, im_end).try_into().unwrap();

        out.push((prec, frag, im));
    }
    out
}
