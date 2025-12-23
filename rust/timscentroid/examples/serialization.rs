use half::f16;
use timscentroid::utils::OptionallyRestricted;
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

fn main() {
    // Set up logger
    tracing_subscriber::fmt()
        .with_span_events(tracing_subscriber::fmt::format::FmtSpan::CLOSE)
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .init();

    const DIA_TEST: &str =
        "/Users/sebastianpaez/data/decompressed_timstof/250225_Desnaux_200ng_Hela_ICC_on_DIA.d/";
    let file = TimsTofPath::new(DIA_TEST).unwrap();

    let centroiding_config = CentroidingConfig {
        max_peaks: 20_000,
        mz_ppm_tol: 5.0,
        im_pct_tol: 3.0,
        early_stop_iterations: 200,
    };

    println!("=== Centroiding & Indexing ===");
    println!("Config: {:#?}", centroiding_config);

    let start_centroid = std::time::Instant::now();
    let (index_original, index_stats) =
        IndexedTimstofPeaks::from_timstof_file(&file, centroiding_config);
    let centroid_time = start_centroid.elapsed();

    println!("Indexing Stats: {}", index_stats);
    println!("Centroiding time: {:?}", centroid_time);

    // Serialize to disk
    println!("\n=== Serialization ===");
    let output_dir = std::path::Path::new("./serialized_peaks");

    // Clean up previous serialization if it exists
    if output_dir.exists() {
        std::fs::remove_dir_all(output_dir).unwrap();
    }

    let start_write = std::time::Instant::now();
    index_original.save_to_directory(output_dir).unwrap();
    let write_time = start_write.elapsed();

    println!("Serialization time: {:?}", write_time);

    // Calculate disk usage
    let disk_usage = calculate_directory_size(output_dir).unwrap();
    println!(
        "Disk usage: {:.2} MB ({} bytes)",
        disk_usage as f64 / 1_048_576.0,
        disk_usage
    );

    // Deserialize from disk
    println!("\n=== Deserialization ===");
    let start_read = std::time::Instant::now();
    let location = StorageLocation::from_path(output_dir);
    let index_loaded = IndexedTimstofPeaks::load_from_storage(location).unwrap();
    let read_time = start_read.elapsed();

    println!("Deserialization time: {:?}", read_time);

    // Compare performance
    println!("\n=== Performance Comparison ===");
    println!("Centroiding:      {:>10.2?}", centroid_time);
    println!("Serialization:    {:>10.2?}", write_time);
    println!("Deserialization:  {:>10.2?}", read_time);
    println!(
        "Speedup (load vs centroid): {:.1}x faster",
        centroid_time.as_secs_f64() / read_time.as_secs_f64()
    );
    compare_results(&index_original, &index_loaded);
}

fn compare_results(index_original: &IndexedTimstofPeaks, index_loaded: &IndexedTimstofPeaks) {
    // Test querying on both to ensure functional equivalence
    println!("\n=== Query Test (Original) ===");
    let original_result = test_querying(&index_original);

    println!("\n=== Query Test (Loaded) ===");
    let loaded_result = test_querying(&index_loaded);

    println!("\n=== Query Result Comparison ===");
    println!(
        "Original: {} peaks, {:.2e} total intensity",
        original_result.npeaks, original_result.tot_int
    );
    println!(
        "Loaded:   {} peaks, {:.2e} total intensity",
        loaded_result.npeaks, loaded_result.tot_int
    );

    let peak_diff = (original_result.npeaks as i64 - loaded_result.npeaks as i64).abs();
    let int_diff = (original_result.tot_int - loaded_result.tot_int).abs();
    let int_rel_diff = int_diff / original_result.tot_int * 100.0;

    println!("Peak count difference: {}", peak_diff);
    println!(
        "Intensity difference: {:.2e} ({:.4}%)",
        int_diff, int_rel_diff
    );

    if peak_diff == 0 && int_rel_diff < 0.001 {
        println!("✓ Query results are equivalent!");
    } else {
        println!("⚠ Query results differ slightly (expected for f16 precision)");
    }
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

#[derive(Debug)]
struct QueryResult {
    npeaks: usize,
    tot_int: f64,
}

fn test_querying(index: &IndexedTimstofPeaks) -> QueryResult {
    let mut tot_int = 0.0;
    let mut nqueries = 0;
    let mut npeaks = 0;
    let st = std::time::Instant::now();

    // Reduced query space for faster testing
    for prec_mz_start_i in 600..650 {
        let prec_start: f32 = prec_mz_start_i as f32;
        let prec_end = prec_start + 0.05;
        let prec = (prec_start, prec_end).try_into().unwrap();

        for frag_mz_start_i in 600..650 {
            let frag_start: f32 = frag_mz_start_i as f32;
            let frag_end = frag_start + 0.05;
            let frag = (frag_start, frag_end).try_into().unwrap();

            for im_start_i in 700..1200 {
                let im_start: f16 = f16::from_f32(im_start_i as f32 / 1000.0);
                let im_end = im_start + (f16::from_f32(0.05f32));
                let im = (im_start, im_end).try_into().unwrap();

                nqueries += 1;
                index
                    .query_peaks_ms2(prec, frag, Unrestricted, Restricted(im))
                    .for_each(|(_wg, iiter)| {
                        iiter.for_each(|pp| {
                            tot_int += pp.intensity as f64;
                            npeaks += 1;
                        });
                    });
            }
        }
    }

    let et = st.elapsed();
    let query_time = et / nqueries;
    let peaks_per_query = npeaks as f64 / nqueries as f64;

    println!(
        "Time querying: {:?} {} queries, time per query: {:?}",
        et, nqueries, query_time
    );
    println!(
        "Total intensity: {:.2e} on {} peaks; peaks per query: {:.2}",
        tot_int, npeaks, peaks_per_query
    );

    QueryResult { npeaks, tot_int }
}
