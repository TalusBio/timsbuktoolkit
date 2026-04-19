//! Real-speclib per-peak iteration bench.
//!
//! Runs `add_query` against the production index + speclib for two
//! aggregators that touch every matching peak with different sink
//! costs:
//!   - `PointIntensityAggregator`: `agg.intensity += peak.intensity as f64`
//!   - `ChromatogramCollector`: scatter-write into cycle-indexed rows
//!
//! Both go through the same `for_each_peak` / `for_each_ms1_peak` scan
//! machinery, so time-delta between them attributes to sink cost. Total
//! time per aggregator attributes the hot path (scan + filter + sink)
//! at realistic peak density — library fragment m/zs cluster around
//! real peptide peaks, unlike a uniform-random mz bench.
//!
//! Run: cargo run -r -p timsseek --example query_bench
//!
//! Env overrides:
//!   BENCH_DOTD    — .d path (default Hela)
//!   BENCH_SPECLIB — speclib path (default asdad)
//!   QB_N         — number of speclib entries to process (default 2000)
//!   QB_ITERS     — outer repeat count (default 1)

use std::path::Path;
use std::time::Instant;
use timsquery::models::aggregators::{
    ChromatogramCollector,
    PointIntensityAggregator,
};
use timsquery::models::tolerance::Tolerance;
use timsquery::serde::{
    IndexedPeaksHandle,
    load_index_auto,
};
use timsquery::traits::queriable_data::QueriableData;
use timsquery::utils::TupleRange;
use timsseek::IonAnnot;
use timsseek::data_sources::speclib::Speclib;
use timsseek::models::DecoyStrategy;

fn env(key: &str, default: &str) -> String {
    std::env::var(key).unwrap_or_else(|_| default.to_string())
}

fn main() {
    let dotd = env(
        "BENCH_DOTD",
        "/Users/sebastianpaez/data/decompressed_timstof/250225_Desnaux_200ng_Hela_ICC_off_DIA.d",
    );
    let speclib_path = env(
        "BENCH_SPECLIB",
        "/Users/sebastianpaez/fasta/asdad.msgpack.zstd",
    );
    let n: usize = env("QB_N", "2000").parse().unwrap();
    let iters: usize = env("QB_ITERS", "1").parse().unwrap();

    eprintln!("Loading index from {dotd}...");
    let t0 = Instant::now();
    let handle = load_index_auto(&dotd, None).expect("load index");
    let index = match handle {
        IndexedPeaksHandle::Eager(eager) => eager,
        IndexedPeaksHandle::Lazy(_) => panic!("lazy index not supported by this bench"),
    };
    eprintln!("  index ready in {:?}", t0.elapsed());

    eprintln!("Loading speclib from {speclib_path}...");
    let t0 = Instant::now();
    let speclib = Speclib::from_file(Path::new(&speclib_path), DecoyStrategy::default())
        .expect("load speclib");
    eprintln!(
        "  speclib ready in {:?} ({} entries)",
        t0.elapsed(),
        speclib.len()
    );

    let tolerance = Tolerance::default();
    eprintln!("Tolerance: {:?}", tolerance);

    // Cycle-mapping range — use the index's full range so no RT filtering
    // prunes peaks away; isolates per-peak iteration cost.
    let cycle_mapping = index.ms1_cycle_mapping();
    let (min_rt, max_rt) = cycle_mapping.range_milis();
    let rt_range = TupleRange::try_new(min_rt, max_rt).expect("rt range valid");

    let items: Vec<_> = speclib.as_slice().iter().take(n).collect();
    eprintln!("Benching {} items × {} iters each", items.len(), iters);

    // ---- Point aggregator ----
    {
        let t0 = Instant::now();
        let mut total_intensity = 0.0f64;
        for _ in 0..iters {
            for item in &items {
                let mut agg = PointIntensityAggregator::<IonAnnot>::new(&item.query);
                index.add_query(&mut agg, &tolerance);
                total_intensity += agg.intensity;
            }
        }
        let e = t0.elapsed();
        let pep_total = items.len() * iters;
        eprintln!(
            "POINT   total={:?} per-pep={:.1}µs total_intensity={:.3e}",
            e,
            e.as_nanos() as f64 / pep_total as f64 / 1000.0,
            total_intensity,
        );
    }

    // ---- Chromatogram aggregator ----
    {
        let t0 = Instant::now();
        let mut total_frag_peaks = 0u64;
        for _ in 0..iters {
            for item in &items {
                let Ok(mut agg) = ChromatogramCollector::<IonAnnot, f32>::new(
                    &item.query,
                    rt_range,
                    cycle_mapping,
                ) else {
                    continue;
                };
                index.add_query(&mut agg, &tolerance);
                total_frag_peaks += agg.n_fragment_peaks_added;
            }
        }
        let e = t0.elapsed();
        let pep_total = items.len() * iters;
        eprintln!(
            "CHROM   total={:?} per-pep={:.1}µs total_ms2_peaks={}",
            e,
            e.as_nanos() as f64 / pep_total as f64 / 1000.0,
            total_frag_peaks,
        );
        eprintln!(
            "  (ns/peak ≈ {:.1})",
            e.as_nanos() as f64 / total_frag_peaks as f64,
        );
    }

    // No-op when built without `query-instr`; prints filter-funnel
    // counts when built `--features timscentroid/query-instr`.
    timscentroid::indexing::dump_for_each_peak_funnel("query_bench");
}
