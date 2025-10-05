use half::f16;
use timscentroid::utils::OptionallyRestricted;
use timscentroid::{
    CentroidingConfig,
    IndexedTimstofPeaks,
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
    println!("Indexing with config: {:#?}", centroiding_config);
    let (index, index_stats) = IndexedTimstofPeaks::from_timstof_file(&file, centroiding_config);
    println!("Indexing Stats: {}", index_stats);

    // test_querying(&index);
}

fn test_querying(index: &IndexedTimstofPeaks) {
    let mut tot_int = 0.0;
    let mut nqueries = 0;
    let mut npeaks = 0;
    let st = std::time::Instant::now();
    for prec_mz_start_i in 600..700 {
        let prec_start: f32 = prec_mz_start_i as f32;
        let prec_end = prec_start + 0.05;
        let prec = (prec_start, prec_end).try_into().unwrap();
        for frag_mz_start_i in 600..700 {
            let frag_start: f32 = frag_mz_start_i as f32;
            let frag_end = frag_start + 0.05;
            let frag = (frag_start, frag_end).try_into().unwrap();
            for im_start_i in 500..1500 {
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

    const PEPTIDES_PER_PROTEOME: usize = 1_854_958;
    const TRANSITIONS_PER_PEPTIDE: usize = 10;

    let total_transitions = PEPTIDES_PER_PROTEOME * TRANSITIONS_PER_PEPTIDE;
    let time_to_proteome = query_time * total_transitions as u32;
    let peaks_per_query = npeaks as f64 / nqueries as f64;

    println!(
        "Time querying: {:?} {} queries, time per query: {:?}",
        et, nqueries, query_time
    );
    println!(
        "Total intensity: {} on {} peaks; peaks per query: {}",
        tot_int, npeaks, peaks_per_query
    );
    println!(
        "Estimated time to query a full proteome ({} transitions): {:?}",
        total_transitions, time_to_proteome
    );
}
