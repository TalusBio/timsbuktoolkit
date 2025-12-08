use std::time::Instant;

use timsquery::{
    IndexedTimstofPeaks,
    Tolerance,
};
use timsseek::errors::Result;
use timsseek::utils::tdf::get_ms1_frame_times_ms;
use timsseek::{
    ScoringPipeline,
    ToleranceHierarchy,
};

pub fn new_index(
    dotd_file_location: std::path::PathBuf,
    tolerance: Tolerance,
) -> Result<ScoringPipeline<IndexedTimstofPeaks>> {
    let st = Instant::now();
    // Can use centroided for faster queries ...
    //
    let file_loc = dotd_file_location.clone();
    let index = timsquery::serde::load_index_caching(file_loc).unwrap();
    let elap_time = st.elapsed();
    println!(
        "Loading index took: {:?} for {}",
        elap_time,
        dotd_file_location.display()
    );

    let tdf_path = &dotd_file_location.clone().join("analysis.tdf");
    let ref_time_ms = get_ms1_frame_times_ms(tdf_path.to_str().unwrap()).unwrap();
    let fragmented_range = index.fragmented_range();

    Ok(ScoringPipeline {
        index_cycle_rt_ms: ref_time_ms,
        index,
        tolerances: ToleranceHierarchy {
            prescore: tolerance.clone(),
            secondary: tolerance.with_rt_tolerance(
                timsquery::models::tolerance::RtTolerance::Minutes((0.1, 0.1)),
            ),
        },
        fragmented_range,
    })
}
