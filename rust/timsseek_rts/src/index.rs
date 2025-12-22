use std::time::Instant;

use timsquery::{
    IndexedTimstofPeaks,
    Tolerance,
};
use timsseek::errors::Result;
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

    let fragmented_range = index.fragmented_range();

    Ok(ScoringPipeline {
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
