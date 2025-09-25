use std::time::Instant;

use timsquery::{
    CentroidingConfig,
    IndexedTimstofPeaks,
    TimsTofPath,
    Tolerance,
};
use timsseek::Scorer;
use timsseek::errors::Result;
use timsseek::utils::tdf::get_ms1_frame_times_ms;

pub fn new_index(
    dotd_file_location: std::path::PathBuf,
    tolerance: Tolerance,
) -> Result<Scorer<IndexedTimstofPeaks>> {
    let st = Instant::now();
    // Can use centroided for faster queries ...
    //
    let file = TimsTofPath::new(
        dotd_file_location
            .clone()
            .to_str()
            .expect("Path is not convertable to string"),
    )?;
    // TODO: make an actual error...
    let centroiding_config = CentroidingConfig {
        max_peaks: 50_000,
        mz_ppm_tol: 10.0,
        im_pct_tol: 5.0,
        early_stop_iterations: 200,
    };
    print!("Loading index... (this might take a min) ");
    let (index, _stats) = IndexedTimstofPeaks::from_timstof_file(&file, centroiding_config);
    let elap_time = st.elapsed();
    println!(
        "Loading index took: {:?} for {}",
        elap_time,
        dotd_file_location.display()
    );

    let tdf_path = &dotd_file_location.clone().join("analysis.tdf");
    let ref_time_ms = get_ms1_frame_times_ms(tdf_path.to_str().unwrap()).unwrap();
    let fragmented_range = index.fragmented_range();

    Ok(Scorer {
        index_cycle_rt_ms: ref_time_ms,
        index,
        tolerance: tolerance.clone(),
        secondary_tolerance: tolerance.with_rt_tolerance(
            timsquery::models::tolerance::RtTolerance::Minutes((0.5, 0.5)),
        ),
        fragmented_range,
    })
}
