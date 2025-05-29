use std::time::Instant;

use timsquery::Tolerance;
use timsquery::models::indices::ExpandedRawFrameIndex;
use timsseek::Scorer;
use timsseek::errors::Result;
use timsseek::utils::tdf::get_ms1_frame_times_ms;

// TODO: replace with a trait ... This works for now though
type IndexUse = ExpandedRawFrameIndex;

pub fn new_index(
    dotd_file_location: std::path::PathBuf,
    tolerance: Tolerance,
) -> Result<Scorer<IndexUse>> {
    let st = Instant::now();
    // Can use centroided for faster queries ...
    let index = ExpandedRawFrameIndex::from_path(
        // let index = ExpandedRawFrameIndex::from_path_centroided(
        dotd_file_location
            .clone()
            .to_str()
            .expect("Path is not convertable to string"),
    )?;
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
