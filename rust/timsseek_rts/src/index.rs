use std::collections::HashMap;
use std::sync::Arc;
use std::time::Instant;

use serde::{
    Deserialize,
    Serialize,
};
use timsquery::models::indices::ExpandedRawFrameIndex;
use timsquery::{
    ElutionGroup,
    Tolerance,
};
use timsseek::errors::Result;
use timsseek::models::DigestSlice;
use timsseek::utils::tdf::get_ms1_frame_times_ms;
use timsseek::{
    ExpectedIntensities,
    IonAnnot,
    QueryItemToScore,
    Scorer,
};

// TODO: replace with a trait ... This works for now though
type IndexUse = ExpandedRawFrameIndex;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InputQuery {
    pub sequence: String,
    pub charge: u8,
    pub elution_group: Arc<ElutionGroup<IonAnnot>>,
    pub expected_intensities: ExpectedIntensities,
}

impl InputQuery {
    pub fn sample() -> Self {
        let tmp = QueryItemToScore::sample();
        Self {
            sequence: String::from(tmp.digest),
            charge: tmp.charge,
            elution_group: tmp.query,
            expected_intensities: tmp.expected_intensity,
        }
    }
}

impl From<InputQuery> for QueryItemToScore {
    fn from(value: InputQuery) -> Self {
        Self {
            digest: DigestSlice::from_string(value.sequence.clone(), false),
            charge: value.charge,
            query: value.elution_group,
            expected_intensity: value.expected_intensities,
        }
    }
}

pub fn new_index(
    dotd_file_location: std::path::PathBuf,
    tolerance: Tolerance,
) -> Result<Scorer<IndexUse>> {
    let st = Instant::now();
    // Can use centroided for faster queries ...
    let index = ExpandedRawFrameIndex::from_path(
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
