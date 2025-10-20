use crate::IonSearchResults;
pub use calibrt::{
    CalibRtError,
    CalibrationCurve as RTCalibration,
    Point,
    calibrate,
};

// It is significantly cheaper to project the spectral library to the observed data
// instead of the other way around. so we need to train the classifier with x = theoretical
// and y = observed.

impl From<IonSearchResults> for Point {
    fn from(val: IonSearchResults) -> Self {
        Point {
            x: val.precursor_rt_query_seconds as f64,
            y: val.obs_rt_seconds as f64,
            weight: val.main_score as f64,
        }
    }
}

// calibrate()
