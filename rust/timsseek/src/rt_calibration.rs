use crate::{
    IonSearchResults,
    Speclib,
};
pub use calibrt::{
    CalibRtError,
    CalibrationCurve as RTCalibration,
    Point,
    calibrate,
};
use tracing::warn;

// It is significantly cheaper to project the spectral library to the observed data
// instead of the other way around. so we need to train the classifier with x = theoretical
// and y = observed.

impl From<&IonSearchResults> for Point {
    fn from(val: &IonSearchResults) -> Self {
        Point {
            x: val.precursor_rt_query_seconds as f64,
            y: val.obs_rt_seconds as f64,
            weight: 1.0,
        }
        // weight: val.main_score.ln_1p().max(0.0) as f64, // very non-good
        // so just weight equally for now.
        // weight: val.main_score as f64, // works ok...
        // In theory we dont want to use target-decoy info here, to prevent biasing.
        // weight: if val.qvalue < 0.1 { 1.0 } else { 0.0 },
        // weight: 1.0,
    }
}

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
pub fn recalibrate_speclib(
    speclib: &mut Speclib,
    calib_data: &[IonSearchResults],
) -> Result<RTCalibration, CalibRtError> {
    let mut min_x = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_y = f64::NEG_INFINITY;

    calib_data.iter().for_each(|spec| {
        let y = spec.obs_rt_seconds as f64;
        if y < min_y {
            min_y = y;
        }
        if y > max_y {
            max_y = y;
        }
        let x = spec.precursor_rt_query_seconds as f64;
        if x < min_x {
            min_x = x;
        }
        if x > max_x {
            max_x = x;
        }
    });

    let cal_res = calibrate(
        calib_data
            .iter()
            .map(Point::from)
            .collect::<Vec<Point>>()
            .as_slice(),
        (min_x, max_x),
        (min_y, max_y),
        50,
    );

    let mut oob_preds = (0, f32::MAX, f32::MIN);
    let mut cool_preds = 0;
    match cal_res {
        Ok(cal_curve) => {
            speclib.elems.iter_mut().for_each(|spec| {
                let pred = cal_curve.predict(spec.query.rt_seconds as f64);
                let pred = match pred {
                    Ok(pred_rt) => {
                        // I can get this by total - oob but this feels safer ... even if its a
                        // hair slower.
                        cool_preds += 1;
                        pred_rt as f32
                    }
                    Err(CalibRtError::OutOfBounds(pred_rt)) => {
                        oob_preds.0 += 1;
                        let query_rt = spec.query.rt_seconds;
                        oob_preds.1 = oob_preds.1.min(query_rt);
                        oob_preds.2 = oob_preds.2.max(query_rt);

                        pred_rt as f32
                    }
                    Err(_) => {
                        panic!("Unexpected error during RT prediction");
                    }
                };
                spec.query.rt_seconds = pred;
            });
            if oob_preds.0 > 0 {
                warn!(
                    "{}/{} out of bounds RT predictions (min RT {}s, max RT {}s)",
                    oob_preds.0,
                    cool_preds + oob_preds.0,
                    oob_preds.1,
                    oob_preds.2
                );
            }

            Ok(cal_curve)
        }

        Err(e) => Err(e),
    }
}

pub fn recalibrate_results(calibration: &RTCalibration, results: &mut [IonSearchResults]) {
    for v in results.iter_mut() {
        let pred_rt = calibration
            .predict(v.precursor_rt_query_seconds as f64)
            .unwrap_or_else(|e| match e {
                CalibRtError::OutOfBounds(x) => x,
                _ => panic!("Unexpected error during RT prediction"),
            });
        v.recalibrated_query_rt = pred_rt as f32;
        v.calibrated_sq_delta_theo_rt = (v.obs_rt_seconds - v.recalibrated_query_rt).powi(2);
        v.delta_theo_rt = v.obs_rt_seconds - v.recalibrated_query_rt;
    }
}
