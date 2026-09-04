use super::TupleRange;

pub fn ppm_tol_range(elem: f64, tol_ppm: f64) -> TupleRange<f64> {
    let utol = elem * (tol_ppm / 1e6);
    let left_e = elem - utol;
    let right_e = elem + utol;
    (left_e, right_e).try_into().unwrap()
}

pub fn pct_tol_range(elem: f64, tol_pct: f64) -> TupleRange<f64> {
    let utol = elem * (tol_pct / 100.0);
    let left_e = elem - utol;
    let right_e = elem + utol;
    (left_e, right_e).try_into().unwrap()
}
