/// Stirling's approximation for log factorial
// Code shamelessly stolen from
// https://github.com/lazear/sage/blob/445f5ab156c24f5c8f21098717692077e3b1d1ee/crates/sage/src/scoring.rs#L149C1-L157C2
//
pub fn lnfact(n: u16) -> f64 {
    if n == 0 {
        1.0
    } else {
        let n = n as f64;
        n * n.ln() - n + 0.5 * n.ln() + 0.5 * (std::f64::consts::PI * 2.0 * n).ln()
    }
}

pub fn lnfact_float(n: f64) -> f64 {
    if n < 1.0 {
        0.0
    } else {
        n * n.ln() - n + 0.5 * n.ln() + 0.5 * (std::f64::consts::PI * 2.0 * n).ln()
    }
}

/// Logarigthmic mean of a slice of values.
pub fn lnmean(vals: &[f64]) -> f64 {
    let mut sum = 0.0;
    for val in vals {
        sum += val.ln();
    }
    (sum / vals.len() as f64).exp()
}
