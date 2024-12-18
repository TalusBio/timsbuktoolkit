use core::f32;
use serde::{
    Deserialize,
    Serialize,
};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MzToleramce {
    Absolute((f64, f64)),
    Ppm((f64, f64)),
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub enum RtTolerance {
    Absolute((f32, f32)),
    Pct((f32, f32)),
    #[default]
    None,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MobilityTolerance {
    Absolute((f32, f32)),
    Pct((f32, f32)),
    None,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum QuadTolerance {
    Absolute((f32, f32)),
}

// TODO: Rename to something that does not use the 'Default'
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DefaultTolerance {
    pub ms: MzToleramce,
    #[serde(default)]
    pub rt: RtTolerance,
    pub mobility: MobilityTolerance,
    pub quad: QuadTolerance,
}

impl Default for DefaultTolerance {
    fn default() -> Self {
        DefaultTolerance {
            ms: MzToleramce::Ppm((20.0, 20.0)),
            rt: RtTolerance::Absolute((5.0, 5.0)),
            mobility: MobilityTolerance::Pct((3.0, 3.0)),
            quad: QuadTolerance::Absolute((0.1, 0.1)),
        }
    }
}

pub trait Tolerance {
    fn mz_range(&self, mz: f64) -> (f64, f64);
    fn rt_range(&self, rt: f32) -> Option<(f32, f32)>;
    fn mobility_range(&self, mobility: f32) -> Option<(f32, f32)>;
    fn quad_range(&self, precursor_mz_range: (f64, f64)) -> (f32, f32);
}

impl Tolerance for DefaultTolerance {
    fn mz_range(&self, mz: f64) -> (f64, f64) {
        match self.ms {
            MzToleramce::Absolute((low, high)) => (mz - low, mz + high),
            MzToleramce::Ppm((low, high)) => {
                let low = mz * low / 1e6;
                let high = mz * high / 1e6;
                (mz - low, mz + high)
            }
        }
    }

    // TODO add an unit ...
    fn rt_range(&self, rt: f32) -> Option<(f32, f32)> {
        match self.rt {
            RtTolerance::Absolute((low, high)) => Some((rt - low, rt + high)),
            RtTolerance::Pct((low, high)) => {
                let low = rt * low / 100.0;
                let high = rt * high / 100.0;
                Some((rt - low, rt + high))
            }
            RtTolerance::None => None,
        }
    }

    fn mobility_range(&self, mobility: f32) -> Option<(f32, f32)> {
        match self.mobility {
            MobilityTolerance::Absolute((low, high)) => Some((mobility - low, mobility + high)),
            MobilityTolerance::Pct((low, high)) => {
                let low = mobility * (low / 100.0);
                let high = mobility * (high / 100.0);
                Some((mobility - low, mobility + high))
            }
            MobilityTolerance::None => None,
        }
    }

    fn quad_range(&self, precursor_mz_range: (f64, f64)) -> (f32, f32) {
        match self.quad {
            QuadTolerance::Absolute((low, high)) => {
                let mz_low = precursor_mz_range.0.min(precursor_mz_range.1) as f32 - low;
                let mz_high = precursor_mz_range.1.max(precursor_mz_range.0) as f32 + high;
                assert!(mz_low <= mz_high);
                assert!(
                    mz_low > 0.0,
                    "Precursor mz is 0 or less, inputs: self: {:?}, precursor_mz_range: {:?}",
                    self,
                    precursor_mz_range,
                );
                (mz_low, mz_high)
            }
        }
    }
}

/// A trait that can be implemented by types that can convert
/// elution groups into queries.
///
/// The elution group here is generic but most regularly it will be
/// an `ElutionGroup` struct.
pub trait ToleranceAdapter<QF, T> {
    fn query_from_elution_group(&self, tol: &dyn Tolerance, elution_group: &T) -> QF;
}
