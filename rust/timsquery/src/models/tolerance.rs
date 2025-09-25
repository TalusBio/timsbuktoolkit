use OptionallyRestricted::{
    Restricted,
    Unrestricted,
};
use core::f32;
use half::f16;
use serde::{
    Deserialize,
    Serialize,
};
use timscentroid::utils::{
    OptionallyRestricted,
    TupleRange,
};
use timsrust::converters::{
    ConvertableDomain,
    Scan2ImConverter,
    Tof2MzConverter,
};

/// Tolerance settings for the search.
///
/// This is meant to encapsulate the different tolerances needed
/// for all the dimensions in a search.
///
/// Example:
/// ```
/// use timsquery::Tolerance;
///
/// let tolerance = Tolerance::default();
/// ```
///
/// Since every dimenion has a different way of defining its tolerance
/// every dimension has its own type (highly encourage you to thech the options
/// for each dimension).
///
/// Convention:
/// In contrast with how some software defines tolerance, here we define the ranges
/// in terms of positive values. For instance, here a tolerance of (1,1) on a value
/// of 10 means a range of (9,11) while in some software the same range would be defined
/// as (-1,1).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Tolerance {
    pub ms: MzTolerance,
    #[serde(default)]
    pub rt: RtTolerance,
    pub mobility: MobilityTolerance,
    pub quad: QuadTolerance,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MzTolerance {
    #[serde(rename = "da")]
    Absolute((f64, f64)),
    #[serde(rename = "ppm")]
    Ppm((f64, f64)),
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub enum RtTolerance {
    #[serde(rename = "minutes")]
    Minutes((f32, f32)),
    #[serde(rename = "percent")]
    Pct((f32, f32)),
    #[default]
    Unrestricted,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MobilityTolerance {
    #[serde(rename = "absolute")]
    Absolute((f32, f32)),
    #[serde(rename = "percent")]
    Pct((f32, f32)),
    Unrestricted,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum QuadTolerance {
    #[serde(rename = "absolute")]
    Absolute((f32, f32)),
}

impl Default for Tolerance {
    fn default() -> Self {
        Tolerance {
            ms: MzTolerance::Ppm((20.0, 20.0)),
            rt: RtTolerance::Minutes((5.0, 5.0)),
            mobility: MobilityTolerance::Pct((3.0, 3.0)),
            quad: QuadTolerance::Absolute((0.1, 0.1)),
        }
    }
}

impl Tolerance {
    pub fn mz_range(&self, mz: f64) -> TupleRange<f64> {
        match self.ms {
            MzTolerance::Absolute((low, high)) => (mz - low, mz + high).try_into().expect(
                "mz tolerance should never result in an invalid range, since low and high are positive",
            ),
            MzTolerance::Ppm((low, high)) => {
                let low = mz * low / 1e6;
                let high = mz * high / 1e6;
                (mz - low, mz + high).try_into().expect(
                    "mz tolerance should never result in an invalid range, since low and high are positive",
                )
            }
        }
    }

    pub fn mz_range_f32(&self, mz: f32) -> TupleRange<f32> {
        let tmp = self.mz_range(mz as f64);
        (tmp.start() as f32, tmp.end() as f32).try_into().unwrap()
    }

    pub fn rt_range_seconds_f16(&self, rt_seconds: f32) -> OptionallyRestricted<TupleRange<f16>> {
        let minutes = rt_seconds / 60.0;
        let tmp = self.rt_range_minutes(minutes);
        let sixty = f16::from_f32(60.0);
        // TODO: make this significantly more efficient

        match tmp {
            Restricted(x) => Restricted(
                (
                    f16::from_f32(x.start()) * sixty,
                    f16::from_f32(x.end()) * sixty,
                )
                    .try_into()
                    .unwrap(),
            ),
            Unrestricted => Unrestricted,
        }
    }

    pub fn rt_range_minutes(&self, rt_minutes: f32) -> OptionallyRestricted<TupleRange<f32>> {
        match self.rt {
            RtTolerance::Minutes((low, high)) => {
                Restricted((rt_minutes - low, rt_minutes + high).try_into().unwrap())
            }
            RtTolerance::Pct((low, high)) => {
                let low = rt_minutes * low / 100.0;
                let high = rt_minutes * high / 100.0;
                Restricted((rt_minutes - low, rt_minutes + high).try_into().unwrap())
            }
            RtTolerance::Unrestricted => Unrestricted,
        }
    }

    pub fn rt_range_as_milis(&self, rt_seconds: f32) -> OptionallyRestricted<TupleRange<u32>> {
        let minutes = rt_seconds / 60.0;
        let tmp = self.rt_range_minutes(minutes);
        match tmp {
            Restricted(x) => {
                let start_seconds = x.start() * 60.0;
                let end_seconds = x.end() * 60.0;
                Restricted(
                    (
                        (start_seconds * 1000.0) as u32,
                        (end_seconds * 1000.0) as u32,
                    )
                        .try_into()
                        .unwrap(),
                )
            }
            Unrestricted => Unrestricted,
        }
    }

    pub fn mobility_range(&self, mobility: f32) -> OptionallyRestricted<TupleRange<f32>> {
        match self.mobility {
            MobilityTolerance::Absolute((low, high)) => {
                Restricted((mobility - low, mobility + high).try_into().unwrap())
            }
            MobilityTolerance::Pct((low, high)) => {
                let low = mobility * (low / 100.0);
                let high = mobility * (high / 100.0);
                Restricted((mobility - low, mobility + high).try_into().unwrap())
            }
            MobilityTolerance::Unrestricted => Unrestricted,
        }
    }

    pub fn mobility_range_f16(&self, mobility: f32) -> OptionallyRestricted<TupleRange<f16>> {
        let tmp = self.mobility_range(mobility);
        match tmp {
            Restricted(x) => Restricted(
                (f16::from_f32(x.start()), f16::from_f32(x.end()))
                    .try_into()
                    .unwrap(),
            ),
            Unrestricted => Unrestricted,
        }
    }

    pub fn quad_range_f32(&self, precursor_mz_range: (f32, f32)) -> TupleRange<f32> {
        let tmp = self.quad_range((precursor_mz_range.0 as f64, precursor_mz_range.1 as f64));
        (tmp.start() as f32, tmp.end() as f32).try_into().unwrap()
    }

    pub fn quad_range(&self, precursor_mz_range: (f64, f64)) -> TupleRange<f64> {
        match self.quad {
            QuadTolerance::Absolute((low, high)) => {
                let mz_low = precursor_mz_range.0.min(precursor_mz_range.1) - (low as f64);
                let mz_high = precursor_mz_range.1.max(precursor_mz_range.0) + (high as f64);
                assert!(mz_low <= mz_high);
                assert!(
                    mz_low > 0.0,
                    "Precursor mz is 0 or less, inputs: self: {:?}, precursor_mz_range: {:?}",
                    self,
                    precursor_mz_range,
                );
                (mz_low, mz_high).try_into().unwrap()
            }
        }
    }

    pub fn indexed_tof_range(&self, mz: f64, converter: &Tof2MzConverter) -> TupleRange<u32> {
        let mz_rng = self.mz_range(mz);
        (
            converter.invert(mz_rng.start()) as u32,
            converter.invert(mz_rng.end()) as u32,
        )
            .try_into()
            .unwrap()
    }

    pub fn indexed_scan_range(
        &self,
        mobility: f64,
        converter: &Scan2ImConverter,
    ) -> OptionallyRestricted<TupleRange<u16>> {
        match self.mobility_range(mobility as f32) {
            Restricted(im_rng) => {
                let im_rng = (im_rng.start() as f64, im_rng.end() as f64);
                Restricted(
                    (
                        converter.invert(im_rng.0) as u16,
                        converter.invert(im_rng.1) as u16,
                    )
                        .try_into()
                        .unwrap(),
                )
            }
            Unrestricted => Unrestricted,
        }
    }

    pub fn with_rt_tolerance(&self, rt: RtTolerance) -> Self {
        Self { rt, ..self.clone() }
    }
}
