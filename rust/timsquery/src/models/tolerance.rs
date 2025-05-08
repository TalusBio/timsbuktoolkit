use crate::utils::tolerance_ranges::IncludedRange;
use core::f32;
use serde::{
    Deserialize,
    Serialize,
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
    None,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MobilityTolerance {
    #[serde(rename = "absolute")]
    Absolute((f32, f32)),
    #[serde(rename = "percent")]
    Pct((f32, f32)),
    None,
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
    pub fn mz_range(&self, mz: f64) -> IncludedRange<f64> {
        match self.ms {
            MzTolerance::Absolute((low, high)) => (mz - low, mz + high).into(),
            MzTolerance::Ppm((low, high)) => {
                let low = mz * low / 1e6;
                let high = mz * high / 1e6;
                (mz - low, mz + high).into()
            }
        }
    }

    // TODO add an unit ...
    pub fn rt_range(&self, rt_minutes: f32) -> Option<IncludedRange<f32>> {
        match self.rt {
            RtTolerance::Minutes((low, high)) => Some((rt_minutes - low, rt_minutes + high).into()),
            RtTolerance::Pct((low, high)) => {
                let low = rt_minutes * low / 100.0;
                let high = rt_minutes * high / 100.0;
                Some((rt_minutes - low, rt_minutes + high).into())
            }
            RtTolerance::None => None,
        }
    }

    pub fn rt_range_as_milis(&self, rt_seconds: f32) -> Option<IncludedRange<u32>> {
        let tmp = self.rt_range(rt_seconds);
        tmp.map(|x| ((x.start() * 1000.0) as u32, (x.end() * 1000.0) as u32).into())
    }

    pub fn mobility_range(&self, mobility: f32) -> Option<IncludedRange<f32>> {
        match self.mobility {
            MobilityTolerance::Absolute((low, high)) => {
                Some((mobility - low, mobility + high).into())
            }
            MobilityTolerance::Pct((low, high)) => {
                let low = mobility * (low / 100.0);
                let high = mobility * (high / 100.0);
                Some((mobility - low, mobility + high).into())
            }
            MobilityTolerance::None => None,
        }
    }

    pub fn quad_range(&self, precursor_mz_range: (f64, f64)) -> IncludedRange<f64> {
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
                (mz_low, mz_high).into()
            }
        }
    }

    pub fn indexed_tof_range(&self, mz: f64, converter: &Tof2MzConverter) -> IncludedRange<u32> {
        let mz_rng = self.mz_range(mz);
        (
            converter.invert(mz_rng.start()) as u32,
            converter.invert(mz_rng.end()) as u32,
        )
            .into()
    }

    pub fn indexed_scan_range(
        &self,
        mobility: f64,
        converter: &Scan2ImConverter,
    ) -> Option<IncludedRange<u16>> {
        let im_rng = self.mobility_range(mobility as f32)?;
        Some(
            (
                converter.invert(im_rng.start() as f64) as u16,
                converter.invert(im_rng.end() as f64) as u16,
            )
                .into(),
        )
    }
}
