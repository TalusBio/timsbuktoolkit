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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MzToleramce {
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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Tolerance {
    pub ms: MzToleramce,
    #[serde(default)]
    pub rt: RtTolerance,
    pub mobility: MobilityTolerance,
    pub quad: QuadTolerance,
}

impl Default for Tolerance {
    fn default() -> Self {
        Tolerance {
            ms: MzToleramce::Ppm((20.0, 20.0)),
            rt: RtTolerance::Minutes((5.0, 5.0)),
            mobility: MobilityTolerance::Pct((3.0, 3.0)),
            quad: QuadTolerance::Absolute((0.1, 0.1)),
        }
    }
}

impl Tolerance {
    pub fn mz_range(&self, mz: f64) -> IncludedRange<f64> {
        match self.ms {
            MzToleramce::Absolute((low, high)) => (mz - low, mz + high).into(),
            MzToleramce::Ppm((low, high)) => {
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

    pub fn rt_range_as_milis(&self, rt_minutes: f32) -> Option<IncludedRange<u32>> {
        let tmp = self.rt_range(rt_minutes);
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

// impl NaturalPrecursorQuery {
//     pub fn as_precursor_query(
//         &self,
//         mz_converter: &Tof2MzConverter,
//         im_converter: &Scan2ImConverter,
//     ) -> PrecursorIndexQuery {
//         PrecursorIndexQuery {
//             frame_index_range: (
//                 im_converter.invert(self.mobility_range.start()).round() as usize,
//                 im_converter.invert(self.mobility_range.end()).round() as usize,
//             )
//                 .into(),
//             rt_range_seconds: self.rt_range,
//             mz_index_ranges: self
//                 .mz_ranges
//                 .iter()
//                 .map(|mz_range| {
//                     (
//                         mz_converter.invert(mz_range.start()).round() as u32,
//                         mz_converter.invert(mz_range.end()).round() as u32,
//                     )
//                         .into()
//                 })
//                 .collect(),
//             mobility_index_range: (
//                 im_converter.invert(self.mobility_range.start()).round() as usize,
//                 im_converter.invert(self.mobility_range.end()).round() as usize,
//             )
//                 .into(),
//             isolation_mz_range: self.isolation_mz_range,
//         }
//     }
// }
