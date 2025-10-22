use timsquery::{
    MzMobilityStatsCollector,
    SpectralCollector,
};

use crate::IonAnnot;
use crate::utils::top_n_array::TopNArray;

#[derive(Debug)]
pub struct MzMobilityOffsets {
    pub(crate) ms1: TopNArray<3, SortableError>,
    pub(crate) ms2: TopNArray<7, SortableError>,
    pub(crate) ref_mobility: f64,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct SortableError {
    pub(crate) intensity: f64,
    pub(crate) mz_err: f32,
    pub(crate) ims_err: f32,
}

impl PartialOrd for SortableError {
    fn partial_cmp(&self, other: &SortableError) -> Option<std::cmp::Ordering> {
        self.intensity.partial_cmp(&other.intensity)
    }
}

impl Default for SortableError {
    fn default() -> Self {
        Self {
            intensity: 0.0,
            mz_err: f32::NAN,
            ims_err: f32::NAN,
        }
    }
}

impl MzMobilityOffsets {
    pub fn new(
        item: &SpectralCollector<IonAnnot, MzMobilityStatsCollector>,
        ref_mobility: f64,
    ) -> Self {
        let mut ms1 = TopNArray::new();
        let mut ms2 = TopNArray::new();

        for ((key, ref_mz), val) in item.iter_precursors() {
            if *key < 0i8 {
                continue;
            }
            let intensity = val.weight();
            let mz_err = (val.mean_mz().unwrap_or(f64::NAN) - ref_mz) as f32;
            // Make PPM
            let mz_err = mz_err / (*ref_mz as f32) * 1e6;
            let ims_err = (val.mean_mobility().unwrap_or(f64::NAN) - ref_mobility) as f32;
            // Make Pct
            let ims_err = ims_err / (ref_mobility as f32) * 1e2;
            ms1.push(SortableError {
                intensity,
                mz_err,
                ims_err,
            });
        }

        for ((_key, ref_mz), val) in item.iter_fragments() {
            let intensity = val.weight();
            let mz_err = (val.mean_mz().unwrap_or(f64::NAN) - ref_mz) as f32;
            let mz_err = mz_err / (*ref_mz as f32) * 1e6;
            let ims_err = (val.mean_mobility().unwrap_or(f64::NAN) - ref_mobility) as f32;
            let ims_err = ims_err / (ref_mobility as f32) * 1e2;
            ms2.push(SortableError {
                intensity,
                mz_err,
                ims_err,
            });
        }

        Self {
            ms1,
            ms2,
            ref_mobility,
        }
    }

    pub fn ms1_mz_errors(&self) -> [f32; 3] {
        let mut out = [0.0; 3];
        let vals = self.ms1.get_values();
        for i in 0..3 {
            out[i] = vals[i].mz_err;
        }
        out
    }

    pub fn ms2_mz_errors(&self) -> [f32; 7] {
        let mut out = [0.0; 7];
        let vals = self.ms2.get_values();
        for i in 0..7 {
            out[i] = vals[i].mz_err;
        }
        out
    }

    pub fn ms1_mobility_errors(&self) -> [f32; 3] {
        let mut out = [0.0; 3];
        let vals = self.ms1.get_values();
        for i in 0..3 {
            out[i] = vals[i].ims_err;
        }
        out
    }

    pub fn ms2_mobility_errors(&self) -> [f32; 7] {
        let mut out = [0.0; 7];
        let vals = self.ms2.get_values();
        for i in 0..7 {
            out[i] = vals[i].ims_err;
        }
        out
    }

    pub fn avg_delta_mobs(&self) -> (MzMobilityStatsCollector, MzMobilityStatsCollector) {
        let mut ms2 = MzMobilityStatsCollector::default();
        let mut ms1 = MzMobilityStatsCollector::default();
        let vals = self.ms2.get_values();
        for v in vals.iter().take(3) {
            if v.ims_err.is_nan() {
                continue;
            }
            ms2.add(v.intensity, v.mz_err as f64, v.ims_err as f64);
        }

        let vals = self.ms1.get_values();
        for v in vals.iter() {
            if v.ims_err.is_nan() {
                continue;
            }
            ms1.add(v.intensity, v.mz_err as f64, v.ims_err as f64);
        }

        (ms1, ms2)
    }
}
