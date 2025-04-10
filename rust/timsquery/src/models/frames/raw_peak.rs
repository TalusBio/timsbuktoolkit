use serde::Serialize;
use std::convert::From;

#[derive(Debug, Clone, Copy, Serialize)]
pub struct RawPeak {
    pub scan_index: usize,
    pub tof_index: u32,
    pub intensity: u32,
    pub retention_time: f32,
}

impl<T> From<(RawPeak, T)> for RawPeak {
    fn from(x: (RawPeak, T)) -> Self {
        x.0
    }
}
