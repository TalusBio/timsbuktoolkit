use serde::Serialize;

#[derive(Debug, Clone, Copy, Serialize)]
pub struct RawPeak {
    pub scan_index: usize,
    pub tof_index: u32,
    pub intensity: u32,
    pub retention_time: f32,
}
