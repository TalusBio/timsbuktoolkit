
#[derive(Debug, Clone, Copy)]
pub struct PeakInQuad {
    pub scan_index: u16,
    pub corrected_intensity: f32,
    pub retention_time_ms: u32,
    pub tof_index: u32,
}

