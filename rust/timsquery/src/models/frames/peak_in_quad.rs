use super::raw_peak::RawPeak;

#[derive(Debug, Clone, Copy)]
pub struct PeakInQuad {
    pub scan_index: usize,
    pub intensity: u32,
    pub retention_time: f32,
    pub tof_index: u32,
}

impl From<PeakInQuad> for RawPeak {
    fn from(peak_in_quad: PeakInQuad) -> Self {
        RawPeak {
            scan_index: peak_in_quad.scan_index,
            tof_index: peak_in_quad.tof_index,
            intensity: peak_in_quad.intensity,
            retention_time: peak_in_quad.retention_time,
        }
    }
}
