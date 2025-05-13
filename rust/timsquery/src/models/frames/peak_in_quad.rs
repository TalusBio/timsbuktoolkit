use timsrust::converters::{
    ConvertableDomain,
    Scan2ImConverter,
    Tof2MzConverter,
};

#[derive(Debug, Clone, Copy)]
pub struct PeakInQuad {
    pub scan_index: u16,
    pub corrected_intensity: f32,
    pub retention_time_ms: u32,
    pub tof_index: u32,
}

#[derive(Debug, Clone, Copy)]
pub struct ResolvedPeakInQuad {
    pub mobility: f32,
    pub corrected_intensity: f32,
    pub retention_time_ms: u32,
    // Q: Should this be an f64?
    pub mz: f32,
}

impl PeakInQuad {
    pub fn resolve(
        self,
        im_converter: &Scan2ImConverter,
        tof_converter: &Tof2MzConverter,
    ) -> ResolvedPeakInQuad {
        let mobility = im_converter.convert(self.scan_index as f64) as f32;
        let mz = tof_converter.convert(self.tof_index as f64) as f32;
        ResolvedPeakInQuad {
            mobility,
            corrected_intensity: self.corrected_intensity,
            retention_time_ms: self.retention_time_ms,
            mz,
        }
    }
}
