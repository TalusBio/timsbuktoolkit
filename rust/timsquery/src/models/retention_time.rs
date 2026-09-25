//! Library coordinates and observed extraction times are different domains.
pub use calibrt::RtAxis;

#[derive(Debug, Clone, Copy)]
pub struct RtCoordinate<'a> {
    pub value: calibrt::LibraryRT<f32>,
    pub axis: &'a RtAxis,
}

impl RtCoordinate<'static> {
    pub fn seconds(value: f32) -> Self {
        Self {
            value: calibrt::LibraryRT(value),
            axis: &RtAxis::Seconds,
        }
    }
}
