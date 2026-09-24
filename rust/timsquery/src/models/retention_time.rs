//! Library coordinates and observed extraction times are different domains.
use serde::{
    Deserialize,
    Serialize,
};

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum RtAxis {
    #[default]
    Absent,
    Seconds,
    NormalizedIndex {
        scale: Option<String>,
    },
    Unspecified,
}

impl std::fmt::Display for RtAxis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Absent => f.write_str("no RT"),
            Self::Seconds => f.write_str("s"),
            Self::NormalizedIndex { scale: Some(scale) } => write!(f, "index ({scale})"),
            Self::NormalizedIndex { scale: None } => f.write_str("index"),
            Self::Unspecified => f.write_str("unspecified units"),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct RtCoordinate<'a> {
    pub value: f32,
    pub axis: &'a RtAxis,
}

impl RtCoordinate<'static> {
    pub fn seconds(value: f32) -> Self {
        Self {
            value,
            axis: &RtAxis::Seconds,
        }
    }
}

/// A borrowed geometry view at a calibrated or observed apex time.
#[derive(Debug, Clone, Copy)]
pub struct AtObservedRt<'a, T: crate::traits::Target + ?Sized> {
    inner: &'a T,
    seconds: f32,
}
impl<'a, T: crate::traits::Target + ?Sized> AtObservedRt<'a, T> {
    pub fn new(inner: &'a T, seconds: f32) -> Self {
        Self { inner, seconds }
    }

    pub fn inner(&self) -> &'a T {
        self.inner
    }
}
impl<T: crate::traits::Target + ?Sized> crate::traits::Target for AtObservedRt<'_, T> {
    type Label = T::Label;

    fn source_id(&self) -> Option<crate::models::SourceId<'_>> {
        self.inner.source_id()
    }

    fn output_id(&self) -> crate::models::SourceId<'_> {
        self.inner.output_id()
    }

    fn mono_precursor_mz(&self) -> f64 {
        self.inner.mono_precursor_mz()
    }

    fn precursor_charge(&self) -> u8 {
        self.inner.precursor_charge()
    }

    fn rt(&self) -> Option<RtCoordinate<'_>> {
        Some(RtCoordinate::seconds(self.seconds))
    }

    fn mobility_ook0(&self) -> f32 {
        self.inner.mobility_ook0()
    }

    fn precursor_mz_limits(&self) -> (f64, f64) {
        self.inner.precursor_mz_limits()
    }

    fn precursor_count(&self) -> usize {
        self.inner.precursor_count()
    }

    fn fragment_count(&self) -> usize {
        self.inner.fragment_count()
    }

    fn iter_precursors(&self) -> impl Iterator<Item = (i8, f64)> {
        self.inner.iter_precursors()
    }

    fn iter_fragments_refs(&self) -> impl Iterator<Item = (&Self::Label, f64)> {
        self.inner.iter_fragments_refs()
    }
}
