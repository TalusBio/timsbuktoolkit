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
