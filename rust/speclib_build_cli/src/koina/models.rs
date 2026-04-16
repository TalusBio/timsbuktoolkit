use serde::{Deserialize, Serialize};

// ── Model registry ───────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FragmentModel {
    Prosit2020IntensityHcd,
    AlphaPeptDeepMs2Generic,
}

impl FragmentModel {
    pub fn from_name(name: &str) -> Result<Self, String> {
        match name {
            "Prosit_2020_intensity_HCD" => Ok(Self::Prosit2020IntensityHcd),
            "AlphaPeptDeep_ms2_generic" => Ok(Self::AlphaPeptDeepMs2Generic),
            other => Err(format!(
                "unknown fragment model {other:?}; valid: Prosit_2020_intensity_HCD, AlphaPeptDeep_ms2_generic"
            )),
        }
    }

    pub fn model_name(&self) -> &str {
        match self {
            Self::Prosit2020IntensityHcd => "Prosit_2020_intensity_HCD",
            Self::AlphaPeptDeepMs2Generic => "AlphaPeptDeep_ms2_generic",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RtModel {
    Prosit2019Irt,
    AlphaPeptDeepRtGeneric,
}

impl RtModel {
    pub fn from_name(name: &str) -> Result<Self, String> {
        match name {
            "Prosit_2019_irt" => Ok(Self::Prosit2019Irt),
            "AlphaPeptDeep_rt_generic" => Ok(Self::AlphaPeptDeepRtGeneric),
            other => Err(format!(
                "unknown RT model {other:?}; valid: Prosit_2019_irt, AlphaPeptDeep_rt_generic"
            )),
        }
    }

    pub fn model_name(&self) -> &str {
        match self {
            Self::Prosit2019Irt => "Prosit_2019_irt",
            Self::AlphaPeptDeepRtGeneric => "AlphaPeptDeep_rt_generic",
        }
    }
}

// ── Triton v2 wire types ─────────────────────────────────────────────────────

#[derive(Debug, Serialize)]
pub struct KoinaRequest {
    pub id: String,
    pub inputs: Vec<KoinaTensor>,
}

#[derive(Debug, Serialize)]
pub struct KoinaTensor {
    pub name: String,
    pub shape: Vec<usize>,
    pub datatype: String,
    pub data: KoinaTensorData,
}

#[derive(Debug, Serialize)]
#[serde(untagged)]
pub enum KoinaTensorData {
    Strings(Vec<String>),
    Ints(Vec<i32>),
    Floats(Vec<f32>),
}

#[derive(Debug, Deserialize)]
pub struct KoinaResponse {
    pub id: String,
    pub outputs: Vec<KoinaOutputTensor>,
}

#[derive(Debug, Deserialize)]
pub struct KoinaOutputTensor {
    pub name: String,
    pub datatype: String,
    pub shape: Vec<usize>,
    pub data: Vec<serde_json::Value>,
}

// ── Domain prediction types ──────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct PredictionInput {
    pub sequence: String,
    pub charge: u8,
    pub nce: f32,
}

#[derive(Debug, Clone)]
pub struct FragmentPrediction {
    /// Ion annotations, e.g. "y3^1", "b5^2".
    pub annotations: Vec<String>,
    pub mzs: Vec<f64>,
    pub intensities: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct RtPrediction {
    pub irt: f32,
}
