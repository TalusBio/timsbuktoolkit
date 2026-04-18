use serde::{
    Deserialize,
    Serialize,
};
use std::fmt;

/// Library reference retention time. Unit-agnostic — could be iRT, minutes,
/// or arbitrary units depending on the spectral library.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Serialize, Deserialize)]
#[serde(transparent)]
pub struct LibraryRT<T>(pub T);

/// Observed retention time from raw instrument data, always in seconds.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Serialize, Deserialize)]
#[serde(transparent)]
pub struct ObservedRTSeconds<T>(pub T);

impl From<LibraryRT<f32>> for LibraryRT<f64> {
    fn from(v: LibraryRT<f32>) -> Self {
        LibraryRT(v.0 as f64)
    }
}

impl From<ObservedRTSeconds<f32>> for ObservedRTSeconds<f64> {
    fn from(v: ObservedRTSeconds<f32>) -> Self {
        ObservedRTSeconds(v.0 as f64)
    }
}

impl<T: fmt::Display> fmt::Display for LibraryRT<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<T: fmt::Display> fmt::Display for ObservedRTSeconds<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}s", self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_library_rt_widening() {
        let narrow = LibraryRT(42.5f32);
        let wide: LibraryRT<f64> = narrow.into();
        assert!((wide.0 - 42.5).abs() < 1e-5);
    }

    #[test]
    fn test_observed_rt_widening() {
        let narrow = ObservedRTSeconds(123.4f32);
        let wide: ObservedRTSeconds<f64> = narrow.into();
        assert!((wide.0 - 123.4).abs() < 1e-3);
    }

    #[test]
    fn test_display() {
        assert_eq!(format!("{}", LibraryRT(42.5f64)), "42.5");
        assert_eq!(format!("{}", ObservedRTSeconds(42.5f64)), "42.5s");
    }

    #[test]
    fn test_serde_transparent() {
        let lib = LibraryRT(42.5f64);
        let json = serde_json::to_string(&lib).unwrap();
        assert_eq!(json, "42.5");
        let back: LibraryRT<f64> = serde_json::from_str(&json).unwrap();
        assert_eq!(back, lib);
    }
}
