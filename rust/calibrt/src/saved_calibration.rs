//! Saved calibration format, version validation, and file I/O.
use crate::{
    CalibrationSnapshot,
    RtAxis,
};

/// The only version [`SavedCalibration`] reads or writes. Named once so the
/// writer and the reader's gate cannot disagree.
pub const CALIBRATION_FORMAT_VERSION: &str = "v4";

/// JSON v4 calibration file format -- shared between CLI and viewer.
///
/// `calibration` is the grid's own snapshot and the only record of the fit: the
/// curve and the ridge widths are recomputed by refitting it, so the file cannot
/// carry a curve that disagrees with the points that produced it.
///
/// `R` is the residual block a writer measures beyond the curve. This crate does
/// not interpret it; a reader that has no use for it can leave it opaque.
#[derive(Debug, serde::Serialize, serde::Deserialize)]
pub struct SavedCalibration<R = serde_json::Value> {
    #[serde(deserialize_with = "deserialize_version")]
    pub version: String,
    pub rt_range_seconds: [f64; 2],
    /// Input-axis descriptor supplied by the library consumer. The fitter maps
    /// numerical library coordinates to observed seconds without interpreting units.
    pub library_rt_axis: RtAxis,
    pub calibration: CalibrationSnapshot,
    /// The uniform RT tolerance. Every writer has one -- it is what a query falls
    /// back to where the grid measured no ridge.
    pub rt_tolerance_minutes: f32,
    /// What a search measured beyond the curve. `None` for a writer that measures
    /// no residuals -- zeros there would read as "measured, and tight".
    ///
    /// The path is named because a bare `default` would demand `R: Default`.
    #[serde(default = "Option::default")]
    pub residuals: Option<R>,
    pub n_scored: usize,
}

impl<R> SavedCalibration<R> {
    /// Assemble a file record. The version is stamped here, not by callers.
    pub fn new(
        rt_range_seconds: [f64; 2],
        calibration: CalibrationSnapshot,
        rt_tolerance_minutes: f32,
        residuals: Option<R>,
        n_scored: usize,
    ) -> Self {
        Self {
            version: CALIBRATION_FORMAT_VERSION.to_string(),
            library_rt_axis: RtAxis::Unspecified,
            rt_range_seconds,
            calibration,
            rt_tolerance_minutes,
            residuals,
            n_scored,
        }
    }

    pub fn with_library_rt_axis(mut self, axis: RtAxis) -> Self {
        self.library_rt_axis = axis;
        self
    }

    /// Serialize to `path` in the layout [`Self::read`] expects.
    pub fn write(&self, path: &std::path::Path) -> Result<(), String>
    where
        R: serde::Serialize,
    {
        let json = serde_json::to_string_pretty(self).map_err(|e| e.to_string())?;
        std::fs::write(path, json).map_err(|e| e.to_string())
    }

    /// Parse a calibration file and check its provenance. The `Option<String>`
    /// is a reason to distrust the file, not an error: a calibration is only
    /// valid for the run it was fit on, and `raw_rt_range` -- the RT span of the
    /// run it is about to be used on -- is the one cheap way to catch the wrong
    /// file. `None` there means nothing verifies it, which also warns.
    pub fn read(
        path: &std::path::Path,
        raw_rt_range: Option<[f64; 2]>,
    ) -> Result<(Self, Option<String>), String>
    where
        R: serde::de::DeserializeOwned,
    {
        let json = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
        let saved: Self = serde_json::from_str(&json).map_err(|e| e.to_string())?;
        let warning = saved.provenance_warning(raw_rt_range);
        Ok((saved, warning))
    }

    /// Number of calibrant points the curve was fit on.
    pub fn n_calibrants(&self) -> usize {
        self.calibration.points.len()
    }

    fn provenance_warning(&self, raw_rt_range: Option<[f64; 2]>) -> Option<String> {
        let Some(raw) = raw_rt_range else {
            return Some(
                "No raw RT range to check the calibration against -- nothing verifies it was \
                 fit on this run"
                    .to_string(),
            );
        };
        let overlap_lo = self.rt_range_seconds[0].max(raw[0]);
        let overlap_hi = self.rt_range_seconds[1].min(raw[1]);
        let overlap = (overlap_hi - overlap_lo).max(0.0);
        let span = self.rt_range_seconds[1] - self.rt_range_seconds[0];
        if span <= 0.0 || overlap / span >= 0.5 {
            return None;
        }
        Some(format!(
            "Calibration RT range [{:.1}, {:.1}]s overlaps the raw file's [{:.1}, {:.1}]s by \
             {:.0}% -- it may have been fit on a different run",
            self.rt_range_seconds[0],
            self.rt_range_seconds[1],
            raw[0],
            raw[1],
            (overlap / span) * 100.0,
        ))
    }
}

fn deserialize_version<'de, D: serde::Deserializer<'de>>(
    deserializer: D,
) -> Result<String, D::Error> {
    let version = <String as serde::Deserialize>::deserialize(deserializer)?;
    if version != CALIBRATION_FORMAT_VERSION {
        return Err(serde::de::Error::custom(format!(
            "Unsupported calibration version: {version} (expected {CALIBRATION_FORMAT_VERSION})"
        )));
    }
    Ok(version)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> serde_json::Value {
        serde_json::json!({
            "version": "v4",
            "rt_range_seconds": [0.0, 1200.0],
            "library_rt_axis": {"kind": "normalized_index", "scale": "anchors"},
            "calibration": {"points": [], "grid_size": 16, "lookback": 4},
            "rt_tolerance_minutes": 1.25,
            "n_scored": 0
        })
    }

    #[test]
    fn typed_axis_preserves_the_v4_wire_format() {
        let json = fixture();
        let saved: SavedCalibration = serde_json::from_value(json.clone()).unwrap();
        assert_eq!(
            saved.library_rt_axis,
            RtAxis::NormalizedIndex {
                scale: Some("anchors".into())
            }
        );
        assert_eq!(
            serde_json::to_value(saved).unwrap()["library_rt_axis"],
            json["library_rt_axis"]
        );
    }

    #[test]
    fn malformed_axis_is_rejected_at_deserialization() {
        for axis in [
            serde_json::json!({"kind": "unknown"}),
            serde_json::json!({"kind": "normalized_index", "scale": 42}),
        ] {
            let mut json = fixture();
            json["library_rt_axis"] = axis;
            assert!(serde_json::from_value::<SavedCalibration>(json).is_err());
        }
    }

    #[test]
    fn direct_deserialization_checks_version_even_without_v4_fields() {
        for json in [
            r#"{"version":"v3","rt_range_seconds":[0,1200]}"#,
            r#"{"rt_range_seconds":[0,1200],"version":"v3"}"#,
        ] {
            let error = serde_json::from_str::<SavedCalibration>(json)
                .unwrap_err()
                .to_string();
            assert!(
                error.contains("Unsupported calibration version: v3 (expected v4)"),
                "{error}"
            );
        }
        let mut json = fixture();
        json["version"] = serde_json::json!("v2");
        assert!(
            serde_json::from_value::<SavedCalibration>(json)
                .unwrap_err()
                .to_string()
                .contains("Unsupported calibration version")
        );
    }
}
