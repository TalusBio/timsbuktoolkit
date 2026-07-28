//! Standalone replay of a saved `calibration.json` — the RT calibration
//! dashboard without a live Phase 1 search behind it.
//!
//! The JSON file's `"calibration"` field is read with `serde_json` alone rather
//! than through `timsseek`'s wrapper type: `calib_dash` must not pull in a whole
//! search engine crate to read three fields back out of its own save format.
//!
//! The loaded points become a single Phase-1-shaped batch (chunk 0, one frame),
//! so the batch scrubber shows one frame and the Convergence tab's history is a
//! single point — a Phase 2 snapshot is one fit, not a run's sequence of batches.

use calib_dash::{
    CalibDash,
    CalibrantPoint,
    REPLAY_BUDGET_BYTES,
};
use calibrt::CalibrationSnapshot;
use std::path::Path;

fn main() {
    let mut args = std::env::args_os().skip(1);
    let Some(path) = args.next() else {
        eprintln!("usage: calib_dash <calibration.json>");
        std::process::exit(2);
    };
    let path = Path::new(&path);

    let snapshot = match load_snapshot(path) {
        Ok(snapshot) => snapshot,
        Err(e) => {
            eprintln!("failed to load {}: {e}", path.display());
            std::process::exit(1);
        }
    };
    if let Err(e) = validate_snapshot(&snapshot) {
        eprintln!("invalid calibration snapshot in {}: {e}", path.display());
        std::process::exit(1);
    }

    // The saved points are `[library_rt, observed_rt, weight]`; the third column is
    // not read because every producer writes `calibrt::CALIBRANT_WEIGHT` there, which
    // is what the fit would use regardless.
    let points: Vec<CalibrantPoint> = snapshot
        .points
        .iter()
        .enumerate()
        .map(|(i, p)| CalibrantPoint {
            library_rt: p[0],
            observed_rt: p[1],
            speclib_index: i,
        })
        .collect();
    let n_calibrants = points.len();

    eprintln!(
        "calib_dash: loaded {n_calibrants} calibrants (grid_size={}, lookback={}) from {}",
        snapshot.grid_size,
        snapshot.lookback,
        path.display(),
    );

    let mut dash = CalibDash::new(
        1,
        n_calibrants,
        snapshot.grid_size,
        snapshot.lookback,
        REPLAY_BUDGET_BYTES,
    );
    dash.on_batch(0, points.into_iter());
    dash.finish();
}

/// Reads `path` and pulls `CalibrationSnapshot` out of its `"calibration"`
/// field without deserializing (or depending on the type of) the rest of
/// the saved file.
fn load_snapshot(path: &Path) -> Result<CalibrationSnapshot, String> {
    let json = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let value: serde_json::Value = serde_json::from_str(&json).map_err(|e| e.to_string())?;
    let calibration = value
        .get("calibration")
        .ok_or_else(|| "missing top-level \"calibration\" field".to_string())?;
    serde_json::from_value(calibration.clone()).map_err(|e| e.to_string())
}

/// Rejects the two shapes a syntactically valid `calibration.json` can still
/// carry that would otherwise reach a bare `.expect`/panic downstream:
/// `calibrt::Grid::new` rejects `bins == 0`, and `CalibrationCurve::new`
/// needs at least 2 points to define a curve.
fn validate_snapshot(snapshot: &CalibrationSnapshot) -> Result<(), String> {
    if snapshot.grid_size == 0 {
        return Err("grid_size is 0 — a calibration grid needs at least 1 bin".to_string());
    }
    if snapshot.points.len() < 2 {
        return Err(format!(
            "only {} calibrant point(s) in \"calibration.points\" — at least 2 are needed to \
             fit a curve",
            snapshot.points.len()
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn a_file_without_a_calibration_field_names_the_missing_field() {
        // The handle is held (not just its path) because dropping it is what
        // deletes the file — including when the test fails.
        let mut f = NamedTempFile::new().expect("a writable temp dir");
        f.write_all(br#"{"version": 1, "other": {}}"#)
            .expect("write the fixture");
        let err = load_snapshot(f.path()).expect_err("there is no \"calibration\" field");
        assert!(
            err.contains("calibration"),
            "message must name the missing field: {err}"
        );
    }

    #[test]
    fn validate_snapshot_accepts_only_a_fittable_snapshot() {
        // (grid_size, n_points, expect_ok)
        let cases = [(10, 2, true), (0, 2, false), (10, 1, false)];
        for (grid_size, n_points, expect_ok) in cases {
            let snapshot = CalibrationSnapshot {
                points: vec![[1.0, 1.0, 1.0]; n_points],
                grid_size,
                lookback: 3,
            };
            match (validate_snapshot(&snapshot), expect_ok) {
                (Ok(()), true) | (Err(_), false) => {}
                (got, _) => panic!("grid_size {grid_size}, {n_points} points: got {got:?}"),
            }
        }
    }
}
