//! Standalone replay of a saved `calibration.json` — the RT calibration
//! dashboard without a live Phase 1 search behind it.
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
            library_id: i as u64,
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

/// Reads `path` and keeps only the snapshot. The residual block is left opaque —
/// nothing here reads it — and the provenance warning is dropped: there is no raw
/// file to check the calibration against.
fn load_snapshot(path: &Path) -> Result<CalibrationSnapshot, String> {
    let (saved, _) = calibrt::SavedCalibration::<serde_json::Value>::read(path, None)?;
    Ok(saved.calibration)
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

    /// A file whose layout parses but whose version is not the current one is
    /// refused: its fields may mean something else, and the snapshot is refit
    /// here rather than validated field by field.
    #[test]
    fn a_file_of_a_foreign_version_is_refused() {
        // The handle is held (not just its path) because dropping it is what
        // deletes the file — including when the test fails.
        let mut f = NamedTempFile::new().expect("a writable temp dir");
        f.write_all(
            br#"{
              "version": "v2",
              "rt_range_seconds": [0.0, 1200.0],
              "calibration": {
                "points": [[0.0, 30.0, 1.0], [100.0, 200.0, 1.0]],
                "grid_size": 16,
                "lookback": 4
              },
              "rt_tolerance_minutes": 1.25,
              "n_scored": 999
            }"#,
        )
        .expect("write the fixture");
        let err = load_snapshot(f.path()).expect_err("v2 is not the current format");
        assert!(
            err.contains("v2") && err.contains(calibrt::CALIBRATION_FORMAT_VERSION),
            "message must name both versions: {err}"
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
