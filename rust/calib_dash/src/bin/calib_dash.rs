//! Standalone replay of a saved `calibration.json` — the RT calibration
//! dashboard without a live Phase 1 search behind it.
//!
//! Loads `calibrt::CalibrationSnapshot` (points, grid size, DP lookback)
//! directly from the JSON file's `"calibration"` field. That field's shape
//! matches `timsseek::rt_calibration::CalibrationResult::save_json`'s output
//! one-for-one, but this binary deserializes it with `serde_json` alone
//! rather than depending on `timsseek` for the wrapper type
//! (`SavedCalibration`) that field lives inside — `calib_dash` must not pull
//! in a whole search engine crate just to read three fields back out of its
//! own save format.
//!
//! The loaded points become a single Phase-1-shaped batch (chunk 0, one
//! frame), so the `<`/`>` batch scrubber (`App`'s `scrub_frame`/
//! `CalibDash::sync_scrub`) shows exactly one frame and the Convergence
//! tab's history is a single point — there was only ever one fit to show,
//! this being a Phase 2 snapshot rather than a live run's sequence of
//! batches.

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

    let points: Vec<CalibrantPoint> = snapshot
        .points
        .iter()
        .enumerate()
        .map(|(i, p)| CalibrantPoint {
            library_rt: p[0],
            observed_rt: p[1],
            score: p[2],
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
    dash.finish(0);
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
/// carry that would otherwise reach a bare `.expect`/panic downstream rather
/// than a message a human can act on:
///
/// - `grid_size == 0`: `calibrt::Grid::new` rejects `bins == 0` before it
///   even looks at the fit range, so this would abort the process inside
///   `CalibDash::new` (which has no `Result` to report it through — see its
///   own doc comment).
/// - fewer than 2 points: `calibrt::CalibrationCurve::new` itself requires
///   at least 2 points to define a curve at all (`CalibRtError::InsufficientPoints`);
///   below that there is nothing to replay, so this is caught here rather
///   than surfacing as an empty dashboard with no explanation.
///
/// `lookback` needs no check: `calibrt::pathfinding::find_optimal_path`
/// only ever uses it via `saturating_sub`, so every value including `0`
/// degrades gracefully (each node considers only itself) rather than
/// panicking or looping unboundedly.
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

    /// Writes `contents` to a fresh temp file, so each test gets an isolated
    /// file `load_snapshot` can read without the tests stepping on each
    /// other's fixtures. The handle is returned (not just its path) because
    /// dropping it is what deletes the file — including when the test fails.
    fn temp_file(contents: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().expect("a writable temp dir");
        f.write_all(contents.as_bytes()).expect("write the fixture");
        f
    }

    #[test]
    fn missing_file_reports_a_message_not_a_panic() {
        let path = Path::new("/definitely/does/not/exist/calibration.json");
        let err = load_snapshot(path).expect_err("the file does not exist");
        assert!(
            !err.is_empty(),
            "must carry a human-readable reason, not just fail silently"
        );
    }

    #[test]
    fn malformed_json_reports_a_message_not_a_panic() {
        let f = temp_file("{ not valid json ");
        let err = load_snapshot(f.path()).expect_err("the file is not valid JSON");
        assert!(!err.is_empty());
    }

    #[test]
    fn a_zero_grid_size_is_rejected_before_it_can_panic() {
        // Syntactically valid — this is exactly the file this task's binary
        // must survive: `CalibDash::new`'s `.expect` has no `Result` to
        // report a zero `bins` through, so `validate_snapshot` must catch it
        // first.
        let snapshot = CalibrationSnapshot {
            points: vec![[1.0, 1.0, 1.0], [2.0, 2.0, 1.0]],
            grid_size: 0,
            lookback: 3,
        };
        let err = validate_snapshot(&snapshot).expect_err("grid_size 0 must be rejected");
        assert!(
            err.contains("grid_size"),
            "message must name the field: {err}"
        );
    }

    #[test]
    fn too_few_points_is_rejected_with_an_actionable_message() {
        let snapshot = CalibrationSnapshot {
            points: vec![[1.0, 1.0, 1.0]],
            grid_size: 10,
            lookback: 3,
        };
        let err = validate_snapshot(&snapshot).expect_err("one point cannot define a curve");
        assert!(
            err.contains('1'),
            "message should state how many points were actually found: {err}"
        );
    }

    #[test]
    fn a_well_formed_snapshot_passes_validation() {
        let snapshot = CalibrationSnapshot {
            points: vec![[1.0, 1.0, 1.0], [2.0, 2.0, 1.0]],
            grid_size: 10,
            lookback: 3,
        };
        assert!(validate_snapshot(&snapshot).is_ok());
    }
}
