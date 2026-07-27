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
//! frame), so the batch scrubber a future task wires in shows exactly one
//! frame and the Convergence tab's history is a single point — there was
//! only ever one fit to show, this being a Phase 2 snapshot rather than a
//! live run's sequence of batches.

use calib_dash::{
    CalibDash,
    CalibrantPoint,
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

    // A single frame never strides, so the budget only has to be big enough
    // for one frame's worth of points; `1 << 20` is generous for any
    // realistic calibrant count.
    let mut dash = CalibDash::new(
        1,
        n_calibrants,
        snapshot.grid_size,
        snapshot.lookback,
        1 << 20,
    );
    dash.on_batch(0, 0..n_calibrants, points.into_iter());
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
