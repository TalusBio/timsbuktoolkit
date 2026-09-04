//! End-to-end half of the Carafe contract in `rust/timsquery/tests/carafe_contract/`.
//!
//! The other two halves check the reader and the serializer in isolation. This
//! one runs the built binary against a generated mzML, so it covers the parts
//! they explicitly cannot: the `-o` directory layout, the `results.json`
//! basename and the exit code.

mod mini_dia;

use std::path::Path;
use std::process::Command;

/// Carafe's target payload, spelled as in the contract README: `precursor`,
/// `fragments` and `fragment_labels` are the aliases it sends.
fn targets_json() -> String {
    format!(
        r#"[
  {{ "id": 7, "mobility": 0.95, "rt_seconds": {rt}, "precursor": {mz},
     "precursor_charge": {charge}, "precursor_isotopes": [0,1,2],
     "fragments": [{f0}, {f1}], "fragment_labels": ["y1","y3^2"] }}
]"#,
        rt = mini_dia::RT_START_S + (mini_dia::N_CYCLES / 2) as f64 * mini_dia::RT_STEP_S,
        mz = mini_dia::PRECURSOR_MZ,
        charge = mini_dia::PRECURSOR_CHARGE,
        f0 = mini_dia::FRAGMENT_MZS[0],
        f1 = mini_dia::FRAGMENT_MZS[1],
    )
}

/// Wide enough that the fixture's exact values are inside every window; this
/// test is about the contract, not about tolerance behaviour.
const TOLERANCES: &str = r#"{
  "ms":       { "ppm":      [50.0, 50.0] },
  "rt":       { "minutes":  [1.0, 1.0] },
  "mobility": { "percent":  [50.0, 50.0] },
  "quad":     { "absolute": [1.0, 1.0] }
}"#;

fn write(dir: &Path, name: &str, contents: &str) -> std::path::PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, contents).expect("write fixture");
    path
}

#[test]
fn carafe_cli_writes_results_json_and_exits_zero() {
    let raw = mini_dia::path();
    let tmp = tempfile::tempdir().expect("tempdir");
    let targets = write(tmp.path(), "psm_query.json", &targets_json());
    let tolerances = write(tmp.path(), "tolerances.json", TOLERANCES);
    let out_dir = tmp.path().join("out");

    let status = Command::new(env!("CARGO_BIN_EXE_timsquery_cli"))
        .args([
            "query-index",
            "--raw-file-path",
            raw.to_str().unwrap(),
            "--tolerance-settings-path",
            tolerances.to_str().unwrap(),
            "--elution-groups-path",
            targets.to_str().unwrap(),
            "--output-path",
            out_dir.to_str().unwrap(),
            "--aggregator",
            "chromatogram-aggregator",
            "--format",
            "ndjson",
        ])
        .output()
        .expect("run timsquery_cli");

    assert!(
        status.status.success(),
        "cli must exit 0; stderr:\n{}",
        String::from_utf8_lossy(&status.stderr)
    );

    // Invariant 5: the basename is exactly `results.json`, inside `-o`.
    let results = out_dir.join("results.json");
    assert!(
        results.is_file(),
        "expected {} to exist; `-o` held: {:?}",
        results.display(),
        std::fs::read_dir(&out_dir).map(|d| d.flatten().map(|e| e.file_name()).collect::<Vec<_>>()),
    );

    let body = std::fs::read_to_string(&results).expect("read results");
    let lines: Vec<&str> = body.lines().filter(|l| !l.trim().is_empty()).collect();
    assert_eq!(lines.len(), 1, "one target with signal -> one result line");

    // Invariant 2: one complete object per line, no array wrapper.
    let row: serde_json::Value = serde_json::from_str(lines[0]).expect("ndjson object per line");
    let row = row.as_object().expect("an object, not an array");

    // Invariant 1: the caller's id, echoed.
    assert_eq!(row["id"], 7, "the id Carafe sent, not the row position");

    // Invariant 3: the chromatogram spelling, not the spectrum one.
    for key in [
        "precursor_mzs",
        "precursor_intensities",
        "fragment_mzs",
        "fragment_labels",
        "fragment_intensities",
        "retention_time_results_seconds",
    ] {
        assert!(
            row.contains_key(key),
            "missing `{key}`; keys: {:?}",
            row.keys().collect::<Vec<_>>()
        );
    }

    // Invariant 4: every intensity row is as long as the RT axis.
    let n_rt = row["retention_time_results_seconds"]
        .as_array()
        .unwrap()
        .len();
    for matrix in ["precursor_intensities", "fragment_intensities"] {
        for (i, ion) in row[matrix].as_array().unwrap().iter().enumerate() {
            assert_eq!(
                ion.as_array().unwrap().len(),
                n_rt,
                "{matrix}[{i}] is not one value per RT point",
            );
        }
    }

    // The intensities the fixture wrote come back unchanged, so a pass here
    // cannot be the empty-chromatogram path in disguise. The first fragment is
    // written at 0.8x the cycle apex.
    let apex = row["fragment_intensities"][0]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .fold(f64::MIN, f64::max);
    assert!(
        (apex - mini_dia::apex_intensity() * 0.8).abs() < 1e-3,
        "extracted apex {apex} is not the intensity the fixture wrote",
    );
}
