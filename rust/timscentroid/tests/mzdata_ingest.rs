//! mzdata ingest validation.
//!
//! The real-file test is `#[ignore]` (data lives outside the repo). Run it
//! against a real Thermo DIA mzML with:
//!
//!   MZML_TEST_PATH=/path/to/file.mzML \
//!     cargo test -p timscentroid --features mzdata --test mzdata_ingest -- --ignored --nocapture

#![cfg(feature = "mzdata")]

use timscentroid::MobilityKind;
use timscentroid::reader::mzdata::from_mzml_file;

#[test]
#[ignore = "requires a local Thermo DIA mzML via MZML_TEST_PATH"]
fn ingests_real_thermo_dia_mzml() {
    let path =
        std::env::var("MZML_TEST_PATH").expect("set MZML_TEST_PATH to a local Thermo DIA .mzML");
    let cfg = timscentroid::IndexingCentroidingConfig::default();
    let idx = from_mzml_file(std::path::Path::new(&path), &cfg)
        .expect("from_mzml_file should succeed on a centroided DIA mzML");

    // Thermo Astral DIA carries FAIMS CV (value 0) → Unsupported, NOT Ook0.
    println!("mobility_kind = {:?}", idx.mobility_kind());
    assert!(
        !matches!(idx.mobility_kind(), MobilityKind::Ook0),
        "mzML must never resolve to a searchable TIMS axis"
    );

    let n_windows = idx.n_ms2_window_groups();
    let n_cycles = idx.ms1_cycle_mapping().len();
    println!("windows={n_windows} cycles={n_cycles}");
    assert!(n_windows > 0, "expected at least one isolation window");
    assert!(n_cycles > 0, "expected at least one cycle");
}
