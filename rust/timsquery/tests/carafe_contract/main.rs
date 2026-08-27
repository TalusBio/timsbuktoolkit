//! Input half of the contract in `README.md`, next to this file.
//!
//! Carafe parses with fastjson: a renamed field is a null on their side, not
//! an error. Payloads are copied from the README and go in through the public
//! entry points, so a change at the boundary fails here.
//!
//! Output half: `timsquery_cli`'s `carafe_output_contract`. Not covered: the
//! `-o` directory layout and the exit code, which need the built binary.

use std::io::Write;
use timsquery::models::tolerance::{
    MobilityTolerance,
    MzTolerance,
    QuadTolerance,
    RtTolerance,
    Tolerance,
};
use timsquery::serde::{
    LibraryArena,
    read_library_file,
};

/// README, "Targets (`-e`, `psm_query.json`)".
const CARAFE_TARGETS: &str = r#"[
  { "id": 0, "mobility": 0.95, "rt_seconds": 1234.5, "precursor": 650.32,
    "precursor_charge": 2, "precursor_isotopes": [0,1,2],
    "fragments": [175.1, 288.2], "fragment_labels": ["y1","y3^2"] }
]"#;

/// README, "Tolerances (`-t`)". Carafe spells these `percent` and `absolute`;
/// the CLI's own templates use `pct` and `da`.
const CARAFE_TOLERANCES: &str = r#"{
  "ms":       { "ppm":      [13.0, 17.0] },
  "rt":       { "minutes":  [0.1, 0.1] },
  "mobility": { "percent":  [3.0, 3.0] },
  "quad":     { "absolute": [0.1, 0.1] }
}"#;

fn write_targets(json: &str) -> tempfile::NamedTempFile {
    let mut f = tempfile::Builder::new()
        .suffix(".json")
        .tempfile()
        .expect("tempfile");
    f.write_all(json.as_bytes()).expect("write targets");
    f
}

/// `precursor`, `fragments` and `fragment_labels` are serde aliases; unifying
/// them with the canonical names would break Carafe and nothing else.
#[test]
fn carafe_target_payload_loads_through_the_public_reader() {
    let f = write_targets(CARAFE_TARGETS);
    let arena = read_library_file(f.path()).expect("Carafe's target JSON must load");

    // Labels must resolve to ion annotations, not to the string arena.
    let LibraryArena::Mzpaf { geom, .. } = arena else {
        panic!("labelled targets must land in the ion-annotated arena");
    };
    assert_eq!(geom.n_rows(), 1);
    assert_eq!(geom.precursor_mz[0], 650.32);
    assert_eq!(geom.charge[0], 2);
    assert_eq!(geom.rt_seconds[0], 1234.5);
    assert_eq!(geom.mobility[0], 0.95);

    // Positionally paired with the m/z values; the `^N` suffix survives.
    assert_eq!(geom.frag_mzs.len(), 2);
    assert_eq!(geom.frag_mzs[0], 175.1);
    assert_eq!(geom.frag_mzs[1], 288.2);
    let labels: Vec<String> = geom.frag_labels.iter().map(|l| l.to_string()).collect();
    assert_eq!(labels, vec!["y1".to_string(), "y3^2".to_string()]);
}

/// `id` is Carafe's map key; defaulting it NPEs downstream instead of erroring.
#[test]
fn carafe_id_field_is_required() {
    // Spelled out, not surgered from the const: malformed JSON would "pass".
    let without_id = r#"[
      { "mobility": 0.95, "rt_seconds": 1234.5, "precursor": 650.32,
        "precursor_charge": 2, "precursor_isotopes": [0,1,2],
        "fragments": [175.1, 288.2], "fragment_labels": ["y1","y3^2"] }
    ]"#;
    let f = write_targets(without_id);
    assert!(
        read_library_file(f.path()).is_err(),
        "a target without `id` must fail rather than default it"
    );
}

/// Dropping either spelling silently breaks one of the two callers.
#[test]
fn carafe_tolerance_spellings_deserialize() {
    let tol: Tolerance =
        serde_json::from_str(CARAFE_TOLERANCES).expect("Carafe's tolerance JSON must deserialize");

    // Asymmetric by design: the edges encode a calibration offset (README).
    // Collapsing them to one symmetric tolerance recentres every window.
    assert_eq!(tol.ms, MzTolerance::Ppm((13.0, 17.0)));
    assert_eq!(tol.rt, RtTolerance::Minutes((0.1, 0.1)));
    assert_eq!(tol.mobility, MobilityTolerance::Pct((3.0, 3.0)));
    assert_eq!(tol.quad, QuadTolerance::Absolute((0.1, 0.1)));
}

/// When Carafe's calibration offset exceeds the half-width the low edge goes
/// negative — a window sitting entirely above the target mass. This already
/// worked; the test pins it so a future "both edges must be positive" check
/// cannot be added without failing. `mz_range`'s invariant is `low + high >= 0`.
#[test]
fn a_negative_low_tolerance_edge_is_accepted() {
    let tol: Tolerance = serde_json::from_str(
        r#"{
      "ms":       { "ppm":      [-2.0, 32.0] },
      "rt":       { "minutes":  [0.1, 0.1] },
      "mobility": { "percent":  [3.0, 3.0] },
      "quad":     { "absolute": [0.1, 0.1] }
    }"#,
    )
    .expect("a negative low edge must deserialize");
    assert_eq!(tol.ms, MzTolerance::Ppm((-2.0, 32.0)));

    let range = tol.mz_range(1000.0);
    assert!(
        range.start() > 1000.0 && range.end() > range.start(),
        "expected a well-formed window above the target, got {:?}..{:?}",
        range.start(),
        range.end()
    );
}
