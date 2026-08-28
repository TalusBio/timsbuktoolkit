//! Input half of the contract in `README.md`, next to this file.
//!
//! Carafe parses with fastjson: a renamed field is a null on their side, not
//! an error. Payloads are copied from the README and go in through the public
//! entry points, so a change at the boundary fails here.
//!
//! Output half: `timsquery_cli`'s `carafe_output_contract`. End-to-end,
//! including the `-o` layout and the exit code: `timsquery_cli`'s
//! `tests/carafe_e2e.rs`, which runs the built binary over a generated mzML.

use std::io::Write;
use timsquery::models::SourceId;
use timsquery::models::tolerance::{
    MobilityTolerance,
    MzTolerance,
    QuadTolerance,
    RtTolerance,
    Tolerance,
};
use timsquery::serde::{
    TargetTable,
    read_targets,
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
    let arena = read_targets(f.path()).expect("Carafe's target JSON must load");

    // Labels must resolve to ion annotations, not to the string arena.
    let TargetTable::Mzpaf { geom, .. } = arena else {
        panic!("labelled targets must land in the ion-annotated arena");
    };
    assert_eq!(geom.n_rows(), 1);
    let row = geom.rows().next().unwrap();
    assert_eq!(geom.precursor_mz(row), 650.32);
    assert_eq!(geom.charge(row), 2);
    assert_eq!(geom.rt_seconds(row), 1234.5);
    assert_eq!(geom.mobility(row), 0.95);

    // Positionally paired with the m/z values; the `^N` suffix survives.
    assert_eq!(geom.frag_mzs(row), [175.1, 288.2]);
    let labels: Vec<String> = geom
        .frag_labels(row)
        .iter()
        .map(|l| l.to_string())
        .collect();
    assert_eq!(labels, vec!["y1".to_string(), "y3^2".to_string()]);
}

/// Invariant 1: `id` is echoed, not renumbered. Ids that are not `0..n-1` in
/// file order used to come back as the arena position, silently relabelling
/// every row.
#[test]
fn non_sequential_ids_survive_the_reader() {
    let f = write_targets(
        r#"[
      { "id": 7, "mobility": 0.95, "rt_seconds": 1234.5, "precursor": 650.32,
        "precursor_charge": 2, "precursor_isotopes": [0],
        "fragments": [175.1], "fragment_labels": ["y1"] },
      { "id": 42, "mobility": 0.80, "rt_seconds": 999.0, "precursor": 500.0,
        "precursor_charge": 2, "precursor_isotopes": [0],
        "fragments": [200.0], "fragment_labels": ["y2"] }
    ]"#,
    );
    let arena = read_targets(f.path()).expect("loads");
    let TargetTable::Mzpaf { geom, .. } = arena else {
        panic!("mzpaf labels")
    };
    let ids: Vec<SourceId<'_>> = geom.rows().map(|r| geom.output_id(r)).collect();
    assert_eq!(
        ids,
        vec![SourceId::Numeric(7), SourceId::Numeric(42)],
        "the caller's ids, not 0..n-1"
    );
}

/// Carafe keys results by `id`, so a repeat makes one row unreachable.
#[test]
fn duplicate_ids_are_rejected_by_the_reader() {
    let f = write_targets(
        r#"[
      { "id": 7, "mobility": 0.95, "rt_seconds": 1234.5, "precursor": 650.32,
        "precursor_charge": 2, "precursor_isotopes": [0],
        "fragments": [175.1], "fragment_labels": ["y1"] },
      { "id": 7, "mobility": 0.80, "rt_seconds": 999.0, "precursor": 500.0,
        "precursor_charge": 2, "precursor_isotopes": [0],
        "fragments": [200.0], "fragment_labels": ["y2"] }
    ]"#,
    );
    assert!(read_targets(f.path()).is_err());
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
        read_targets(f.path()).is_err(),
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
