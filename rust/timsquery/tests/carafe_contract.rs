//! Executable form of `docs/CARAFE_CONTRACT.md`.
//!
//! Carafe drives `timsquery_cli` as a subprocess and parses its output with
//! fastjson, with no field remapping and no schema negotiation. A renamed
//! field or a dropped serde alias fails loudly on neither side — Carafe just
//! gets a null and NPEs somewhere unrelated. These tests pin the mechanically
//! checkable parts, using the literal payloads from the document so the two
//! cannot drift apart silently.
//!
//! Everything goes through the public entry points Carafe actually reaches
//! (a file path into `read_library_file`, a JSON blob into `Tolerance`) rather
//! than internal types, so a refactor that keeps the internals working but
//! changes the boundary still fails here.
//!
//! This file covers the INPUT direction. The output direction — result field
//! names, the `-a`/`-f` flag values and the `results.json` basename — is
//! pinned by `timsquery_cli`'s `carafe_output_contract` module, which has to
//! live inside that crate because it has no library target.
//!
//! NOT covered, because it needs the built binary and a real `.d`: the `-o`
//! directory layout and the process exit code.

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

/// Verbatim from the contract's "Targets (`-e`, `psm_query.json`)" section.
const CARAFE_TARGETS: &str = r#"[
  { "id": 0, "mobility": 0.95, "rt_seconds": 1234.5, "precursor": 650.32,
    "precursor_charge": 2, "precursor_isotopes": [0,1,2],
    "fragments": [175.1, 288.2], "fragment_labels": ["y1","y3^2"] }
]"#;

/// Verbatim from the contract's "Tolerances (`-t`)" section. Note `percent`
/// and `absolute`, not the `pct`/`da` spellings the CLI's own templates use —
/// both must deserialize.
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
    f.flush().expect("flush");
    f
}

/// Every field name Carafe emits must land somewhere. `precursor`,
/// `fragments` and `fragment_labels` are serde aliases of names used elsewhere
/// in the codebase, so a cleanup that "unifies" them would break Carafe
/// without failing any other test.
#[test]
fn carafe_target_payload_loads_through_the_public_reader() {
    let f = write_targets(CARAFE_TARGETS);
    let arena = read_library_file(f.path()).expect("Carafe's target JSON must load");

    // String fragment labels must resolve to ion annotations, not to the
    // string-labelled arena, which carries no ion chemistry.
    let LibraryArena::Mzpaf { geom, .. } = arena else {
        panic!("labelled targets must land in the ion-annotated arena");
    };
    assert_eq!(geom.n_rows(), 1);
    assert_eq!(geom.precursor_mz[0], 650.32);
    assert_eq!(geom.charge[0], 2);
    assert_eq!(geom.rt_seconds[0], 1234.5);
    assert_eq!(geom.mobility[0], 0.95);

    // Labels are positionally paired with their m/z values, and the `^N`
    // charge suffix survives.
    assert_eq!(geom.frag_mzs.len(), 2);
    assert_eq!(geom.frag_mzs[0], 175.1);
    assert_eq!(geom.frag_mzs[1], 288.2);
    let labels: Vec<String> = geom.frag_labels.iter().map(|l| l.to_string()).collect();
    assert_eq!(labels, vec!["y1".to_string(), "y3^2".to_string()]);
}

/// `id` is echoed back and used as Carafe's map key. A missing `id` NPEs
/// downstream rather than erroring here, so it must not silently default.
#[test]
fn carafe_id_field_is_required() {
    // Build the payload without `id` rather than string-surgering the const,
    // so the test cannot pass because the edit produced malformed JSON.
    let without_id = r#"[
      { "mobility": 0.95, "rt_seconds": 1234.5, "precursor": 650.32,
        "precursor_charge": 2, "precursor_isotopes": [0,1,2],
        "fragments": [175.1, 288.2], "fragment_labels": ["y1","y3^2"] }
    ]"#;
    serde_json::from_str::<serde_json::Value>(without_id).expect("still valid JSON");
    let f = write_targets(without_id);
    assert!(
        read_library_file(f.path()).is_err(),
        "a target without `id` must fail rather than default it"
    );
}

/// Carafe writes `percent` and `absolute`; the CLI's own templates write `pct`
/// and `da`. Dropping either spelling silently breaks one caller.
#[test]
fn carafe_tolerance_spellings_deserialize() {
    let tol: Tolerance =
        serde_json::from_str(CARAFE_TOLERANCES).expect("Carafe's tolerance JSON must deserialize");

    // The `ms` window is `[itol - itol_shift, itol + itol_shift]`: a +-itol
    // window recentred on a measured calibration offset, read as "13 ppm light
    // to 17 ppm heavy". Both edges must stay distinct — collapsing them to a
    // symmetric tolerance would quietly recentre every extraction window.
    assert_eq!(tol.ms, MzTolerance::Ppm((13.0, 17.0)));
    assert_eq!(tol.rt, RtTolerance::Minutes((0.1, 0.1)));
    assert_eq!(tol.mobility, MobilityTolerance::Pct((3.0, 3.0)));
    assert_eq!(tol.quad, QuadTolerance::Absolute((0.1, 0.1)));
}

/// Carafe encodes a systematic calibration offset as
/// `[itol - itol_shift, itol + itol_shift]`, so when the offset exceeds the
/// half-width the LOW edge is negative — a window that is entirely to the
/// heavy side of the target mass. That is a supported input, not a bug, and
/// `mz_range`'s invariant is `low + high >= 0` rather than "both positive".
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

    // Both edges land above the target mass, in ascending order.
    let range = tol.mz_range(1000.0);
    assert!(
        range.start() > 1000.0 && range.end() > range.start(),
        "expected a well-formed window above the target, got {:?}..{:?}",
        range.start(),
        range.end()
    );
}
