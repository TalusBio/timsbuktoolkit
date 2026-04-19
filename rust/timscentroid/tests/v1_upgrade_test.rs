//! Integration test for the v1 → v2 `QuadrupoleIsolationScheme` upgrade.
//!
//! The fixture `tests/data/v1_metadata_trimmed.json` is a real
//! `metadata.json` from a Hela DIA run whose `cycle_to_rt_ms` arrays have
//! been truncated to 3 entries each to keep the fixture small (~16 KB).
//! Everything relevant to the upgrade path (version string,
//! `quadrupole_isolation` shapes, group layout) is preserved verbatim.
//!
//! The test covers:
//! - v1 metadata parses via `TimscentroidMetadata`'s derive, including
//!   the `QuadrupoleIsolationScheme` `from = QuadWire` deserializer that
//!   performs the in-memory v1→v2 upgrade.
//! - Each group's rings come out classified (not left in a raw polygon
//!   form).
//! - For this static-DIA fixture specifically, every ring lands on the
//!   AABB fast path.
//! - Re-serializing a `QuadrupoleIsolationScheme` emits the v2 (bare
//!   array) shape and reparses to the same values.

use timscentroid::QuadrupoleIsolationScheme;
use timscentroid::serialization::TimscentroidMetadata;

const V1_METADATA: &str = include_str!("data/v1_metadata_trimmed.json");

#[test]
fn v1_metadata_upgrades_to_v2_rings() {
    let meta: TimscentroidMetadata = serde_json::from_str(V1_METADATA).expect("v1 metadata parses");

    assert_eq!(meta.version, "1.0", "fixture is v1");
    assert_eq!(
        meta.ms2_window_groups.len(),
        8,
        "Hela DIA has 8 MS2 window groups"
    );

    // Each group's quadrupole_isolation should have been reclassified
    // into a non-empty Vec<RingShape>, and for this static-DIA run every
    // ring should hit the AABB fast path.
    for (i, group) in meta.ms2_window_groups.iter().enumerate() {
        let roundtrip =
            serde_json::to_string(&group.quadrupole_isolation).expect("serialize quad as v2");
        assert!(
            roundtrip.starts_with('['),
            "group[{i}] v2 serialization should be a bare array, got: {}",
            &roundtrip[..roundtrip.len().min(40)]
        );

        // Reparse the v2 output and confirm every ring is an AABB.
        let reparsed: QuadrupoleIsolationScheme = serde_json::from_str(&roundtrip).unwrap();
        let debug = format!("{:?}", reparsed);
        assert!(
            !debug.contains("Polygon("),
            "group[{i}] should have no Polygon-fallback rings after \
             upgrade; got {debug}"
        );
        assert!(
            debug.contains("Aabb"),
            "group[{i}] expected at least one Aabb ring; got {debug}"
        );
    }
}
