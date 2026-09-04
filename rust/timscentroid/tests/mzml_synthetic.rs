//! CI coverage for the mzML ingest happy path, without a checked-in binary
//! fixture: build a tiny centroided DIA mzML in-process (2 cycles × 2 windows,
//! no ion mobility), then run it through `from_mzml_file` and assert the derived
//! structure (cycles, windows, mobility kind, RT).

#![cfg(feature = "mzdata")]

use base64::Engine;
use timscentroid::MobilityKind;
use timscentroid::reader::mzdata::from_mzml_file;

/// base64 of a little-endian f64 array (mzML `64-bit float`, no compression).
fn b64_f64(vals: &[f64]) -> String {
    let mut bytes = Vec::with_capacity(vals.len() * 8);
    for v in vals {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    base64::engine::general_purpose::STANDARD.encode(bytes)
}

fn binary_array(accession: &str, name: &str, vals: &[f64]) -> String {
    let enc = b64_f64(vals);
    format!(
        r#"<binaryDataArray encodedLength="{len}">
          <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>
          <cvParam cvRef="MS" accession="MS:1000576" name="no compression"/>
          <cvParam cvRef="MS" accession="{accession}" name="{name}"/>
          <binary>{enc}</binary>
        </binaryDataArray>"#,
        len = enc.len(),
    )
}

/// One `<spectrum>`. `iso` = Some((lower_mz, upper_mz)) for an MS2 window.
fn spectrum(
    index: usize,
    id: &str,
    ms_level: u8,
    rt_min: f64,
    iso: Option<(f64, f64)>,
    peaks: &[(f64, f64)],
) -> String {
    let mzs: Vec<f64> = peaks.iter().map(|(m, _)| *m).collect();
    let ints: Vec<f64> = peaks.iter().map(|(_, i)| *i).collect();
    let precursor = match iso {
        Some((lo, hi)) => {
            let target = (lo + hi) / 2.0;
            let half = (hi - lo) / 2.0;
            format!(
                r#"<precursorList count="1"><precursor>
              <isolationWindow>
                <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{target}"/>
                <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="{half}"/>
                <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="{half}"/>
              </isolationWindow>
            </precursor></precursorList>"#
            )
        }
        None => String::new(),
    };
    format!(
        r#"<spectrum index="{index}" id="{id}" defaultArrayLength="{n}">
      <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="{ms_level}"/>
      <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum"/>
      <scanList count="1"><scan>
        <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="{rt_min}" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>
      </scan></scanList>
      {precursor}
      <binaryDataArrayList count="2">
        {mz_arr}
        {int_arr}
      </binaryDataArrayList>
    </spectrum>"#,
        n = peaks.len(),
        mz_arr = binary_array("MS:1000514", "m/z array", &mzs),
        int_arr = binary_array("MS:1000515", "intensity array", &ints),
    )
}

fn synthetic_dia_mzml() -> String {
    // 2 cycles; each = MS1 then two MS2 windows [400,402] and [402,404].
    let mut specs = String::new();
    let mut idx = 0;
    for (cycle, rt) in [(0usize, 0.10_f64), (1, 0.13)] {
        specs.push_str(&spectrum(
            idx,
            &format!("scan=cyc{cycle}_ms1"),
            1,
            rt,
            None,
            &[(500.0, 1000.0), (600.0, 800.0)],
        ));
        idx += 1;
        for (lo, hi) in [(400.0, 402.0), (402.0, 404.0)] {
            specs.push_str(&spectrum(
                idx,
                &format!("scan=cyc{cycle}_w{lo}"),
                2,
                rt,
                Some((lo, hi)),
                &[(300.0, 500.0), (350.0, 250.0)],
            ));
            idx += 1;
        }
    }
    format!(
        r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">
  <run id="synthetic">
    <spectrumList count="{count}">
      {specs}
    </spectrumList>
  </run>
</mzML>"#,
        count = idx,
    )
}

#[test]
fn ingests_synthetic_dia_mzml() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("synthetic.mzML");
    std::fs::write(&path, synthetic_dia_mzml()).unwrap();

    let cfg = timscentroid::IndexingCentroidingConfig::default();
    let idx = from_mzml_file(&path, &cfg).expect("synthetic DIA mzML should ingest");

    // No ion-mobility axis declared → Absent.
    assert_eq!(idx.mobility_kind(), &MobilityKind::Absent);
    // Two distinct isolation windows, two MS1 cycles.
    assert_eq!(idx.n_ms2_window_groups(), 2, "expected 2 windows");
    assert_eq!(idx.ms1_cycle_mapping().len(), 2, "expected 2 cycles");
    // RT stored in ms, monotonic (0.10 min → 6000 ms, 0.13 min → 7800 ms).
    let rts = idx.ms1_cycle_mapping().clone().unpack();
    assert!(
        rts.windows(2).all(|w| w[0] < w[1]),
        "cycle RTs must increase: {rts:?}"
    );
    assert_eq!(rts[0], 6000);
}
