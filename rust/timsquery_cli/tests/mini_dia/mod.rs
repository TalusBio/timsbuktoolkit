//! Builder for the mini DIA mzML the end-to-end test reads.
//!
//! The builder is version-controlled, the mzML is not: it is written into
//! `CARGO_TARGET_TMPDIR` on first use and reused afterwards, so a clone
//! generates it once instead of carrying a binary blob in git.

use base64::Engine;
use std::path::{
    Path,
    PathBuf,
};

/// Cycles in the run. Each is one MS1 followed by two DIA windows.
pub const N_CYCLES: usize = 5;
/// Seconds between cycles; the run starts at `RT_START_S`.
pub const RT_STEP_S: f64 = 1.0;
pub const RT_START_S: f64 = 100.0;

/// The one analyte the fixture contains. Charge 2, so its isotopes are one
/// half-neutron apart.
pub const PRECURSOR_MZ: f64 = 650.32;
pub const PRECURSOR_CHARGE: u8 = 2;
pub const ISOTOPE_STEP: f64 = 1.00335 / 2.0;
/// Fragment m/z values, matching the `fragment_labels` the test asks for.
pub const FRAGMENT_MZS: [f64; 2] = [175.119, 288.203];

/// The two DIA isolation windows. The analyte falls in the first.
const WINDOWS: [(f64, f64); 2] = [(650.0, 25.0), (700.0, 25.0)];

/// The peak height at the apex cycle, so the test can assert the extracted
/// values are the ones written here.
pub fn apex_intensity() -> f64 {
    cycle_intensity(N_CYCLES / 2)
}

/// A chromatographic peak apexing at the middle cycle, so extraction has a
/// shape to find rather than a flat line.
fn cycle_intensity(cycle: usize) -> f64 {
    let apex = (N_CYCLES / 2) as f64;
    let d = cycle as f64 - apex;
    1000.0 * (-d * d / 2.0).exp()
}

fn b64_f64(values: &[f64]) -> String {
    let mut bytes = Vec::with_capacity(values.len() * 8);
    for v in values {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    base64::engine::general_purpose::STANDARD.encode(bytes)
}

fn binary_array(values: &[f64], accession: &str, name: &str, unit: &str) -> String {
    let encoded = b64_f64(values);
    format!(
        r#"          <binaryDataArray encodedLength="{len}">
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="{accession}" name="{name}" value="" {unit}/>
            <binary>{encoded}</binary>
          </binaryDataArray>
"#,
        len = encoded.len(),
    )
}

fn spectrum(
    index: usize,
    id: &str,
    ms_level: u8,
    rt_s: f64,
    mzs: &[f64],
    intens: &[f64],
    precursor: Option<(f64, f64)>,
) -> String {
    let mut s = format!(
        r#"      <spectrum index="{index}" id="{id}" defaultArrayLength="{n}">
        <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="{ms_level}"/>
        <cvParam cvRef="MS" accession="{spectrum_type_acc}" name="{spectrum_type}" value=""/>
        <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>
        <cvParam cvRef="MS" accession="MS:1000130" name="positive scan" value=""/>
        <scanList count="1">
          <cvParam cvRef="MS" accession="MS:1000795" name="no combination" value=""/>
          <scan>
            <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="{rt_s}" unitCvRef="UO" unitAccession="UO:0000010" unitName="second"/>
          </scan>
        </scanList>
"#,
        n = mzs.len(),
        spectrum_type_acc = if ms_level == 1 {
            "MS:1000579"
        } else {
            "MS:1000580"
        },
        spectrum_type = if ms_level == 1 {
            "MS1 spectrum"
        } else {
            "MSn spectrum"
        },
    );

    if let Some((target, half_width)) = precursor {
        s.push_str(&format!(
            r#"        <precursorList count="1">
          <precursor>
            <isolationWindow>
              <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{target}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
              <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="{half_width}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
              <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="{half_width}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
            </isolationWindow>
            <selectedIonList count="1">
              <selectedIon>
                <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="{target}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
              </selectedIon>
            </selectedIonList>
            <activation>
              <cvParam cvRef="MS" accession="MS:1000133" name="collision-induced dissociation" value=""/>
            </activation>
          </precursor>
        </precursorList>
"#
        ));
    }

    s.push_str("        <binaryDataArrayList count=\"2\">\n");
    s.push_str(&binary_array(
        mzs,
        "MS:1000514",
        "m/z array",
        r#"unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z""#,
    ));
    s.push_str(&binary_array(
        intens,
        "MS:1000515",
        "intensity array",
        r#"unitCvRef="MS" unitAccession="MS:1000131" unitName="number of detector counts""#,
    ));
    s.push_str("        </binaryDataArrayList>\n      </spectrum>\n");
    s
}

fn build() -> String {
    let mut spectra = String::new();
    let mut index = 0usize;

    for cycle in 0..N_CYCLES {
        let rt = RT_START_S + cycle as f64 * RT_STEP_S;
        let inten = cycle_intensity(cycle);

        // MS1: the analyte's first three isotopes.
        let ms1_mzs: Vec<f64> = (0..3)
            .map(|i| PRECURSOR_MZ + i as f64 * ISOTOPE_STEP)
            .collect();
        let ms1_intens = vec![inten, inten * 0.5, inten * 0.2];
        spectra.push_str(&spectrum(
            index,
            &format!("scan={}", index + 1),
            1,
            rt,
            &ms1_mzs,
            &ms1_intens,
            None,
        ));
        index += 1;

        for (w, (target, half_width)) in WINDOWS.iter().enumerate() {
            // Only the window containing the precursor carries its fragments;
            // the other is present so the run has a real DIA cycle structure.
            let (mzs, intens): (Vec<f64>, Vec<f64>) = if w == 0 {
                (FRAGMENT_MZS.to_vec(), vec![inten * 0.8, inten * 0.4])
            } else {
                (vec![300.0], vec![10.0])
            };
            spectra.push_str(&spectrum(
                index,
                &format!("scan={}", index + 1),
                2,
                rt,
                &mzs,
                &intens,
                Some((*target, *half_width)),
            ));
            index += 1;
        }
    }

    format!(
        r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd" id="mini_dia" version="1.1.0">
  <cvList count="2">
    <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" version="4.1.92" URI="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"/>
    <cv id="UO" fullName="Unit Ontology" version="09:04:2014" URI="https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo"/>
  </cvList>
  <fileDescription>
    <fileContent>
      <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>
      <cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" value=""/>
    </fileContent>
  </fileDescription>
  <softwareList count="1">
    <software id="tbtk" version="0.0.0">
      <cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool" value=""/>
    </software>
  </softwareList>
  <instrumentConfigurationList count="1">
    <instrumentConfiguration id="IC1">
      <cvParam cvRef="MS" accession="MS:1000031" name="instrument model" value=""/>
    </instrumentConfiguration>
  </instrumentConfigurationList>
  <dataProcessingList count="1">
    <dataProcessing id="DP1">
      <processingMethod order="0" softwareRef="tbtk">
        <cvParam cvRef="MS" accession="MS:1000544" name="Conversion to mzML" value=""/>
      </processingMethod>
    </dataProcessing>
  </dataProcessingList>
  <run id="mini_dia" defaultInstrumentConfigurationRef="IC1">
    <spectrumList count="{count}" defaultDataProcessingRef="DP1">
{spectra}    </spectrumList>
  </run>
</mzML>
"#,
        count = index,
    )
}

/// Path to the fixture, writing it if this clone has not built it yet.
pub fn path() -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    std::fs::create_dir_all(dir).expect("create target tmpdir");
    let path = dir.join("mini_dia.mzML");
    if !path.exists() {
        std::fs::write(&path, build()).expect("write mini DIA mzML");
    }
    path
}
