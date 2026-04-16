use crate::koina::models::{
    FragmentPrediction, KoinaRequest, KoinaTensor, KoinaTensorData, KoinaResponse, PredictionInput,
    RtPrediction,
};

// ── Request builders ─────────────────────────────────────────────────────────

/// Build a fragment-intensity inference request.
///
/// Tensors sent:
/// - `peptide_sequences` (BYTES)
/// - `precursor_charges` (INT32)
/// - `collision_energies` (FP32)
pub fn build_fragment_request(inputs: &[PredictionInput], request_id: &str) -> KoinaRequest {
    let n = inputs.len();
    let sequences: Vec<String> = inputs.iter().map(|i| i.sequence.clone()).collect();
    let charges: Vec<i32> = inputs.iter().map(|i| i.charge as i32).collect();
    let nces: Vec<f32> = inputs.iter().map(|i| i.nce).collect();

    KoinaRequest {
        id: request_id.to_string(),
        inputs: vec![
            KoinaTensor {
                name: "peptide_sequences".to_string(),
                shape: vec![n, 1],
                datatype: "BYTES".to_string(),
                data: KoinaTensorData::Strings(sequences),
            },
            KoinaTensor {
                name: "precursor_charges".to_string(),
                shape: vec![n, 1],
                datatype: "INT32".to_string(),
                data: KoinaTensorData::Ints(charges),
            },
            KoinaTensor {
                name: "collision_energies".to_string(),
                shape: vec![n, 1],
                datatype: "FP32".to_string(),
                data: KoinaTensorData::Floats(nces),
            },
        ],
    }
}

/// Build an RT inference request.
///
/// Tensor sent:
/// - `peptide_sequences` (BYTES)
pub fn build_rt_request(inputs: &[PredictionInput], request_id: &str) -> KoinaRequest {
    let n = inputs.len();
    let sequences: Vec<String> = inputs.iter().map(|i| i.sequence.clone()).collect();

    KoinaRequest {
        id: request_id.to_string(),
        inputs: vec![KoinaTensor {
            name: "peptide_sequences".to_string(),
            shape: vec![n, 1],
            datatype: "BYTES".to_string(),
            data: KoinaTensorData::Strings(sequences),
        }],
    }
}

// ── Response parsers ─────────────────────────────────────────────────────────

/// Parse a fragment-intensity Koina response into per-peptide predictions.
///
/// Expects outputs named `annotation`, `mz`, and `intensities`.
/// Zero-intensity ions are filtered out from the returned predictions.
pub fn parse_fragment_response(
    response: KoinaResponse,
    batch_size: usize,
) -> Result<Vec<FragmentPrediction>, String> {
    let annotations_tensor = find_output(&response, "annotation")?;
    let mz_tensor = find_output(&response, "mz")?;
    let intensities_tensor = find_output(&response, "intensities")?;

    // Each tensor has shape [batch_size, num_ions]; flatten then chunk.
    let ions_per_peptide = total_len(annotations_tensor) / batch_size;

    // Koina returns annotations like "y1+1" (plus for charge), but
    // IonAnnot expects "y1^1" (caret). Convert.
    let annotations_flat: Vec<String> = annotations_tensor
        .data
        .iter()
        .map(|v| {
            let s = v.as_str().unwrap_or("");
            koina_annotation_to_mzpaf(s)
        })
        .collect();
    let mzs_flat: Vec<f64> = mz_tensor
        .data
        .iter()
        .map(|v| v.as_f64().unwrap_or(0.0))
        .collect();
    let intensities_flat: Vec<f32> = intensities_tensor
        .data
        .iter()
        .map(|v| v.as_f64().unwrap_or(0.0) as f32)
        .collect();

    let mut predictions = Vec::with_capacity(batch_size);
    for i in 0..batch_size {
        let start = i * ions_per_peptide;
        let end = start + ions_per_peptide;

        let mut anns = Vec::new();
        let mut mzs = Vec::new();
        let mut ints = Vec::new();

        for j in start..end {
            let intensity = intensities_flat[j];
            if intensity > 0.0 {
                anns.push(annotations_flat[j].clone());
                mzs.push(mzs_flat[j]);
                ints.push(intensity);
            }
        }

        predictions.push(FragmentPrediction {
            annotations: anns,
            mzs,
            intensities: ints,
        });
    }

    Ok(predictions)
}

/// Parse an RT Koina response into per-peptide iRT predictions.
///
/// Expects an output named `irt`.
pub fn parse_rt_response(
    response: KoinaResponse,
    batch_size: usize,
) -> Result<Vec<RtPrediction>, String> {
    let irt_tensor = find_output(&response, "irt")?;

    if irt_tensor.data.len() != batch_size {
        return Err(format!(
            "irt tensor has {} entries, expected {batch_size}",
            irt_tensor.data.len()
        ));
    }

    let predictions = irt_tensor
        .data
        .iter()
        .map(|v| RtPrediction {
            irt: v.as_f64().unwrap_or(0.0) as f32,
        })
        .collect();

    Ok(predictions)
}

// ── Helpers ──────────────────────────────────────────────────────────────────

fn find_output<'a>(
    response: &'a KoinaResponse,
    name: &str,
) -> Result<&'a crate::koina::models::KoinaOutputTensor, String> {
    response
        .outputs
        .iter()
        .find(|o| o.name == name)
        .ok_or_else(|| format!("Koina response missing output tensor {name:?}"))
}

fn total_len(tensor: &crate::koina::models::KoinaOutputTensor) -> usize {
    tensor.data.len()
}

/// Convert Koina annotation format to mzPAF format.
///
/// Koina: "y1+1" (series, ordinal, '+', charge)
/// mzPAF: "y1^1" (series, ordinal, '^', charge)
fn koina_annotation_to_mzpaf(s: &str) -> String {
    // Find the last '+' that separates ordinal from charge
    // Annotations look like: y1+1, y1+2, y1+3, b12+2
    if let Some(pos) = s.rfind('+') {
        let (prefix, charge) = s.split_at(pos);
        format!("{prefix}^{}", &charge[1..])
    } else {
        s.to_string()
    }
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::koina::models::{KoinaOutputTensor, KoinaResponse, KoinaTensorData, PredictionInput};

    fn sample_inputs() -> Vec<PredictionInput> {
        vec![
            PredictionInput {
                sequence: "PEPTIDE".to_string(),
                charge: 2,
                nce: 0.3,
            },
            PredictionInput {
                sequence: "LAGER".to_string(),
                charge: 3,
                nce: 0.25,
            },
        ]
    }

    #[test]
    fn test_build_fragment_request() {
        let inputs = sample_inputs();
        let req = build_fragment_request(&inputs, "req-1");

        assert_eq!(req.id, "req-1");
        assert_eq!(req.inputs.len(), 3);

        let seq_tensor = &req.inputs[0];
        assert_eq!(seq_tensor.name, "peptide_sequences");
        assert_eq!(seq_tensor.datatype, "BYTES");
        assert_eq!(seq_tensor.shape, vec![2, 1]);
        if let KoinaTensorData::Strings(ref seqs) = seq_tensor.data {
            assert_eq!(seqs, &["PEPTIDE", "LAGER"]);
        } else {
            panic!("expected Strings data");
        }

        let charge_tensor = &req.inputs[1];
        assert_eq!(charge_tensor.name, "precursor_charges");
        assert_eq!(charge_tensor.datatype, "INT32");
        assert_eq!(charge_tensor.shape, vec![2, 1]);
        if let KoinaTensorData::Ints(ref charges) = charge_tensor.data {
            assert_eq!(charges, &[2, 3]);
        } else {
            panic!("expected Ints data");
        }

        let nce_tensor = &req.inputs[2];
        assert_eq!(nce_tensor.name, "collision_energies");
        assert_eq!(nce_tensor.datatype, "FP32");
        assert_eq!(nce_tensor.shape, vec![2, 1]);
        if let KoinaTensorData::Floats(ref nces) = nce_tensor.data {
            assert!((nces[0] - 0.3).abs() < 1e-6);
            assert!((nces[1] - 0.25).abs() < 1e-6);
        } else {
            panic!("expected Floats data");
        }
    }

    #[test]
    fn test_build_rt_request() {
        let inputs = sample_inputs();
        let req = build_rt_request(&inputs, "req-rt-1");

        assert_eq!(req.id, "req-rt-1");
        assert_eq!(req.inputs.len(), 1);

        let seq_tensor = &req.inputs[0];
        assert_eq!(seq_tensor.name, "peptide_sequences");
        assert_eq!(seq_tensor.datatype, "BYTES");
        assert_eq!(seq_tensor.shape, vec![2, 1]);
        if let KoinaTensorData::Strings(ref seqs) = seq_tensor.data {
            assert_eq!(seqs, &["PEPTIDE", "LAGER"]);
        } else {
            panic!("expected Strings data");
        }
    }

    #[test]
    fn test_parse_fragment_response_filters_zeros() {
        // 2 peptides, 3 ions each (flat: 6 entries per tensor)
        // Peptide 0: ions [1.0, 0.0, 0.5] → only 2 survive
        // Peptide 1: ions [0.0, 0.0, 0.8] → only 1 survives
        let annotations_data: Vec<serde_json::Value> = vec![
            "y1^1", "b2^1", "y3^1", "y1^1", "b2^1", "y3^1",
        ]
        .into_iter()
        .map(|s| serde_json::Value::String(s.to_string()))
        .collect();

        let mz_data: Vec<serde_json::Value> = vec![100.0_f64, 200.0, 300.0, 110.0, 210.0, 310.0]
            .into_iter()
            .map(serde_json::Value::from)
            .collect();

        let intensity_data: Vec<serde_json::Value> =
            vec![1.0_f64, 0.0, 0.5, 0.0, 0.0, 0.8]
                .into_iter()
                .map(serde_json::Value::from)
                .collect();

        let response = KoinaResponse {
            id: "test".to_string(),
            outputs: vec![
                KoinaOutputTensor {
                    name: "annotation".to_string(),
                    datatype: "BYTES".to_string(),
                    shape: vec![2, 3],
                    data: annotations_data,
                },
                KoinaOutputTensor {
                    name: "mz".to_string(),
                    datatype: "FP32".to_string(),
                    shape: vec![2, 3],
                    data: mz_data,
                },
                KoinaOutputTensor {
                    name: "intensities".to_string(),
                    datatype: "FP32".to_string(),
                    shape: vec![2, 3],
                    data: intensity_data,
                },
            ],
        };

        let preds = parse_fragment_response(response, 2).unwrap();
        assert_eq!(preds.len(), 2);

        // Peptide 0: ions at index 0 (1.0) and 2 (0.5) survive; index 1 (0.0) is dropped.
        assert_eq!(preds[0].intensities.len(), 2);
        assert_eq!(preds[0].annotations, vec!["y1^1", "y3^1"]);
        assert!((preds[0].mzs[0] - 100.0).abs() < 1e-6);
        assert!((preds[0].mzs[1] - 300.0).abs() < 1e-6);

        // Peptide 1: only index 2 (0.8) survives.
        assert_eq!(preds[1].intensities.len(), 1);
        assert_eq!(preds[1].annotations, vec!["y3^1"]);
        assert!((preds[1].mzs[0] - 310.0).abs() < 1e-6);
    }
}
