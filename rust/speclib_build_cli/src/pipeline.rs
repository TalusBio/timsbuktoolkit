use std::io::{
    BufRead,
    Write,
};

use indicatif::{
    ProgressBar,
    ProgressFinish,
    ProgressStyle,
};
use timsseek::ProteinSlice;
use timsseek::data_sources::speclib::SpeclibWriter;
use timsseek::digest::digestion::{
    DigestionEnd,
    DigestionParameters,
    DigestionPattern,
};
use timsseek::protein::fasta::ProteinSequenceCollection;

use crate::config::SpeclibBuildConfig;
use crate::decoys::{
    DecoyMode,
    generate_decoy,
};
use crate::dedup::PeptideDedup;
use crate::entry::{
    EntryFilters,
    build_entry,
};
use crate::koina::KoinaClient;
use crate::koina::models::{
    FragmentModel,
    PredictionInput,
    RtModel,
};
use crate::mods::{
    Modification,
    apply_fixed_mods,
    expand_variable_mods,
};

// ── BatchItem ───────────────────────────────────────────────────────────────

struct BatchItem {
    sequence: String,
    charge: u8,
    decoy: bool,
    decoy_group: u32,
}

// ── flush_batch ─────────────────────────────────────────────────────────────

async fn flush_batch(
    koina: &KoinaClient,
    writer: &mut SpeclibWriter<std::fs::File>,
    batch: &mut Vec<BatchItem>,
    filters: &EntryFilters,
    nce: f32,
    entry_id: &mut u32,
    progress: &ProgressBar,
) -> Result<(), Box<dyn std::error::Error>> {
    if batch.is_empty() {
        return Ok(());
    }

    // Koina/Prosit accepts UNIMOD notation — convert our short form [U:N] → [UNIMOD:N]
    let inputs: Vec<PredictionInput> = batch
        .iter()
        .map(|item| PredictionInput {
            sequence: item.sequence.replace("[U:", "[UNIMOD:"),
            charge: item.charge,
            nce,
        })
        .collect();

    let (fragments, rts) = tokio::try_join!(
        async {
            koina
                .predict_fragments(&inputs)
                .await
                .map_err(|e| -> Box<dyn std::error::Error> { e.into() })
        },
        async {
            koina
                .predict_rt(&inputs)
                .await
                .map_err(|e| -> Box<dyn std::error::Error> { e.into() })
        },
    )?;

    for ((item, fragment), rt) in batch.drain(..).zip(fragments.iter()).zip(rts.iter()) {
        if let Some(elem) = build_entry(
            &item.sequence,
            item.charge,
            item.decoy,
            item.decoy_group,
            *entry_id,
            fragment,
            rt,
            filters,
        ) {
            writer
                .append(&elem)
                .map_err(|e| format!("write error: {e:?}"))?;
            *entry_id += 1;
        }
        progress.inc(1);
    }

    Ok(())
}

// ── run ─────────────────────────────────────────────────────────────────────

pub async fn run(config: &SpeclibBuildConfig) -> Result<(), Box<dyn std::error::Error>> {
    // ── Phase 1: Get base peptides ──────────────────────────────────────────

    let base_peptides: Vec<ProteinSlice> = if let Some(fasta_uri) = &config.fasta {
        tracing::info!("Reading FASTA: {}", fasta_uri);
        let collection = if is_remote_uri(fasta_uri) {
            // Download remote FASTA to a tempfile (preserving extension for any
            // downstream sniffing), then parse as usual.
            let mut reader = tims_stage::open_reader(fasta_uri)
                .map_err(|e| format!("open fasta {fasta_uri}: {e}"))?;
            let ext = std::path::Path::new(fasta_uri.trim_end_matches('/'))
                .extension()
                .and_then(|s| s.to_str())
                .unwrap_or("fasta");
            let mut tf = tempfile::Builder::new()
                .prefix("speclib-fasta-")
                .suffix(&format!(".{ext}"))
                .tempfile()?;
            std::io::copy(&mut reader, tf.as_file_mut())?;
            tf.as_file_mut().flush()?;
            ProteinSequenceCollection::from_fasta_file(tf.path())?
        } else {
            ProteinSequenceCollection::from_fasta_file(std::path::Path::new(fasta_uri))?
        };
        let total_aa: usize = collection.sequences.iter().map(|p| p.sequence.len()).sum();
        tracing::info!(
            "Read {} proteins ({} amino acids)",
            collection.sequences.len(),
            total_aa,
        );

        let params = DigestionParameters {
            min_length: config.digestion.min_length,
            max_length: config.digestion.max_length,
            pattern: DigestionPattern::trypsin(),
            digestion_end: DigestionEnd::CTerm,
            max_missed_cleavages: config.digestion.missed_cleavages,
        };

        // Warn about proteins containing non-standard amino acids
        for prot in &collection.sequences {
            if prot
                .sequence
                .chars()
                .any(|c| matches!(c, 'U' | 'B' | 'J' | 'Z' | 'X'))
            {
                tracing::warn!(
                    "Protein {} contains non-standard amino acids (U/B/J/Z/X) — affected peptides will be skipped",
                    prot.description,
                );
            }
        }

        let protein_seqs: Vec<std::sync::Arc<str>> = collection
            .sequences
            .iter()
            .map(|p| p.sequence.clone())
            .collect();
        let raw = params.digest_multiple(&protein_seqs);
        tracing::info!("Digested into {} raw peptides", raw.len());

        let estimated =
            PeptideDedup::estimate_from_proteins(total_aa, config.digestion.missed_cleavages);
        let deduped = PeptideDedup::dedup(raw, estimated);
        let (deduped, skipped) = filter_nonstandard_aa(deduped);
        if skipped > 0 {
            tracing::warn!("Skipped {skipped} peptides containing non-standard amino acids");
        }
        tracing::info!("Deduplicated to {} unique peptides", deduped.len());
        deduped
    } else if let Some(list_uri) = &config.peptide_list {
        tracing::info!("Reading peptide list: {}", list_uri);
        let raw_reader = tims_stage::open_reader(list_uri)
            .map_err(|e| format!("open peptide list {list_uri}: {e}"))?;
        let reader = std::io::BufReader::new(raw_reader);
        let mut slices: Vec<ProteinSlice> = Vec::new();
        for (i, line) in reader.lines().enumerate() {
            let line = line?;
            let seq = line.trim().to_string();
            if seq.is_empty() || seq.starts_with('#') {
                continue;
            }
            slices.push(ProteinSlice::from_string(seq, false, i as u32));
        }
        tracing::info!("Read {} peptides from list", slices.len());

        let estimated = slices.len();
        let deduped = PeptideDedup::dedup(slices, estimated);
        let (deduped, skipped) = filter_nonstandard_aa(deduped);
        if skipped > 0 {
            tracing::warn!("Skipped {skipped} peptides containing non-standard amino acids");
        }
        tracing::info!("Deduplicated to {} unique peptides", deduped.len());
        deduped
    } else {
        return Err("No input provided (--fasta or --peptide-list)".into());
    };

    if base_peptides.is_empty() {
        return Err("No peptides after digestion/dedup — nothing to do".into());
    }

    // ── Phase 2: Set up prediction ──────────────────────────────────────────

    let fragment_model = FragmentModel::from_name(&config.prediction.fragment_model)?;
    let rt_model = RtModel::from_name(&config.prediction.rt_model)?;
    let koina = KoinaClient::new(
        config.prediction.koina_url.clone(),
        fragment_model,
        rt_model,
    );

    let fixed_mods: Vec<Modification> = config
        .modifications
        .fixed
        .iter()
        .map(|s| Modification::parse(s))
        .collect::<Result<Vec<_>, _>>()?;

    let var_mods: Vec<Modification> = config
        .modifications
        .variable
        .iter()
        .map(|s| Modification::parse(s))
        .collect::<Result<Vec<_>, _>>()?;

    let decoy_mode = DecoyMode::from_str(&config.decoys.strategy)?;

    let filters = EntryFilters {
        max_ions: config.filters.max_ions,
        min_mz: config.filters.min_mz,
        max_mz: config.filters.max_mz,
        min_ion_mz: config.filters.min_ion_mz,
        max_ion_mz: config.filters.max_ion_mz,
        min_ions: config.filters.min_ions,
    };

    let nce = config.prediction.nce;
    let batch_size = config.prediction.batch_size;
    let request_delay = if config.prediction.request_delay_ms > 0 {
        Some(std::time::Duration::from_millis(
            config.prediction.request_delay_ms,
        ))
    } else {
        None
    };
    let charges: Vec<u8> = (config.charges.min..=config.charges.max).collect();
    let max_var_mods = config.modifications.max_variable;

    // ── Phase 3: Streaming expansion + prediction ───────────────────────────

    let output_uri = &config.output;
    let remote_output = is_remote_uri(output_uri);
    // Keep the tempfile handle alive until after upload so the file isn't
    // dropped out from under us.
    let output_tempfile: Option<tempfile::NamedTempFile> = if remote_output {
        let ext = std::path::Path::new(output_uri.trim_end_matches('/'))
            .extension()
            .and_then(|s| s.to_str())
            .unwrap_or("msgpack.zst");
        let tf = tempfile::Builder::new()
            .prefix("speclib-out-")
            .suffix(&format!(".{ext}"))
            .tempfile()?;
        Some(tf)
    } else {
        // Ensure parent directory exists for local output.
        if let Some(parent) = std::path::Path::new(output_uri).parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent)?;
            }
        }
        None
    };
    let working_path: std::path::PathBuf = if let Some(tf) = &output_tempfile {
        tf.path().to_path_buf()
    } else {
        std::path::PathBuf::from(output_uri)
    };

    let out_file = std::fs::File::create(&working_path)?;
    let mut writer = SpeclibWriter::new_msgpack_zstd(out_file)?;

    // Estimate total items for progress bar: peptides * mod_variants * charges * (1 + decoy)
    let decoy_mult: u64 = if decoy_mode != DecoyMode::None { 2 } else { 1 };
    let estimated_total = base_peptides.len() as u64 * charges.len() as u64 * decoy_mult;
    let progress = ProgressBar::new(estimated_total).with_style(
        ProgressStyle::with_template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({per_sec}, ETA {eta})",
        )
        .unwrap()
        .progress_chars("#>-"),

        ).with_finish(ProgressFinish::AndLeave);

    let mut batch: Vec<BatchItem> = Vec::with_capacity(batch_size);
    let mut entry_id: u32 = 0;
    let mut decoy_group: u32 = 0;

    for digest_slice in &base_peptides {
        let base_seq = digest_slice.as_str();

        // Apply fixed modifications.
        let fixed_seq = if fixed_mods.is_empty() {
            base_seq.to_string()
        } else {
            apply_fixed_mods(base_seq, &fixed_mods)
        };

        // Expand variable modifications.
        let mod_variants = if var_mods.is_empty() {
            vec![fixed_seq]
        } else {
            expand_variable_mods(&fixed_seq, &var_mods, max_var_mods)
        };

        for variant in &mod_variants {
            for &charge in &charges {
                // Target
                batch.push(BatchItem {
                    sequence: variant.clone(),
                    charge,
                    decoy: false,
                    decoy_group,
                });

                // Decoy
                if decoy_mode != DecoyMode::None {
                    let decoy_seq = generate_decoy(variant, decoy_mode);
                    batch.push(BatchItem {
                        sequence: decoy_seq,
                        charge,
                        decoy: true,
                        decoy_group,
                    });
                }

                if batch.len() >= batch_size {
                    flush_batch(
                        &koina,
                        &mut writer,
                        &mut batch,
                        &filters,
                        nce,
                        &mut entry_id,
                        &progress,
                    )
                    .await?;
                    if let Some(delay) = request_delay {
                        tokio::time::sleep(delay).await;
                    }
                }
            }
        }

        decoy_group += 1;
    }

    // Flush remaining items.
    flush_batch(
        &koina,
        &mut writer,
        &mut batch,
        &filters,
        nce,
        &mut entry_id,
        &progress,
    )
    .await?;

    progress.finish_with_message("done");
    writer.finish()?;

    // Upload tempfile for remote outputs. Hold tempfile alive until upload
    // completes, then let it drop (which removes the local file).
    if remote_output {
        tims_stage::upload_file(&working_path, output_uri)
            .map_err(|e| format!("upload {output_uri}: {e}"))?;
    }
    drop(output_tempfile);

    tracing::info!("Wrote {} entries to {}", entry_id, output_uri);

    Ok(())
}

use tims_stage::is_remote_uri;

const NONSTANDARD_AA: &[char] = &['U', 'B', 'J', 'Z', 'X'];

/// Filter out peptides containing non-standard amino acids (U, B, J, Z, X).
/// Returns (kept, skipped_count).
fn filter_nonstandard_aa(peptides: Vec<ProteinSlice>) -> (Vec<ProteinSlice>, usize) {
    let before = peptides.len();
    let kept: Vec<ProteinSlice> = peptides
        .into_iter()
        .filter(|p| !p.as_str().chars().any(|c| NONSTANDARD_AA.contains(&c)))
        .collect();
    let skipped = before - kept.len();
    (kept, skipped)
}
