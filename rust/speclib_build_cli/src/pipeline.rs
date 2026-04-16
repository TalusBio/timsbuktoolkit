use std::io::BufRead;

use indicatif::{ProgressBar, ProgressStyle};
use timsseek::DigestSlice;
use timsseek::data_sources::speclib::SpeclibWriter;
use timsseek::digest::digestion::{DigestionEnd, DigestionParameters, DigestionPattern};
use timsseek::protein::fasta::ProteinSequenceCollection;

use crate::config::SpeclibBuildConfig;
use crate::decoys::{DecoyMode, generate_decoy};
use crate::dedup::PeptideDedup;
use crate::entry::{EntryFilters, build_entry, strip_mods};
use crate::koina::KoinaClient;
use crate::koina::models::{FragmentModel, PredictionInput, RtModel};
use crate::mods::{Modification, apply_fixed_mods, expand_variable_mods};

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

    // Koina (Prosit) expects bare AA sequences — strip bracket mods
    let inputs: Vec<PredictionInput> = batch
        .iter()
        .map(|item| PredictionInput {
            sequence: strip_mods(&item.sequence),
            charge: item.charge,
            nce,
        })
        .collect();

    let (fragments, rts) = tokio::try_join!(
        async { koina.predict_fragments(&inputs).await.map_err(|e| -> Box<dyn std::error::Error> { e.into() }) },
        async { koina.predict_rt(&inputs).await.map_err(|e| -> Box<dyn std::error::Error> { e.into() }) },
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
            writer.append(&elem).map_err(|e| format!("write error: {e:?}"))?;
            *entry_id += 1;
        }
        progress.inc(1);
    }

    Ok(())
}

// ── run ─────────────────────────────────────────────────────────────────────

pub async fn run(config: &SpeclibBuildConfig) -> Result<(), Box<dyn std::error::Error>> {
    // ── Phase 1: Get base peptides ──────────────────────────────────────────

    let base_peptides: Vec<DigestSlice> = if let Some(fasta_path) = &config.fasta {
        tracing::info!("Reading FASTA: {}", fasta_path.display());
        let collection = ProteinSequenceCollection::from_fasta_file(fasta_path)?;
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

        let protein_seqs: Vec<std::sync::Arc<str>> = collection
            .sequences
            .iter()
            .map(|p| p.sequence.clone())
            .collect();
        let raw = params.digest_multiple(&protein_seqs);
        tracing::info!("Digested into {} raw peptides", raw.len());

        let estimated = PeptideDedup::estimate_from_proteins(total_aa, config.digestion.missed_cleavages);
        let deduped = PeptideDedup::dedup(raw, estimated);
        tracing::info!("Deduplicated to {} unique peptides", deduped.len());
        deduped
    } else if let Some(list_path) = &config.peptide_list {
        tracing::info!("Reading peptide list: {}", list_path.display());
        let file = std::fs::File::open(list_path)?;
        let reader = std::io::BufReader::new(file);
        let mut slices: Vec<DigestSlice> = Vec::new();
        for (i, line) in reader.lines().enumerate() {
            let line = line?;
            let seq = line.trim().to_string();
            if seq.is_empty() || seq.starts_with('#') {
                continue;
            }
            slices.push(DigestSlice::from_string(seq, false, i as u32));
        }
        tracing::info!("Read {} peptides from list", slices.len());

        let estimated = slices.len();
        let deduped = PeptideDedup::dedup(slices, estimated);
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
        Some(std::time::Duration::from_millis(config.prediction.request_delay_ms))
    } else {
        None
    };
    let charges: Vec<u8> = (config.charges.min..=config.charges.max).collect();
    let max_var_mods = config.modifications.max_variable;

    // ── Phase 3: Streaming expansion + prediction ───────────────────────────

    let out_file = std::fs::File::create(&config.output)?;
    let mut writer = SpeclibWriter::new_msgpack_zstd(out_file)?;

    // Estimate total items for progress bar: peptides * mod_variants * charges * (1 + decoy)
    let decoy_mult: u64 = if decoy_mode != DecoyMode::None { 2 } else { 1 };
    let estimated_total =
        base_peptides.len() as u64 * charges.len() as u64 * decoy_mult;
    let progress = ProgressBar::new(estimated_total);
    progress.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({per_sec}, ETA {eta})",
        )
        .unwrap()
        .progress_chars("#>-"),
    );

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

    tracing::info!(
        "Wrote {} entries to {}",
        entry_id,
        config.output.display(),
    );

    Ok(())
}
