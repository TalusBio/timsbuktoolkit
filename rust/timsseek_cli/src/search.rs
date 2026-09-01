//! The search command: validate, stage, and drive the pipeline once per raw
//! input, writing one subdirectory per sample plus the run-level artifacts.

use tims_stage::is_remote_uri;
use timsquery::utils::TupleRange;
use timsquery::{
    IndexedTimstofPeaks,
    IndexingCentroidingConfig,
    load_index,
};
use timsseek::scoring::Scorer;
use timsseek::scoring::timings::TimedStep;
use tracing::{
    error,
    info,
    info_span,
};

use crate::cli::SearchArgs;
use crate::config::{
    self,
    Config,
    OutputConfig,
    load_config,
};
use crate::logging::init_tracing;
use crate::output_sink::{
    OutputSink,
    sample_name_from_uri,
};
use crate::run_inputs::{
    LibrarySource,
    ResolvedInputs,
    resolve_run_inputs,
};
use crate::{
    artifacts,
    build_library,
    errors,
    processing,
};

/// Probe the filesystem for everything the run is about to touch, so a missing
/// input or a colliding artifact fails before the heavy analysis rather than
/// after it. Every value it reads was resolved by [`resolve_run_inputs`].
fn validate_inputs(resolved: &ResolvedInputs) -> std::result::Result<(), errors::CliError> {
    info!("Validating inputs and outputs before processing...");

    let ResolvedInputs {
        raw_inputs,
        library,
        calib_lib_uri,
        output_uri,
        overwrite,
    } = resolved;

    match library {
        // Local only: remote resolution happens at open.
        LibrarySource::File(uri) => {
            if !is_remote_uri(uri) && !std::path::Path::new(uri).exists() {
                return Err(errors::CliError::Io {
                    source: "Speclib file does not exist".to_string(),
                    path: Some(uri.clone()),
                });
            }
            info!("✓ Speclib URI: {}", uri);
        }
        // Probed here rather than left to msspeculator, which opens the FASTA
        // after the model has loaded.
        LibrarySource::Fasta(fasta) => {
            if !fasta.exists() {
                return Err(errors::CliError::Io {
                    source: "Sequence database does not exist".to_string(),
                    path: Some(fasta.to_string_lossy().to_string()),
                });
            }
            info!("✓ FASTA to predict a library from: {}", fasta.display());
        }
    }

    if let Some(uri) = calib_lib_uri {
        if !is_remote_uri(uri) && !std::path::Path::new(uri).exists() {
            return Err(errors::CliError::Io {
                source: "Calibration library file does not exist".to_string(),
                path: Some(uri.clone()),
            });
        }
        info!("✓ Calibration library URI: {}", uri);
    }

    // Local-path existence only; remote URIs are resolved by tims_stage at
    // staging time.
    for raw_uri in raw_inputs {
        if !is_remote_uri(raw_uri) && !std::path::Path::new(raw_uri.as_str()).exists() {
            return Err(errors::CliError::Io {
                source: "Raw file does not exist".to_string(),
                path: Some(raw_uri.clone()),
            });
        }
    }
    info!("✓ All {} raw input(s) validated", raw_inputs.len());

    // Remote outputs skip the writability probe -- the upload is the write test.
    if !is_remote_uri(output_uri) {
        let output_dir_path = std::path::Path::new(output_uri);

        match std::fs::create_dir_all(output_dir_path) {
            Ok(_) => {
                info!("✓ Output directory ready");
            }
            Err(e) => {
                return Err(errors::CliError::Io {
                    source: format!("Failed to create output directory: {}", e),
                    path: Some(output_uri.clone()),
                });
            }
        };

        let test_file = output_dir_path.join(".write_test");
        match std::fs::write(&test_file, "test") {
            Ok(_) => {
                let _ = std::fs::remove_file(&test_file);
                info!("✓ Output directory is writable");
            }
            Err(e) => {
                return Err(errors::CliError::Io {
                    source: format!("Output directory is not writable: {}", e),
                    path: Some(output_uri.clone()),
                });
            }
        }
    }

    // Probe every artifact up-front so a collision fails before the heavy
    // analysis rather than after it.
    if !overwrite {
        let collisions = artifacts::probe_collisions(output_uri, raw_inputs)?;
        if !collisions.is_empty() {
            let list = collisions
                .iter()
                .take(8)
                .map(|s| format!("  - {s}"))
                .collect::<Vec<_>>()
                .join("\n");
            let more = if collisions.len() > 8 {
                format!("\n  ... and {} more", collisions.len() - 8)
            } else {
                String::new()
            };
            return Err(errors::CliError::Config {
                source: format!(
                    "{} output artifact(s) already exist; pass --overwrite/-O to replace or remove them:\n{list}{more}",
                    collisions.len()
                ),
            });
        }
        info!("✓ No output collisions (checked local + remote artifacts)");
    } else {
        info!("✓ --overwrite set; existing output artifacts will be replaced");
    }

    info!("All validations passed! Starting processing...");

    Ok(())
}

fn get_frag_range_from_index(
    index: &IndexedTimstofPeaks,
) -> Result<TupleRange<f64>, errors::CliError> {
    // `fragmented_range` panics only when there are no ms2 window groups --
    // catch it and surface a readable "not a DIA run" error instead.
    let result =
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| index.fragmented_range()));
    match result {
        Ok(r) => Ok(r),
        Err(_) => Err(errors::CliError::DataReading {
            source: "Index has no MS2 window groups -- is this a DIA run?".to_string(),
        }),
    }
}

fn process_single_file(
    raw_uri: &str,
    backend: &tims_stage::PerRunTempdir,
    save_sidecar: bool,
    speclib: &timsseek::data_sources::reference_library::ReferenceLibrary,
    calib_lib: Option<&timsseek::data_sources::reference_library::ReferenceLibrary>,
    config: &Config,
    sink: &OutputSink,
    overwrite: bool,
    max_qvalue: f32,
    no_feature_stats: bool,
) -> std::result::Result<timsseek::scoring::PipelineReport, errors::CliError> {
    let file_name = std::path::Path::new(raw_uri)
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or(raw_uri);
    let _file_span = info_span!("file", name = file_name).entered();
    info!("Processing raw input: {}", raw_uri);

    let step = TimedStep::begin("Loading index");
    let (index, index_source) = load_index(
        raw_uri,
        backend,
        save_sidecar,
        IndexingCentroidingConfig::default(),
    )
    .map_err(|e| errors::CliError::Io {
        source: format!("load_index({raw_uri}): {e}"),
        path: Some(raw_uri.to_string()),
    })?;
    // Cache reuse vs raw build is otherwise invisible to the user.
    let load_index_ms = step.finish_with(format_args!("{index_source}")).as_millis() as u64;
    match index.mobility_kind() {
        timscentroid::MobilityKind::Ook0 => info!("Mobility axis: TIMS 1/K0 (searchable)"),
        timscentroid::MobilityKind::Absent => {
            info!("Mobility axis: none -- mobility filter + score features disabled")
        }
        timscentroid::MobilityKind::Unsupported(label) => info!(
            "Mobility axis: unsupported [{label}] -- mobility filter + score features disabled"
        ),
    }
    alloc_track::snap!("Loading index");

    // On-disk caches use `bucket_size=4096`, too large for the tight mz
    // tolerances of Phase 1 / Phase 3 (measured: ~−24% wall at 256).
    const REBUCKET_LEN: usize = 256;
    let step = TimedStep::begin(format_args!("Rebucket at {}", REBUCKET_LEN));
    let index = index.rebucket(REBUCKET_LEN);
    step.finish();

    let fragmented_range = get_frag_range_from_index(&index)?;

    let pipeline = Scorer {
        index,
        broad_tolerance: config.analysis.tolerance.clone(),
        fragmented_range,
    };

    let file_stem = sample_name_from_uri(raw_uri).ok_or_else(|| errors::CliError::Io {
        source: "Unable to derive sample name from URI".to_string(),
        path: Some(raw_uri.to_string()),
    })?;
    let file_output_dir = sink.sample_dir(&file_stem);

    std::fs::create_dir_all(&file_output_dir).map_err(|e| errors::CliError::Io {
        source: format!("Failed to create output subdirectory: {}", e),
        path: Some(file_output_dir.to_string_lossy().to_string()),
    })?;

    if overwrite {
        sink.clear_existing(&file_stem)?;
    }

    let file_output_config = OutputConfig {
        uri: file_output_dir.to_string_lossy().to_string(),
    };

    let report = processing::run_pipeline(
        speclib,
        calib_lib,
        &pipeline,
        config.analysis.chunk_size,
        &file_output_config,
        max_qvalue,
        load_index_ms,
        &config.calibration,
        no_feature_stats,
        config.analysis.rescore_model,
    )
    .map_err(|e| errors::CliError::DataReading {
        source: format!("{}", e),
    })?;

    info!("Successfully processed {}", raw_uri);
    Ok(report)
}

/// Load a library from a local path or remote URI.
///
/// Remote URIs are streamed (not buffered -- speclibs can be multi-GB) to a
/// tempfile keeping the original basename, because `ReferenceLibrary::from_file`
/// sniffs the format by extension.
fn speclib_from_uri(
    uri: &str,
    decoy_strategy: timsseek::DecoyPolicy,
) -> Result<
    (
        timsseek::data_sources::reference_library::ReferenceLibrary,
        Option<tempfile::TempDir>,
    ),
    errors::CliError,
> {
    if !is_remote_uri(uri) {
        let path = std::path::Path::new(uri);
        let lib = timsseek::data_sources::reference_library::ReferenceLibrary::from_file(
            path,
            decoy_strategy,
        )
        .map_err(|e| errors::CliError::Config {
            source: format!("Failed to load speclib {uri}: {:?}", e),
        })?;
        return Ok((lib, None));
    }

    let trimmed = uri.trim_end_matches('/');
    let fname = trimmed
        .rsplit('/')
        .next()
        .filter(|s| !s.is_empty())
        .unwrap_or("speclib.bin");
    let td = tempfile::Builder::new()
        .prefix("timsseek-speclib-")
        .tempdir()
        .map_err(|e| errors::CliError::Io {
            source: format!("speclib tempdir: {e}"),
            path: None,
        })?;
    let local = td.path().join(fname);
    tims_stage::download_to_file(uri, &local).map_err(|e| errors::CliError::Io {
        source: format!("download speclib {uri}: {e}"),
        path: Some(uri.to_string()),
    })?;
    let lib = timsseek::data_sources::reference_library::ReferenceLibrary::from_file(
        &local,
        decoy_strategy,
    )
    .map_err(|e| errors::CliError::Config {
        source: format!("Failed to load speclib {uri}: {:?}", e),
    })?;
    Ok((lib, Some(td)))
}

pub(crate) fn search(args: &SearchArgs) -> std::result::Result<(), errors::CliError> {
    if args.print_default_config {
        print!("{}", config::DEFAULT_CONFIG_TOML);
        return Ok(());
    }
    if let Some(ref path) = args.write_default_config {
        std::fs::write(path, config::DEFAULT_CONFIG_TOML).map_err(|e| errors::CliError::Io {
            source: e.to_string(),
            path: Some(path.to_string_lossy().to_string()),
        })?;
        eprintln!("Wrote default config to {}", path.display());
        return Ok(());
    }

    let config = load_config(args.config.as_deref())?;
    let (mut config, validated) = resolve_run_inputs(args, config)?;

    // Held in `search()`'s scope so the instrumentation flush guard drops after
    // all work completes.
    let _tracing = init_tracing(args, &validated);

    info!("Parsed configuration: {:#?}", config.clone());
    alloc_track::snap!("start");

    validate_inputs(&validated)?;

    // The stale-tempdir sweep runs inside `PerRunTempdir::new`.
    let staging_cfg = config.staging.clone().unwrap_or_default();
    let save_sidecar_flag = staging_cfg.save_sidecar;
    let backend = tims_stage::PerRunTempdir::new(tims_stage::StagingConfig {
        tempdir_root: staging_cfg.tempdir_root.clone(),
        max_prefix_keys: staging_cfg.max_prefix_keys,
        stale_sweep_age_hours: staging_cfg.stale_sweep_age_hours,
    })
    .map_err(|e| errors::CliError::Io {
        source: format!("staging backend: {e}"),
        path: None,
    })?;

    let sink = OutputSink::new(&validated.output_uri)?;

    let config_output_path = sink.root().join(artifacts::CONFIG_USED);

    // Local-only; a remote destination just overwrites on upload.
    if validated.overwrite && config_output_path.exists() {
        std::fs::remove_file(&config_output_path).map_err(|e| errors::CliError::Io {
            source: format!("Failed to remove existing config file: {}", e),
            path: Some(config_output_path.to_string_lossy().to_string()),
        })?;
    }

    let mut run_report = timsseek::scoring::RunReport {
        overwrite: validated.overwrite,
        ..Default::default()
    };
    let mut failed_files: Vec<(String, errors::CliError)> = Vec::new();
    let mut successful_files: Vec<String> = Vec::new();

    // Obtained once, shared across all files.
    info!(
        "Decoy generation strategy: {}",
        config.analysis.decoy_strategy
    );
    let step = TimedStep::begin(match &validated.library {
        LibrarySource::File(_) => "Loading speclib",
        LibrarySource::Fasta(_) => "Predicting speclib",
    });
    let (speclib, _speclib_td, provenance) = match &validated.library {
        LibrarySource::File(uri) => {
            info!("Building database from speclib URI {uri}");
            let (lib, td) = speclib_from_uri(uri, config.analysis.decoy_strategy)?;
            (lib, td, None)
        }
        LibrarySource::Fasta(fasta) => {
            let prediction = build_library::resolve_search_prediction(
                fasta.clone(),
                config.library.as_ref(),
                config.analysis.decoy_strategy,
            );
            let predicted =
                build_library::predict_in_memory(&prediction, config.analysis.decoy_strategy)?;
            (predicted.library, None, Some(predicted.provenance))
        }
    };
    let load_speclib_ms = step
        .finish_with(format_args!(
            "{} entries, {:.1} frags/entry",
            speclib.len(),
            speclib.avg_fragments()
        ))
        .as_millis() as u64;
    alloc_track::snap!("Loading speclib");

    let (calib_lib, _calib_td, load_calib_lib_ms) = match &validated.calib_lib_uri {
        Some(uri) => {
            let step = TimedStep::begin("Loading calib lib");
            info!("Loading calibration library from {}", uri);
            let (lib, td) = speclib_from_uri(uri, config.analysis.decoy_strategy)?;
            let ms = step
                .finish_with(format_args!(
                    "{} entries, {:.1} frags/entry",
                    lib.len(),
                    lib.avg_fragments()
                ))
                .as_millis() as u64;
            alloc_track::snap!("Loading calib lib");
            (Some(lib), td, ms)
        }
        None => (None, None, 0),
    };

    run_report.load_speclib_ms = load_speclib_ms;
    run_report.load_calib_lib_ms = load_calib_lib_ms;
    run_report.speclib_entries = speclib.len();
    run_report.calib_lib_entries = calib_lib.as_ref().map_or(0, |l| l.len());

    // Written once the library is in hand, not before: a run with no library
    // file is traceable only through the provenance of what it predicted, and
    // there is nothing to record until the prediction has run.
    config.library_provenance = provenance;
    let config_json =
        serde_json::to_string_pretty(&config).map_err(|e| errors::CliError::ParseError {
            msg: format!("Failed to serialize config: {}", e),
        })?;
    std::fs::write(&config_output_path, config_json).map_err(|e| errors::CliError::Io {
        source: e.to_string(),
        path: Some(config_output_path.to_string_lossy().to_string()),
    })?;
    info!("Wrote final configuration to {:?}", config_output_path);

    let total_files = validated.raw_inputs.len();
    info!("Processing {} raw input(s)", total_files);

    for (idx, raw_uri) in validated.raw_inputs.iter().enumerate() {
        info!(
            "Processing input {} of {}: {}",
            idx + 1,
            total_files,
            raw_uri
        );

        let sample_name = match sample_name_from_uri(raw_uri) {
            Some(s) => s,
            None => {
                let e = errors::CliError::Io {
                    source: "Unable to derive sample name from URI".to_string(),
                    path: Some(raw_uri.clone()),
                };
                error!("Failed to process {}: {}", raw_uri, e);
                failed_files.push((raw_uri.clone(), e));
                continue;
            }
        };

        // Header/footer pair; `processing::run_pipeline` phase output lands
        // between them, so batched runs still show per-input wall time.
        println!("=== [{}/{}] {} ===", idx + 1, total_files, sample_name);
        let file_start = std::time::Instant::now();
        let sample_dest = sink.dest_uri_for(&sample_name);

        match process_single_file(
            raw_uri,
            &backend,
            save_sidecar_flag,
            &speclib,
            calib_lib.as_ref(),
            &config,
            &sink,
            validated.overwrite,
            args.max_qvalue,
            args.no_feature_stats,
        ) {
            Ok(report) => {
                if let Err(e) = sink.finalize_sample(&sample_name) {
                    error!("Failed to finalize sample {}: {}", sample_name, e);
                    println!(
                        "=== [{}/{}] {} failed upload after {:?} ===",
                        idx + 1,
                        total_files,
                        sample_name,
                        file_start.elapsed()
                    );
                    run_report.status = timsseek::scoring::timings::RunStatus::Aborted;
                    run_report.abort_reason =
                        Some(format!("upload failure on sample {sample_name}: {e}"));
                    failed_files.push((raw_uri.clone(), e));
                    error!("Aborting batch due to upload failure");
                    break;
                }
                println!("Output: {sample_dest}");
                println!(
                    "=== [{}/{}] {} done in {:?} ===",
                    idx + 1,
                    total_files,
                    sample_name,
                    file_start.elapsed()
                );
                successful_files.push(raw_uri.clone());
                let mut outputs = vec![format!("{sample_dest}/{}", artifacts::RESULTS_PARQUET)];
                if !args.no_feature_stats {
                    outputs.push(format!("{sample_dest}/{}", artifacts::FEATURE_STATS_TSV));
                    outputs.push(format!(
                        "{sample_dest}/{}",
                        artifacts::FEATURE_IMPORTANCE_TSV
                    ));
                }
                run_report.files.push(timsseek::scoring::FileReport {
                    file_name: raw_uri.clone(),
                    pipeline: report,
                    outputs,
                });
            }
            Err(e) => {
                error!("Failed to process {}: {}", raw_uri, e);
                println!(
                    "=== [{}/{}] {} FAILED after {:?}: {} ===",
                    idx + 1,
                    total_files,
                    sample_name,
                    file_start.elapsed(),
                    e
                );
                // I/O errors are likely systemic (disk full, permissions) --
                // abort the batch instead of failing every remaining file.
                if matches!(e, errors::CliError::Io { .. }) {
                    run_report.status = timsseek::scoring::timings::RunStatus::Aborted;
                    run_report.abort_reason = Some(format!("I/O error on {raw_uri}: {e}"));
                    failed_files.push((raw_uri.clone(), e));
                    error!("Aborting batch due to I/O error");
                    break;
                }
                failed_files.push((raw_uri.clone(), e));
            }
        }
    }

    let finalize_step = TimedStep::begin("Finalize run");
    // Must happen before serialization so the report self-describes where to
    // fetch everything.
    run_report.artifacts = artifacts::RUN_ARTIFACTS
        .iter()
        .map(|artifact| sink.dest_uri_for(artifact))
        .collect();
    let run_report_path = sink.root().join(artifacts::RUN_REPORT);
    if let Ok(json) = serde_json::to_string_pretty(&run_report) {
        let _ = std::fs::write(&run_report_path, json);
        info!("Wrote run report to {:?}", run_report_path);
    }

    sink.finalize_run(artifacts::RUN_ARTIFACTS)?;
    finalize_step.finish();

    info!("Successfully processed {} file(s)", successful_files.len());
    if !failed_files.is_empty() {
        error!("Failed to process {} file(s):", failed_files.len());
        for (file, err) in &failed_files {
            error!("  {}: {}", file, err);
        }
        return Err(errors::CliError::Config {
            source: format!("Failed to process {} file(s)", failed_files.len()),
        });
    }

    Ok(())
}
