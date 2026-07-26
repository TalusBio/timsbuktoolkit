# Scoring Pipeline Data Flow

Orientation only — the code is the spec. Paths are relative to `rust/`.

## Shape

A global two-pass search: score everything broadly, calibrate off the best
hits, then re-extract everything with the calibrated windows.

| Phase | Entry point | Does |
|---|---|---|
| 1. Prescore | `phase1_prescore` — `timsseek_cli/src/processing.rs` | Broad, uncalibrated extraction of the whole library; keeps the top-N best-scoring peptides (default 2000) as calibrants. |
| 2. Calibrate | `calibrate_from_phase1` — `processing.rs` | Fits the iRT→RT curve from the calibrants, measures m/z and mobility error, derives the narrow tolerances. |
| 3. Score | `phase3_score` — `processing.rs` | Re-extracts every peptide at its calibrated RT with the derived tolerances and computes the full feature set. |
| 4. Compete | `target_decoy_compete` — `processing.rs` | Dedups by sequence, competes within decoy groups, adds `delta_group`. |
| 5. Rescore | `ml::rescore_with` — `timsseek/src/ml/mod.rs` | Fits a discriminant over the feature matrix and assigns q-values. |

## Entry point

`main()` (`timsseek_cli/src/main.rs`) → `process_single_file()`, which loads the
raw index, builds a `Scorer { index, broad_tolerance, fragmented_range }`, and
calls `processing::run_pipeline()`. That is a thin wrapper over
`execute_pipeline()`, which runs all five phases and writes
`performance_report.json`.

The library is columnar: `Speclib` is an alias for `ReferenceLibrary`
(`timsseek/src/data_sources/reference_library.rs`), and the scorer iterates it
through the `RefQuery<'_>` flyweight returned by `item_at()` rather than owning
per-peptide structs. With `--calib-lib`, a second library supplies the
calibrants; `check_rt_scale_compatibility()` warns when the two RT scales
barely overlap.

Parallelism is compile-time gated: everything goes through
`timsseek/src/utils/maybe_par.rs`, so `--no-default-features` (dropping
`rayon`) gives the serial path. There is no `serial_scoring` feature.

## Apex scoring

Both prescore and full scoring live in `timsseek/src/scoring/apex_finding.rs`
and share `compute_traces()`, which produces `ElutionTraces` — six per-cycle
vectors (cosine, lazyscore, scribe, log-intensity, MS1 precursor trace, and the
composite `apex_profile`). `suggest_apex()` argmaxes the apex profile.

* Phase 1 (`find_apex_location`) stops there and returns `ApexLocation`, whose
  `score` is the apex profile peak value itself — not an evidence term. That
  value is what ranks calibrants.
* Phase 3 (`find_apex`) additionally computes `ApexEvidence` (extraction-global,
  so computed once) and calls `score_at()` at the chosen cycle for the 11
  `ApexFeatures` and the delta-scored runner-up peaks, returning `ApexBlocks`.
  `main_score` is `compute_weighted_score` in `scoring/blocks/apex_features.rs`:
  `apex_evidence × Π(offset_k + scale_k × feature_k)` over those 11 features.

Phase 3 then runs `execute_secondary_query()` (a two-pass refinement: query at
the apex to get observed mobility, then re-query tightly around it for main and
isotope patterns) and `finalize_results()`, which assembles `ScoringFields`.

## Rescoring

Three interchangeable rescorers in `timsseek/src/ml/qvalues.rs`, selected by the
`rescore_model` config field / `--rescore-model` flag (CLI wins). All share the
same input, the same canonical sort + seeded shuffle, and the same FDR tail
(`q = cummin(decoys / targets)`); only the discriminant differs.

| Model | Function | Feature lane | Fit |
|---|---|---|---|
| `gbm` (default) | `rescore` | ALL (131 cols) | 3-fold CV gradient boosting (forust) |
| `lda` | `rescore_lda` | LINEAR (102 cols) | one closed-form shrinkage LDA, no CV |
| `hybrid` | `rescore_hybrid` | NONLINEAR (29 cols) + `lda_score` | per-fold cross-fit LDA feeding the GBM |

Feature extraction is a **lane walk**, not a per-result method: a flat row-major
`Vec<f64>` (`feat[i*ncols + j]`) built by `build_{linear,nonlinear,all}_matrix`
over the already-shuffled slice, so row *i* aligns with candidate *i*. The
matrix MUST be built after the shuffle — fold assignment and labels are
positional. Each row concatenates the contributing blocks' lane arrays, so
`ncols` is a compile-time constant; column names come from the same walk order
via `*_feature_name_set`, cross-checked by `lane_matrix_widths_match_name_sets`.

Walk order per lane: scoring blocks (composition order) → `ResultMeta` →
`Derived` → (nonlinear only) `sequence_counts`. `sequence_counts` is
unconditional: a peptide with no parsed sequence contributes `f64::NAN` for all
22 of its features (forust reads NaN as missing) rather than a shorter row.
Those AA counts have no linear lane, so LDA never sees them.

For `hybrid`, the cross-fit uses `fold(i) = i % 3`, matching the GBM's internal
fold assignment, so a candidate's `lda_score` never comes from an LDA that saw
it.

## Key types

| Type | Where |
|---|---|
| `Speclib` = `ReferenceLibrary`, `RefQuery<'_>` | `data_sources/reference_library.rs` |
| `Scorer` | `scoring/pipeline.rs` |
| `ElutionTraces` | `scoring/apex_finding.rs` |
| `ApexLocation` (Phase 1) / `ApexBlocks` (Phase 3) | `scoring/apex_finding.rs` |
| `ApexEvidence` (9) / `ApexFeatures` (11) | `scoring/blocks/` |
| `CalibrationResult` | `rt_calibration.rs` |
| `ScoringFields` — composed from 13 ordered blocks | `scoring/results.rs` |
| `ScoredCandidate` → `CompetedCandidate` → `FinalResult` | `scoring/results.rs` |

The candidate walks that last chain: Phase 3 output, then `+ delta_group` after
competition, then `+ discriminant_score / qvalue` after rescoring, then parquet
via `ResultParquetWriter`.
