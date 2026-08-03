# Scoring Pipeline Data Flow

The flow, not the details — counts, defaults and hyperparameters live in the
code and are expected to drift. Paths are relative to `rust/`.

## Shape

A global two-pass search: score everything broadly, calibrate off the best
hits, then re-extract everything with the calibrated windows.

| Phase | Entry point | Does |
|---|---|---|
| 1. Prescore | `phase1_prescore` — `timsseek_cli/src/processing.rs` | Broad, uncalibrated extraction of the whole library; keeps the best-scoring peptides as calibrants. |
| 2. Calibrate | `calibrate_from_phase1` | Fits the iRT→RT curve from the calibrants, measures m/z and mobility error, derives the narrow tolerances. |
| 3. Score | `phase3_score` | Re-extracts every peptide at its calibrated RT with the derived tolerances and computes the full feature set. |
| 4. Compete | `target_decoy_compete` | Dedups by sequence, competes within decoy groups, adds log-space group-difference features. |
| 5. Rescore | `ml::rescore_with` — `timsseek/src/ml/mod.rs` | Fits a discriminant over feature projections and assigns q-values. LDA/MLP stream raw rows; GBM retains its input frame. |
| 6. Write | `execute_pipeline` — `timsseek_cli/src/processing.rs` | Writes filtered results and optional feature-statistics sidecars, then completes the performance report. |

`main()` → `process_single_file()` loads the raw index, builds a `Scorer`, and
calls `processing::run_pipeline()` → `execute_pipeline()`, which runs all six
phases and writes a performance report.

The library is columnar: the scorer iterates `ReferenceLibrary` through a
borrowed flyweight rather than owning per-peptide structs. `--calib-lib` lets a
second library supply the calibrants, with a warning when the two RT scales
barely overlap.

Parallelism is compile-time gated through `utils/maybe_par.rs`, so
`--no-default-features` (dropping `rayon`) gives the serial path.

## Apex scoring

Prescore and full scoring share `compute_traces()`
(`scoring/apex_finding.rs`), which produces the per-cycle elution traces and a
composite apex profile; the apex is that profile's argmax.

* **Phase 1** stops there. Its score is the apex profile peak itself, *not* an
  evidence term — that value only has to rank calibrants.
* **Phase 3** additionally computes `ApexEvidence` (extraction-global, so
  computed once) and evaluates the apex-local features at the chosen cycle.
  `main_score` combines them multiplicatively — evidence scaled by a product
  over the apex features.

Phase 3 then refines in two passes (query at the apex for observed mobility,
re-query tightly around it for main and isotope patterns) and assembles
`ScoringFields`.

## Rescoring

Four interchangeable rescorers in `ml/qvalues.rs`, selected by the
`rescore_model` config field / `--rescore-model` flag (CLI wins). All share the
same input, the same canonical sort + seeded shuffle, and the same FDR tail
(`q = cummin(decoys / targets)`); only the discriminant differs.

| Model | Feature lane | Fit |
|---|---|---|
| `gbm` | all | cross-validated gradient boosting (forust) |
| `lda` | linear | closed-form shrinkage LDA, cross-fit |
| `hybrid` | nonlinear + `lda_score` | per-fold cross-fit LDA feeding the GBM |
| `mlp` (default) | all | cross-validated MLP (`ml/mlp.rs`) |

Features are extracted by a **lane walk**, not a per-result method. GBM retains
the resulting flat row-major frame; LDA and MLP stream the same row projection
through scratch buffers. Projection happens over the already-shuffled slice so
row *i* aligns with candidate *i*, because fold assignment and labels are
positional. Each row concatenates the contributing blocks' fixed-size lane
arrays, so the width is a compile-time constant; names come from the same walk
order.

Walk order per lane: scoring blocks (composition order) → `ResultMeta` →
`Derived` → (nonlinear only) `sequence_counts`. `sequence_counts` is
unconditional — an unparsed peptide contributes NaN rather than a shorter row,
so every row is the same width. Those AA counts have no linear lane, so LDA
never sees them.

The hybrid cross-fits its extra column under the same fold ASSIGNMENT the GBM
uses internally, so a candidate's `lda_score` never comes from a
model that saw it. The two train/score SPLITS deliberately differ (the cross-fit
fits on all-but-one fold, the GBM scorer on one fold at a time) and both are
leak-free; what has to agree is only which rows land in which fold, and both
sides read that from `RowMajorDataset::get_fold`.

## The candidate chain

`ScoredCandidate` (Phase 3) → `CompetedCandidate` (`+ delta_group_ln1p_diff`, `+ delta_group_ln1p_ratio`) →
`FinalResult` (`+ discriminant_score`, `qvalue`) → parquet. All three live in
`scoring/results.rs`; `ScoringFields` is composed there from an ordered block
list, and that order is what parquet columns and the ML lanes both follow.
