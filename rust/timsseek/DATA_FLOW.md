# Scoring Pipeline Data Flow

## Overview

The timsseek scoring pipeline uses a global two-pass architecture:

```
Phase 1: Broad prescore → collect top-2000 calibrants
Phase 2: Calibrate iRT→RT + derive tolerances from calibrant errors
Phase 3: Narrow calibrated extraction → full scoring → results
Post:    Target/decoy competition → ML rescoring → parquet output
```

### Conceptual Decomposition

The pipeline optimizes three independent concerns per peptide:

1. **Apex finding** — locate the correct elution peak in time. Uses the composite apex profile (`cos³ × I × (0.5 + S_norm)`) for peak-picking. Quality measured by whether the detected apex RT matches the true elution time.

2. **Feature extraction** — at the detected apex, compute features that discriminate true peptides from false matches. The 11 apex features, apex evidence, lazyscore baseline stats, m/z/mobility errors, and relative intensities all serve this purpose.

3. **Post-processing** — refine target/decoy assignments using inter-peptide relationships. Currently: target-decoy competition (dedup by sequence, compete within decoy groups). Potential: overlapping fragment ion competition (if two candidates share fragments and RT, the weaker one is likely false).

## Entry Point

**`main()`** — `timsseek_cli/src/main.rs`

1. Parse CLI args (`Cli` struct — `cli.rs`)
   - `--speclib-file` — main spectral library
   - `--calib-lib` — optional calibration library (different RT scale OK)
   - `--dotd-files` — raw instrument files
   - `--output-dir`, `--overwrite`, `--config`, `--decoy-strategy`
2. Load/merge config (JSON config + CLI overrides)
3. Validate inputs (speclib, raw files, output dir, calib lib if provided)
4. For each .d file → `process_single_file()`:
   - Load raw data as `IndexedTimstofPeaks` via `load_index_auto()`
   - Extract DIA fragmentation m/z range from frame reader
   - Build `ScoringPipeline { index, tolerances, fragmented_range }`
   - Call `process_speclib()`

**`process_speclib()`** — `timsseek_cli/src/processing.rs`

1. Load main speclib: `Speclib::from_file(path, decoy_strategy)` → `Vec<QueryItemToScore>`
2. Optionally load calibration library (`--calib-lib`)
3. If calib lib provided: `check_rt_scale_compatibility()` — warns if RT ranges have < 50% overlap
4. Call `main_loop(speclib, calib_lib, pipeline, chunk_size, output)`

## Phase 1: Broad Prescore

**Goal:** Find the top-2000 best-scoring peptides for calibration.

```
main_loop()
  → phase1_prescore()                          processing.rs
    → pipeline.prescore_batch()                pipeline.rs
      → pipeline.prescore()                    pipeline.rs        [per peptide]
        → build_candidate_context()            pipeline.rs
        → buffer.find_apex_location()          apex_finding.rs
```

Both `prescore_batch` and `score_calibrated_batch` have serial paths gated behind `--features serial_scoring` for instrumentation/debugging.

### `build_candidate_context()`
**In:** `QueryItemToScore` (peptide with iRT, m/z, fragments, expected intensities)
**Out:** `(PeptideMetadata, ScoringContext<IonAnnot>)`

1. Compute RT range from `self.tolerances.prescore` (broad: ±5 min or unrestricted)
2. Create `ChromatogramCollector` — allocates cycle×ion intensity arrays
3. `index.add_query(&mut agg, &prescore_tolerance)` — extract peaks from raw data
4. `filter_zero_intensity_ions()` — drop ions with no signal
5. `select_top_n_fragments(n=8)` — keep only top-8 by predicted intensity

### `find_apex_location()`
**In:** `ScoringContext` (chromatogram data + expected intensities)
**Out:** `ApexLocation { score, retention_time_ms, apex_cycle, rising_cycles, falling_cycles }`

1. **`compute_pass_1()`** — single pass over all fragments × cycles:
   - **Cosine**: `dot(obs, sqrt(expected)) / (||obs|| × ||sqrt(expected)||)` per cycle
   - **Scribe**: SSE of sqrt-normalized observed vs predicted distributions → `-ln(SSE)`
   - **Lazyscore**: `lnfact(Σ ln(intensity))` per cycle
   - **Log-intensity**: `ln(1 + Σ raw_intensity)` per cycle
   - **Precursor trace**: summed MS1 precursor intensity (keys ≥ 0)

2. **`compute_main_score_trace()`** — apex profile:
   ```
   C(t) = cos(t)³ × I(t)
   S(t) = scribe(t) × I(t)
   S_norm = (S - min(S)) / (max(S) - min(S))
   apex_profile(t) = C(t) × (0.5 + S_norm(t))
   ```

3. **Peak-pick** on apex profile → argmax (`suggest_apex`). The returned
   `ApexLocation.score` is the **apex profile peak value itself** — no evidence
   or feature term is computed here. That value is what ranks calibrants in
   Phase 1. (`ApexEvidence` is a Phase-3-only concern; see below.)

### `prescore_batch()`
- Parallel (or serial with `--features serial_scoring`) iteration
- Per-thread bounded min-heap (`CalibrantHeap`, capacity=2000)
- Pushes `CalibrantCandidate { score, apex_rt_seconds, speclib_index }`
- Merge heaps across threads → top-2000 globally

## Phase 2: Calibration

**Goal:** Fit iRT→RT curve, measure instrument errors, derive tolerances.

```
main_loop()
  → [build precursor+fragment lookup if calib lib]    processing.rs
  → calibrate_from_phase1()                           processing.rs
```

### Calib lib matching (when `--calib-lib` is used)

Builds a lookup: `HashMap<(quantized_mz_0.01Da, charge), Vec<(rt, sorted_fragment_mzs)>>` from the main speclib. For each calibrant from the calib lib, matches by:
1. Same precursor m/z (within 0.01 Da) + same charge
2. ≥ 5 shared fragment masses (within 0.01 Da, sorted merge)
3. Among matches, pick the one with most shared fragments (break ties by closest RT)

This maps the calibrant's observed apex RT to the **main speclib's iRT** for curve fitting.

### Step A: Fit iRT → RT curve

For each calibrant:
- `x = iRT` (from main speclib if using calib lib, else from phase1 lib)
- `y = observed apex RT` (from Phase 1)

Call `calibrate_with_ranges()` (calibrt crate):
- 100×100 grid, AND-intersection NMS
- Non-suppressed cells use **weighted centroids** (not bin centers)
- Pathfinding: optimal ascending path with bounded lookback L=30
- Result: piecewise-linear `CalibrationCurve`

### Step B: Measure m/z and mobility errors

For each calibrant, re-query the index at apex RT with broad tolerance:
```rust
SpectralCollector<IonAnnot, MzMobilityStatsCollector>
```
Extract observed precursor m/z and mobility:
- `mz_error_ppm = (obs_mz - expected_mz) / expected_mz × 1e6`
- `mobility_error_pct = (obs_mob - expected_mob) / expected_mob × 100`

### Step C: Derive tolerances

- **RT**: `tolerance = rt_sigma_factor × MAD(|residuals|) / 60` (min 0.5 min)
- **m/z**: `asymmetric_tolerance(errors, sigma=2.0)` → `(left_ppm, right_ppm)`
- **mobility**: `asymmetric_tolerance(errors, sigma=3.0)` → `(left_pct, right_pct)`

Asymmetric formula: `left = max(min_val, -(mean - σ×std))`, `right = max(min_val, mean + σ×std)`

**Output:** `CalibrationResult { cal_curve, rt_tolerance_minutes, mz_tolerance_ppm, mobility_tolerance_pct }`

## Phase 3: Calibrated Scoring

**Goal:** Re-extract every peptide with calibrated RT + derived tolerances, compute full features.

```
main_loop()
  → phase3_score()                             processing.rs
    → pipeline.score_calibrated_batch()        pipeline.rs
      → score_calibrated_extraction()          pipeline.rs        [per peptide]
        → build_calibrated_context()           pipeline.rs
        → buffer.find_apex()                   apex_finding.rs
        → execute_secondary_query()            pipeline.rs
        → finalize_results()                   pipeline.rs
```

### `build_calibrated_context()`
**In:** `QueryItemToScore` + `CalibrationResult`
**Out:** `(PeptideMetadata, ScoringContext<IonAnnot>)`

Key difference from Phase 1:
- RT = `calibration.convert_irt(original_irt)` (not the raw iRT)
- Tolerance = `calibration.get_tolerance(mz, mobility, rt)` (narrow, asymmetric)

### `find_apex()` — full scoring
**In:** `ScoringContext` (narrow extraction)
**Out:** `ApexBlocks`

Same as `find_apex_location()` (compute_traces + suggest_apex), then additionally:

1. **`compute_apex_evidence()`** — extraction-global (cycle-invariant), so it is
   computed once and handed to `score_at`: `ApexEvidence` with 9 component fields.
2. **`score_at()`** at the suggested apex cycle:
   - Peak-pick with delta scoring (mask primary, find 2nd/3rd peaks)
   - **11 apex features** at the apex (`compute_apex_features()`):

     | Feature | Description |
     |---------|------------|
     | peak_shape | 0.5×symmetry + 0.5×sharpness around apex |
     | ratio_cv | 1/(1+CV) of obs/predicted ratios at apex |
     | centered_apex | 1 - |apex - center| / (n_cycles/2) |
     | precursor_coelution | Pearson(precursor, summed_fragments) ±10 cycles |
     | fragment_coverage | Fraction of fragments with intensity > 0 |
     | precursor_apex_match | 0.5×proximity + 0.5×fraction |
     | xic_quality | Min XIC signal quality across fragments |
     | fragment_apex_agreement | Pearson(expected, observed per-fragment peaks) |
     | isotope_correlation | Correlation between precursor isotope patterns |
     | gaussian_correlation | Correlation vs Gaussian peak shape |
     | per_frag_gaussian_corr | Mean per-fragment Gaussian correlation |

   - **Weighted product score**:
     ```
     score = apex_evidence × Π(offset_k + scale_k × feature_k)
     ```
     With 11 (offset, scale) pairs from `SCORING_WEIGHTS`

### `execute_secondary_query()` — two-pass spectral refinement
**In:** `QueryItemToScore` + `ApexBlocks`
**Out:** `(SpectralCollector<MzMobilityStatsCollector>, SpectralCollector<f32>)`

1. **Pass 1**: Query at apex RT with secondary tolerance → observed mobility
2. **Pass 2**: Query at apex RT + observed mobility with 3% mobility tolerance
   - Main pattern → `MzMobilityStatsCollector` (m/z, mobility, intensity stats)
   - Isotope pattern (+1 Da offset) → `f32` intensities

### `finalize_results()`
**In:** metadata + `ApexBlocks` + secondary collectors
**Out:** `ScoredCandidate` (wrapping `ScoringFields`)

Via `ScoringFields::compute(FinalizeInputs { .. })`:
1. `MzMobilityOffsets::new()` — top-3 MS1 + top-7 MS2 m/z/mobility errors
2. `RelativeIntensityCollector::new()` — log-normalized MS1/MS2 intensities
3. `compute_secondary_lazyscores()` — main + isotope lazyscores + ratio
4. Assemble the ~90 fields from the `ApexBlocks` blocks (evidence, features,
   primary, counts, intensities, apex_lazy) plus the offsets/intensity blocks

The candidate then walks `ScoredCandidate` → `CompetedCandidate` (after
target/decoy competition) → `FinalResult` (after rescoring, written to parquet).

## Post-Processing

**`target_decoy_compete()`** — `processing.rs`
1. Deduplicate by (sequence, charge, m/z) — keep best score
2. Sort by (decoy_group_id, charge, score desc)
3. Compute delta_group scores between target/decoy pairs
4. Keep one winner per (decoy_group_id, charge)

**Phase 5: rescoring** — `timsseek/src/ml/qvalues.rs`

Three interchangeable rescorers, selected by the `rescore_model` config field /
`--rescore-model` CLI flag (CLI wins). All three share the same input
(`Vec<CompetedCandidate>`), the same canonical sort + seeded shuffle
(`RESCORE_SHUFFLE_SEED = 42`), and the same FDR tail
(`assign_qval`, `q = cummin(decoys / targets)`) — only the discriminant-score
source differs.

| `rescore_model` | Function | Feature lane | Model |
|---|---|---|---|
| `gbm` (default) | `rescore()` | ALL (`build_all_frame`) | 3-fold CV gradient boosting (forust) |
| `lda` | `rescore_lda()` | LINEAR (`build_lane_frame(Linear)`) | single closed-form shrinkage LDA, no CV (~100× cheaper) |
| `hybrid` | `rescore_hybrid()` | NONLINEAR + `lda_score` | per-fold cross-fit LDA on the linear lane, its score pushed as one extra column into the GBM frame |

Feature extraction is the **lane walk**, not a per-result feature method: a
batch `FeatFrame` is built by `build_lane_frame` / `build_all_frame` over the
already-shuffled candidate slice, so frame row *i* aligns with candidate *i*
(the frame MUST be built after the shuffle — fold assignment and labels are
positional). Column names come from the same walk via
`linear_feature_name_set` / `nonlinear_feature_name_set` /
`all_feature_name_set`, which is why columns and names cannot desync (asserted
with `debug_assert_eq!` at each call site).

Walk order per lane: scoring blocks (composition order) → `ResultMeta` →
`Derived` → (nonlinear only) `sequence_counts`, the last of which is gated on
the speclib having parsed sequences.

For `hybrid`, the cross-fit uses `fold(i) = i % 3`, matching the GBM's internal
fold assignment, so a candidate's `lda_score` never comes from an LDA that saw
it.

**Output:** Parquet file via `ResultParquetWriter`

## Feature Vector (combined lane feature set)

The union of the LINEAR and NONLINEAR lanes — i.e. what `all_feature_name_set` /
`build_all_frame` produce, the GBM's input. The LDA sees the linear subset; the
hybrid GBM sees the nonlinear subset plus `lda_score`.

**Context features:**
- `precursor_mz` (binned to 5 Da), `charge`, `mobility_query`, `rt_query` (rounded)
- `nqueries` (fragment count after filtering)

**Primary scores:**
- `main_score`, `main_score/delta_next`, `delta_next`, `delta_second_next`

**RT/mobility deltas:**
- `obs_rt`, `obs_mobility`, `delta_theo_rt`, `sq_delta_theo_rt`
- `delta_ms1_ms2_mobility`, `sq_delta_ms1_ms2_mobility`
- `calibrated_sq_delta_theo_rt`, `recalibrated_query_rt`

**Peak shape:**
- `rising_cycles`, `falling_cycles`

**MS2 scores:**
- `npeaks`, `apex_lazyerscore`, `ln(ms2_summed_intensity)`
- `ms2_lazyerscore`, `ms2_isotope_lazyerscore`, `ms2_isotope_ratio`
- `lazyscore_z`, `lazyscore_vs_baseline`

**Apex evidence & apex features (19 values):**
- `ln(apex_evidence)`, `ln(cosine_au)`, `ln(scribe_au)`
- `cosine_cg`, `scribe_cg` (coelution-gradient, per profile)
- `cosine_weighted_coelution`, `cosine_gradient_consistency`
- `scribe_weighted_coelution`, `scribe_gradient_consistency`
- 11 apex features (peak_shape through per_frag_gaussian_corr)

**Per-ion errors (20 values):**
- 7× `ms2_mz_errors_{i}`, 7× `ms2_mobility_errors_{i}`
- 3× `ms1_mz_errors_{i}`, 3× `ms1_mobility_errors_{i}`

**Intensities:**
- `ln(ms1_summed_precursor_intensity)`
- 3× `ms1_inten_ratio`, 7× `ms2_inten_ratio`

**Target-decoy competition:**
- `delta_group`, `delta_group_ratio`

**Interaction features:**
- `main_score × delta_next`
- `apex_evidence × fragment_coverage`

**Summary error features:**
- Mean |ms2_mz_error|, Mean |ms2_mobility_error|
- Mean |ms1_mz_error|, Mean |ms1_mobility_error|

**Derived features:**
- Max fragment intensity ratio (dominance of strongest fragment)

## Key Types

```
QueryItemToScore
├── digest: Peptide                  raw sequence + optional parsed residues/mods
│   ├── raw: Arc<str>                on-disk modified-sequence string
│   ├── parsed: Option<ParsedSequence>  None iff speclib gate=off (all-or-none)
│   ├── decoy: DecoyMarking          Target / ReversedDecoy
│   └── decoy_group: u32
├── query: TimsElutionGroup<IonAnnot> precursor/fragment m/z, RT, mobility
└── expected_intensity: ExpectedIntensities<IonAnnot>

ChromatogramCollector<IonAnnot, f32>
├── fragments: MzMajorIntensityArray  [n_fragments × n_cycles]
├── precursors: MzMajorIntensityArray [n_precursors × n_cycles]
└── eg: TimsElutionGroup<IonAnnot>

ScoreTraces                           6 per-cycle vectors
├── ms2_cosine_ref_sim, ms2_lazyscore, ms2_scribe
├── ms2_log_intensity, ms1_precursor_trace
└── main_score                        composite apex profile

ApexLocation                          Phase 1 result (lightweight)
├── score: f32                        apex_profile peak value (NOT apex_evidence)
├── retention_time_ms: u32
├── apex_cycle: usize
└── rising_cycles, falling_cycles: u8

ApexBlocks                            Phase 3 result (full)
├── joint_apex_cycle: usize
├── retention_time_ms: u32            also the secondary-query bridge input
├── evidence: ApexEvidence            9 component scores
├── features: ApexFeatures            11 apex-local features
├── primary: PrimaryScores            main_score, delta_next, delta_second_next
├── counts: ApexCounts                rising/falling cycles, npeaks
├── intensities: Intensities          ms1/ms2 summed intensity at apex
└── apex_lazy: ApexLazyScores         lazyscore_z, lazyscore_vs_baseline

CalibrationResult                     immutable, speclib not mutated
├── cal_curve: RTCalibration          iRT → calibrated RT
├── rt_tolerance_minutes: f32
├── mz_tolerance_ppm: (f64, f64)     asymmetric
└── mobility_tolerance_pct: (f32, f32) asymmetric

ScoringFields                         ~90 fields, composed from ordered blocks
  ScoredCandidate                     Phase 3 output
  → CompetedCandidate                 + delta_group / delta_group_ratio
  → FinalResult                       + discriminant_score / qvalue → parquet
```

## Constants

| Constant | Value | Location |
|----------|-------|----------|
| TOP_N_FRAGMENTS | 8 | pipeline.rs |
| NUM_MS1_IONS | 3 | scoring/mod.rs |
| NUM_MS2_IONS | 7 | scoring/mod.rs |
| SCRIBE_FLOOR | -100.0 | scores/scribe.rs |
| n_calibrants | 2000 | CalibrationConfig |
| grid_size | 100 | CalibrationConfig |
| rt_sigma_factor | 3.0 | CalibrationConfig |
| min_rt_tolerance_minutes | 0.5 | CalibrationConfig |
| mz_sigma | 3.0 | CalibrationConfig |
| mobility_sigma | 3.0 | CalibrationConfig |
| MIN_SHARED_FRAGMENTS | 5 | processing.rs (calib lib matching) |
| Area-uniqueness hw | 5 | blocks/apex_evidence.rs |
| Coelution hw | 20 | blocks/apex_evidence.rs |
| Gradient hw | 10 | blocks/apex_evidence.rs |
| GBM iterations | 1000 | cv.rs |
| GBM learning_rate | 0.1 | cv.rs |
| GBM max_depth | 6 | cv.rs |
| GBM min_leaf_weight | 5.0 | cv.rs |
| GBM subsample | 0.8 | cv.rs |
| GBM colsample_bytree | 0.8 | cv.rs |
| GBM early_stopping | None | cv.rs |
| CV folds | 3 | qvalues.rs |
| LDA shrinkage (`DEFAULT_LDA_SHRINKAGE`) | 1e-2 | ml/lda.rs |
| LDA row-chunk (`CHUNK_ROWS`, determinism) | 65536 | ml/lda.rs |
| Rescore shuffle seed | 42 | qvalues.rs |
| Pathfinding lookback | 30 | pathfinding.rs |
