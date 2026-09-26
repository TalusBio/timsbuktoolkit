# Entrapment exploration notebook

Base: PR #152, commit `dff06cc`. This branch records exploratory changes;
follow-up branches can re-implement retained ideas cleanly.

## Fixed experiment

- Raw file: `Ast_20240130_Bo_AI_30_2mz_HeLa01.mzML`
- Library: `shitshit/hela_entrapment/entrapment.mzspeclib.txt.gz`
- Pair file: `shitshit/hela_entrapment/peptide_pairs.tsv`
- Search: `--max-qvalue 1 --decoy-strategy never`, MLP rescorer.
- Primary metric: maximum unique original-target peptide/charge IDs at a
  q-value cutoff whose one-fold paired entrapment FDP is at most 1%.
- Record the reported 1% q-value cutoff and its paired FDP as diagnostics.
- Exploration uses this one raw file without validation splits. Results select
  hypotheses; they do not establish general calibration.

## Trial log

| Trial | Change | Empirical 1% FDP target IDs | Reported 1% q target IDs | Paired FDP at reported 1% q | Decision |
|---|---|---:|---:|---:|---|
| Baseline | Current MLP, target/decoy competition | 23,048 | 35,663 | 4.21% | Reference |
| E1 | Omit both decoy-group margin features; keep competition | 31,012 | 34,070 | 1.29% | Keep (+7,964) |
| E2 | E1 plus no target/decoy competition | 30,430 | 31,368 | 1.10% | Roll back (−582 vs E1) |
| E3 | E1 with LDA rescoring | 20,086 | 21,375 | 1.18% | Reject (−10,926 vs E1) |
| E4 | E1 with GBM rescoring | 30,035 | 31,384 | 1.15% | Reject (−977 vs E1) |
| E5 | E1 with hybrid LDA→GBM rescoring | 21,187 | 22,946 | 1.29% | Reject (−9,825 vs E1) |
| E6 | E1 MLP score − 0.10 × ln(1 + main score) | 31,810 | 33,059 | 1.17% | Keep (+798 vs E1) |
| E7 | E6 score + 0.40 × held-out GBM score | 32,864 | 34,816 | 1.28% | Keep (+1,054 vs E6) |
| E8 | E7 score − 1.10 × held-out LDA score | 33,214 | 34,854 | 1.26% | Keep (+350 vs E7) |
| E9 | Cached direct-field additions to E8 | 33,257* | — | — | Reject (only +43 offline) |
| E10 | Widen MLP hidden layers from 32/16 to 64/32 | 31,330 | 34,030 | 1.25% | Roll back (−1,884 vs E8) |

E1 reused the earlier full search in `shitshit/hela_entrapment/search_no_group`.
It scored the same 1,753,338 candidates with identical raw scores as the
baseline; only the MLP inputs changed. The top cutoff satisfying paired FDP
≤1% is reported q=0.007557, with paired FDP 0.996% and 31,012 original targets.

E2 retained both target and decoy candidates after sequence deduplication.
All E1 winners were present with identical raw scores. The paired FDP curve
moved closer to reported q, but E2 lost 582 IDs at the empirical 1% FDP cutoff.
The competition bypass was temporary and is not in this branch.

E3 used `--rescore-model lda` with the E1 binary. It lost substantial target
yield, despite a slightly smaller paired FDP at the reported 1% q cutoff.
E4 used `--rescore-model gbm`; it approached E1 but still lost 977 IDs at the
empirical cutoff.
E5 used `--rescore-model hybrid`; it also lost substantial target yield.

E6 rescored the saved E1 candidates offline across nearby weights, then ran
the best weight through the full engine. The offline code reproduced all
1,753,338 stored q-values at weight zero and predicted exactly 31,810 IDs at
weight −0.10. The full run matched that count. Neighboring weights −0.095 and
−0.105 yielded 31,769 and 31,795 IDs offline. This is a one-file exploratory
weight and should be re-estimated on independent data before production use.

E7 joined saved E6 and E4 candidate scores by unique library ID. The cached
scan reached 32,864 IDs at GBM weight 0.40, with nearby weights 0.35 and 0.45
yielding 32,645 and 32,817. The dual-model full run exactly reproduced 32,864.
It uses the same shuffled candidates and fold assignment for both models.
The reported q≤1% cutoff still exceeds 1% paired FDP (1.28%), so the empirical
cutoff is required for this exploratory result. Adding a further ±0.025 to
±0.15 main-score weight did not improve E7. A saved hybrid-score blend also
failed to improve it.

E8 tried the standalone LDA score as an additive signal. A negative weight
near −1.1 improved the cached E7 result; neighboring weights −0.9 and −1.25
yielded 33,092 and 32,961 IDs. The integrated three-model run reproduced the
cached 33,214 exactly. The negative sign is empirical on this one file and
should not be interpreted as a general property of LDA.

E9 screened standardized direct fields on saved E8 candidates: fragment
coverage, absolute RT error, cosine, scribe, isotope correlation, ratio CV,
fragment apex agreement, MS2 lazy scores, scored fragment count, and logged MS2
intensity. Most best weights were zero. Adding 0.15× standardized scored
fragment count reached 33,257 IDs (*offline only*), a 43-ID gain that did not
justify another full search or a dataset-wide normalization step. No engine
change was retained.

E10 widened only the MLP while keeping E8's fixed score weights. The full run
lost 1,884 IDs at empirical 1% paired FDP, so the width change was reverted.

Run an offline score scan with:

```sh
uv run --group interactive python -m experiments.entrapment_20260925.offline_rank \
  --results shitshit/hela_entrapment/search_no_group/Ast_20240130_Bo_AI_30_2mz_HeLa01.mzML/results.parquet \
  --pairs shitshit/hela_entrapment/peptide_pairs.tsv \
  --alpha 0 -0.095 -0.1 -0.105
```

For each new trial: make one change, run the fixed experiment, append its
metrics and observation here. Commit improvements to the primary metric;
restore unsuccessful code edits while preserving this log.

Run a model/config trial with:

```sh
uv run --group interactive python -m experiments.entrapment_20260925.run_trial \
  --name e3_lda \
  --rescore-model lda \
  --raw /Users/sebastianpaez/data/astral_data/Ast_20240130_Bo_AI_30_2mz_HeLa01.mzML
```

This writes search and FDP files under `shitshit/hela_entrapment/trials/<name>`
and a small metric JSON under this directory's `results/`.
