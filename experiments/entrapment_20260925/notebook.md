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

E1 reused the earlier full search in `shitshit/hela_entrapment/search_no_group`.
It scored the same 1,753,338 candidates with identical raw scores as the
baseline; only the MLP inputs changed. The top cutoff satisfying paired FDP
≤1% is reported q=0.007557, with paired FDP 0.996% and 31,012 original targets.

For each new trial: make one change, run the fixed experiment, append its
metrics and observation here. Commit improvements to the primary metric;
restore unsuccessful code edits while preserving this log.
