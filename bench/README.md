# bench

Fixture-driven bench harness for `timsseek`. Each fixture is a TOML in `bench/fixtures/` pointing at S3 URIs.

## Run a fixture

    uv run --group bench python -m bench.wandb_bench hela_iccoff_gt20peps
    uv run --group bench python -m bench.wandb_bench --all
    uv run --group bench python -m bench.wandb_bench --match 'hela*'

Outputs land under `bench_out/` (gitignored): `logs/<name>-<ts>/`, `parquets/<name>-<ts>-classified.parquet`, `plots/<name>-fdr_curve-<ts>.png`, `plots/<name>-mainscore_hist-<ts>.png`. Wandb runs go to `jspaezp/timsseek`.

Fixtures with `entrapment_peptides` set automatically run entrapment classification, emit lower-bound + combined FDR estimators (and matched FDR if `pairing` is set), and a 4-group `main_score` histogram (target/decoy × class=target/entrap).

## Push a new fixture

Requires `aws` CLI (auth via env / profile).

### Foreign-species entrapment (Algorithm 2 of Noble et al, FDRBench paper)

    uv run --group bench python -m bench.push_fixture \
      --name hela_iccoff_human_yeast \
      --bucket terraform-workstations-bucket --prefix jspaezp/timsseek_fixtures \
      --db UP000005640 \
      --entrap-db UP000002311 \
      --raw ~/data/decompressed_timstof/250225_Desnaux_200ng_Hela_ICC_off_DIA.d \
      --config bench/configs/default.toml \
      --entrap-ratio 1.0 \
      --request-delay-ms 250

Pipeline: digest target + entrap (trypsin, 1 missed cleavage), filter to length 7-30, drop entrap peptides that also appear in target, randomly subsample entrap to `r × |target|` (seed=42; override with `--seed`). Uploads `target.peptides.txt`, `entrap.peptides.txt`, `database.peptides.txt` (union) and builds the speclib via `speclib_build_cli --peptides s3://.../database.peptides.txt`.

Records `entrapment_mode = "foreign"` and the actual achieved `entrapment_ratio` on the fixture.

### Shuffled entrapment (Algorithm 1 — paired estimator)

    uv run --group bench python -m bench.push_fixture \
      --name hela_iccoff_shuffled \
      --bucket terraform-workstations-bucket --prefix jspaezp/timsseek_fixtures \
      --db UP000005640 \
      --entrap-db SHUFFLED \
      --raw ~/data/decompressed_timstof/250225_Desnaux_200ng_Hela_ICC_off_DIA.d \
      --config bench/configs/default.toml \
      --entrap-ratio 1.0 \
      --request-delay-ms 250

Pipeline: digest target, length-filter, then for each surviving target peptide generate r distinct shuffles (interior permuted, C-term residue fixed). Targets that can't produce r unique shuffles are dropped. With `--entrap-ratio 1.0` the runner also emits `pairing.tsv` enabling the matched FDP estimator.

Records `entrapment_mode = "shuffled"` plus `pairing` URI when r=1.

### Common flags

`--db` (and `--entrap-db`, `--calib-db`) accept: local `*.fasta(.gz)`, local `*.txt` accession list, `s3://...` URI, `UPxxxxxxxxx` proteome ID, bare UniProt accession, or the literal `SHUFFLED` for `--entrap-db` only.

`fetch_proteome` defaults to `reviewed:true` (Swiss-Prot only). Pass full proteomes via fasta or accession list if you need TrEMBL.

Other flags: `--peptide-min-len 7`, `--peptide-max-len 30`, `--missed-cleavages 1`, `--seed 42`. `--request-delay-ms` throttles speclib_build_cli's koina calls (default 500).

After upload, hand-edit the generated `bench/fixtures/<name>.toml` to add a description, then `git add bench/fixtures/<name>.toml`.

Re-running `push_fixture` is idempotent by default: existing S3 objects are skipped (single files via `aws s3 ls`; `.d` directory via `aws s3 sync`). Pass `--force` to re-upload everything.

## FDP estimators

Per Noble et al, FDRBench paper (Table S2). When entrapment is configured:

| Estimator | Formula | Available when |
|---|---|---|
| Lower bound | `n_e / (n_e + n_t)` | always |
| Combined (avg upper bound) | `n_e × (1 + 1/r) / (n_e + n_t)` | always |
| Matched (k=1, avg upper bound) | `(n_e + n_p_s_t + 2·n_p_t_s) / (n_e + n_t)` | shuffled mode + r=1 |

`r = entrapment_ratio` is recorded on the fixture from the actual ratio achieved at push time. Counts are walked over `is_target=True` rows only (post-competition target winners; decoy wins are TDA-style FPs, separate from entrapment FPs).

`compute_fdr_curve` emits all available estimator columns; `plot_fdr_curve` overlays them.

## Stage a fixture for offline / repeated runs

When iterating on a fixture, pull its inputs to a local cache once, then run against the staged copy:

    uv run --group bench python -m bench.stage_fixture hela_iccoff_human_yeast
    uv run --group bench python -m bench.wandb_bench --fixtures-dir bench_out/staged hela_iccoff_human_yeast

Defaults: cache root `bench_out/cache/<name>/` (override via `--cache-dir` or `BENCH_CACHE_DIR` env), output TOML `bench_out/staged/<name>.toml` (override via `--out`). Already-cached files are skipped on re-stage; `--force` re-downloads. Inputs that are already absolute local paths are referenced as-is (no copy).

## Schema

See `bench/_fixture_schema.py` for the canonical TOML. Required: `inputs.target_peptides`, `inputs.speclib`, `inputs.raw`. Optional: `inputs.entrapment_peptides` (with `entrapment_ratio` and `entrapment_mode`), `inputs.pairing` (only valid when `entrapment_mode = "shuffled"`), `inputs.calibration_speclib`. URIs accept `s3://` or absolute local paths (the latter for staged fixtures; `push_fixture` always emits `s3://`).
