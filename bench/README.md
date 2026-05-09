# bench

Fixture-driven bench harness for `timsseek`. Each fixture is a TOML in `bench/fixtures/` pointing at S3 URIs.

## Run a fixture

    uv run --group bench python -m bench.wandb_bench hela_iccoff_gt20peps
    uv run --group bench python -m bench.wandb_bench --all
    uv run --group bench python -m bench.wandb_bench --match 'hela*'

Outputs land under `bench_out/` (gitignored): `logs/<name>-<ts>/`, `parquets/<name>-<ts>-classified.parquet`, `plots/<name>-fdr_curve-<ts>.png`. Wandb runs go to `jspaezp/timsseek`.

Fixtures with `entrapment_fasta` set automatically run the entrapment classification + FDR-curve step.

## Push a new fixture

Requires `aws` CLI (auth via env / profile).

    uv run --group bench python -m bench.push_fixture \
      --name hela_iccoff_gt20peps \
      --bucket terraform-workstations-bucket --prefix jspaezp/timsseek_fixtures \
      --db ~/fasta/hela_gt20peps.fasta \
      --raw ~/data/decompressed_timstof/250225_Desnaux_200ng_Hela_ICC_off_DIA.d \
      --config bench/configs/default.toml \
      --koina-url http://localhost:8501/v2/models   # omit for public Koina

`--db` (and `--entrap-db`, `--calib-db`) are repeatable and accept any of: local `*.fasta(.gz)` path, local `*.txt` accession list, `s3://...` URI, `UPxxxxxxxxx` proteome ID, bare uniprot accession. After upload, hand-edit the generated `bench/fixtures/<name>.toml` to add a description, then `git add bench/fixtures/<name>.toml`.

Re-running `push_fixture` is idempotent by default: existing S3 objects are skipped (single files via `aws s3 ls` check; `.d` directory via `aws s3 sync`). Pass `--force` to re-upload everything.

## Stage a fixture for offline / repeated runs

When iterating on a fixture, pull its inputs to a local cache once, then run against the staged copy:

    uv run --group bench python -m bench.stage_fixture hela_iccoff_gt20peps
    uv run --group bench python -m bench.wandb_bench --fixtures-dir bench_out/staged hela_iccoff_gt20peps

`stage_fixture` defaults: cache root `bench_out/cache/<name>/` (override via `--cache-dir` or `BENCH_CACHE_DIR` env), output TOML `bench_out/staged/<name>.toml` (override via `--out`). Already-cached files are skipped on re-stage; pass `--force` to re-download. Inputs that are already absolute local paths are referenced as-is (no copy).

## Schema

See `bench/_fixture_schema.py` for the canonical TOML schema. Inputs accept `s3://` URIs or absolute local paths (the latter for staged fixtures only — `push_fixture` always emits `s3://`).
