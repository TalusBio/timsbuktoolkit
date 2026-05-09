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
      --bucket timsbukto-bench --prefix fixtures \
      --db ~/fasta/hela_gt20peps.fasta \
      --raw ~/data/decompressed_timstof/250225_Desnaux_200ng_Hela_ICC_off_DIA.d \
      --config bench/configs/default.toml \
      --koina-url http://localhost:8501/v2/models   # omit for public Koina

`--db` (and `--entrap-db`, `--calib-db`) are repeatable and accept any of: local `*.fasta(.gz)` path, local `*.txt` accession list, `s3://...` URI, `UPxxxxxxxxx` proteome ID, bare uniprot accession. After upload, hand-edit the generated `bench/fixtures/<name>.toml` to add a description, then `git add bench/fixtures/<name>.toml`.

## Schema

See `bench/_fixture_schema.py` for the canonical TOML schema. Inputs must be `s3://` URIs (project rule: bench data lives in S3 only).
