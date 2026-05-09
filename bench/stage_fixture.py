"""Stage a fixture for offline / repeated runs.

Reads `bench/fixtures/<name>.toml`, downloads its S3 inputs to a local cache,
and writes a new fixture TOML pointing at those local paths. The new TOML
can be run via `python -m bench.wandb_bench --fixtures-dir <dir> <name>`.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from loguru import logger

from bench._fixture_schema import Fixture, load_fixture
from bench._s3 import s3_download_file, s3_sync_dir
from bench.push_fixture import build_fixture_toml


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("name", help="Fixture name (matches a TOML in --fixtures-dir)")
    p.add_argument(
        "--fixtures-dir",
        type=Path,
        default=Path("bench/fixtures"),
        help="Directory containing the source fixture (default: bench/fixtures/)",
    )
    p.add_argument(
        "--cache-dir",
        type=Path,
        default=Path(os.environ.get("BENCH_CACHE_DIR", "bench_out/cache")),
        help="Local cache root (env: BENCH_CACHE_DIR; default: bench_out/cache/)",
    )
    p.add_argument(
        "--out",
        type=Path,
        help="Output TOML path (default: bench_out/staged/<name>.toml)",
    )
    p.add_argument(
        "--overwrite", action="store_true", help="Replace existing output TOML"
    )
    p.add_argument(
        "--force", action="store_true", help="Re-download files even if cached"
    )
    args = p.parse_args(argv)
    if args.out is None:
        args.out = Path("bench_out/staged") / f"{args.name}.toml"
    return args


def _stage_one_file(uri: str, dst: Path, force: bool) -> str:
    """Resolve `uri` to a local path; return path string for the staged TOML.

    - If `uri` is already an absolute local path, return it unchanged (no copy).
    - If `uri` is `s3://...`, download to `dst` (skip if exists, unless `force`).
    """
    if not uri.startswith("s3://"):
        return uri  # already local; reference as-is
    if dst.exists() and not force:
        logger.info("stage: cached {} (skip)", dst)
        return str(dst)
    s3_download_file(uri, str(dst))
    return str(dst)


def _stage_one_dir(uri: str, dst: Path, force: bool) -> str:  # noqa: ARG001
    """Sync `uri` (s3 prefix) into `dst`. Returns the path string."""
    if not uri.startswith("s3://"):
        return uri
    # `aws s3 sync` is itself idempotent; --force just forces a re-sync
    # which has the same observable result, so we always call it.
    s3_sync_dir(uri, str(dst))
    return str(dst)


def stage(
    *,
    name: str,
    fixtures_dir: Path,
    cache_dir: Path,
    out: Path,
    overwrite: bool,
    force: bool,
) -> None:
    """Stage one fixture for offline use."""
    src_toml = fixtures_dir / f"{name}.toml"
    fx: Fixture = load_fixture(src_toml)

    if out.exists() and not overwrite:
        raise FileExistsError(f"staged TOML already exists: {out} (pass --overwrite)")

    cache_root = cache_dir / name
    cache_root.mkdir(parents=True, exist_ok=True)

    fasta_local = _stage_one_file(
        fx.inputs.fasta, cache_root / "proteome.fasta", force
    )
    speclib_local = _stage_one_file(
        fx.inputs.speclib, cache_root / "lib.msgpack.zst", force
    )
    raw_local = _stage_one_dir(fx.inputs.raw, cache_root / "sample.d", force)

    entrap_local: str | None = None
    if fx.inputs.entrapment_fasta is not None:
        entrap_local = _stage_one_file(
            fx.inputs.entrapment_fasta, cache_root / "entrap.fasta", force
        )

    calib_local: str | None = None
    if fx.inputs.calibration_speclib is not None:
        calib_local = _stage_one_file(
            fx.inputs.calibration_speclib, cache_root / "calib.msgpack.zst", force
        )

    body = _build_staged_toml(
        name=name,
        description=fx.description,
        config=fx.config,
        fasta_uri=fasta_local,
        speclib_uri=speclib_local,
        raw_uri=raw_local,
        entrapment_fasta_uri=entrap_local,
        calibration_speclib_uri=calib_local,
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(body)
    logger.info("staged fixture written: {}", out)
    logger.info(
        "run with: python -m bench.wandb_bench --fixtures-dir {} {}",
        out.parent, name,
    )


def _build_staged_toml(
    *,
    name: str,
    description: str,
    config: dict,
    fasta_uri: str,
    speclib_uri: str,
    raw_uri: str,
    entrapment_fasta_uri: str | None,
    calibration_speclib_uri: str | None,
) -> str:
    """Emit a staged-fixture TOML body. Mirrors push_fixture.build_fixture_toml's
    layout but takes the [config] table as a dict (already loaded) instead of a
    file path."""
    lines: list[str] = []
    lines.append(f'name = "{name}"')
    desc = description.replace('"', '\\"')
    lines.append(f'description = "{desc}"')
    lines.append("")
    lines.append("[inputs]")
    lines.append(f'fasta = "{fasta_uri}"')
    lines.append(f'speclib = "{speclib_uri}"')
    lines.append(f'raw = "{raw_uri}"')
    if entrapment_fasta_uri is not None:
        lines.append(f'entrapment_fasta = "{entrapment_fasta_uri}"')
    if calibration_speclib_uri is not None:
        lines.append(f'calibration_speclib = "{calibration_speclib_uri}"')
    lines.append("")
    lines.append("# === embedded timsseek config ===")
    lines.extend(_emit_config(config))
    lines.append("")
    return "\n".join(lines)


def _emit_config(config: dict, prefix: str = "config") -> list[str]:
    """Render a nested dict as TOML lines under `[<prefix>]` (and sub-tables)."""
    lines: list[str] = []
    scalars: dict = {}
    sub_tables: dict = {}
    for k, v in config.items():
        if isinstance(v, dict):
            sub_tables[k] = v
        else:
            scalars[k] = v
    if scalars:
        lines.append(f"[{prefix}]")
        for k, v in scalars.items():
            lines.append(f"{k} = {_toml_value(v)}")
        lines.append("")
    for k, v in sub_tables.items():
        lines.extend(_emit_config(v, f"{prefix}.{k}"))
    return lines


def _toml_value(v: object) -> str:
    """Minimal TOML serializer for scalars + simple lists + inline tables."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return str(v)
    if isinstance(v, str):
        escaped = v.replace('"', '\\"')
        return f'"{escaped}"'
    if isinstance(v, list):
        return "[" + ", ".join(_toml_value(x) for x in v) + "]"
    if isinstance(v, dict):
        items = ", ".join(f"{k} = {_toml_value(val)}" for k, val in v.items())
        return "{" + items + "}"
    raise ValueError(f"unsupported TOML value: {v!r}")


# Re-export so the name is reachable from tests that import from this module;
# suppresses any "imported but unused" lint on the import above.
__all__ = [
    "parse_args",
    "stage",
    "build_fixture_toml",
]


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    stage(
        name=args.name,
        fixtures_dir=args.fixtures_dir,
        cache_dir=args.cache_dir,
        out=args.out,
        overwrite=args.overwrite,
        force=args.force,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
