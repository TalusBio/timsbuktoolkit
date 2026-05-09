"""Bench runner: load named fixtures, run timsseek, log to wandb."""

from __future__ import annotations

import argparse
import fnmatch
import json
import os
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path

import wandb
from loguru import logger

from bench._fixture_schema import Fixture, load_fixture
from bench._s3 import s3_download_file
from bench.entrapment import analyse

DEFAULT_FIXTURES_DIR = Path("bench/fixtures")
DEFAULT_OUT_ROOT = Path("bench_out")
ENTITY = "jspaezp"
PROJECT = "timsseek"


def _list_fixtures(fixtures_dir: Path) -> list[str]:
    return sorted(p.stem for p in fixtures_dir.glob("*.toml"))


def select_fixtures(
    names: list[str],
    all_: bool,
    match: str | None,
    fixtures_dir: Path = DEFAULT_FIXTURES_DIR,
) -> list[Fixture]:
    selectors = sum([bool(names), all_, match is not None])
    if selectors == 0:
        avail = _list_fixtures(fixtures_dir)
        sys.stderr.write(
            "no fixture selected. available: " + ", ".join(avail or ["(none)"]) + "\n"
        )
        raise SystemExit(2)
    if selectors > 1:
        sys.stderr.write(
            "--all, --match, and positional names are mutually exclusive\n"
        )
        raise SystemExit(2)

    if all_:
        chosen = _list_fixtures(fixtures_dir)
    elif match is not None:
        chosen = [n for n in _list_fixtures(fixtures_dir) if fnmatch.fnmatch(n, match)]
    else:
        chosen = list(names)

    out: list[Fixture] = []
    avail = set(_list_fixtures(fixtures_dir))
    for n in chosen:
        if n not in avail:
            raise SystemExit(f"fixture not found: {n!r} (in {fixtures_dir})")
        out.append(load_fixture(fixtures_dir / f"{n}.toml"))
    return out


def _git_short_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True
        ).strip()
    except Exception:
        return "unknown"


def _git_branch() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"], text=True
        ).strip()
    except Exception:
        return "unknown"


def _flatten_config(d: dict, prefix: str = "config") -> dict:
    out = {}
    for k, v in d.items():
        key = f"{prefix}.{k}"
        if isinstance(v, dict):
            out.update(_flatten_config(v, key))
        else:
            out[key] = v
    return out


def _materialize_timsseek_config(config_dict: dict, target_path: Path) -> None:
    """Write the [config] subtree to a JSON file (timsseek accepts both
    TOML and JSON — JSON is simpler to write here)."""
    target_path.parent.mkdir(parents=True, exist_ok=True)
    target_path.write_text(json.dumps(config_dict))


def run_one(
    fx: Fixture,
    out_root: Path = DEFAULT_OUT_ROOT,
    notes: str | None = None,
    dry_run: bool = False,
) -> None:
    """Execute one fixture: timsseek + wandb logging + optional entrapment."""
    ts = datetime.now().strftime("%Y%m%d-%H%M%S")
    run_dir = out_root / "logs" / f"{fx.name}-{ts}"
    res_dir = run_dir / "res"

    # The raw stem is what timsseek uses for its results subdir.
    raw_stem = Path(fx.inputs.raw).name
    if raw_stem.endswith(".d") or raw_stem.endswith(".idx"):
        raw_stem = Path(raw_stem).stem
    elif raw_stem.endswith(".d.tar"):
        raw_stem = raw_stem[: -len(".d.tar")]

    plan_msg = f"fixture={fx.name} run_dir={run_dir} raw_stem={raw_stem}"
    logger.info("plan: {}", plan_msg)
    if dry_run:
        return

    run_dir.mkdir(parents=True, exist_ok=True)
    res_dir.mkdir(parents=True, exist_ok=True)
    config_path = run_dir / "config.json"
    _materialize_timsseek_config(fx.config, config_path)

    sha = _git_short_sha()
    branch = _git_branch()
    tags = [fx.name, branch]
    if fx.has_entrapment():
        tags.append("entrapment")

    wandb_config = {
        "fixture": fx.name,
        "git.sha": sha,
        "git.branch": branch,
        "host": os.uname().nodename,
        **{
            f"inputs.{k}": v
            for k, v in fx.inputs.model_dump().items()
            if v is not None
        },
        **_flatten_config(fx.config),
    }
    run = wandb.init(
        entity=ENTITY,
        project=PROJECT,
        name=f"{fx.name}-{sha}",
        tags=tags,
        notes=notes,
        config=wandb_config,
    )
    try:
        cmd = [
            "cargo", "run", "--release", "--bin", "timsseek", "--",
            "--overwrite",
            "--config", str(config_path),
            "--speclib-file", fx.inputs.speclib,
            "--output-dir", str(res_dir),
            "--dotd-files", fx.inputs.raw,
        ]
        if fx.has_calibration_speclib():
            assert fx.inputs.calibration_speclib is not None
            cmd.extend(["--calib-lib", fx.inputs.calibration_speclib])
        stdout_log = run_dir / "timsseek_stdout.log"
        stderr_log = run_dir / "timsseek_stderr.log"
        logger.info("$ {}", " ".join(cmd))
        with stdout_log.open("w") as so, stderr_log.open("w") as se:
            subprocess.run(cmd, stdout=so, stderr=se, check=True)

        perf = res_dir / raw_stem / "performance_report.json"
        if perf.exists():
            run.log(json.loads(perf.read_text()))
        else:
            logger.warning("performance_report.json missing at {}", perf)

        if fx.has_entrapment():
            assert fx.inputs.entrapment_fasta is not None
            with tempfile.TemporaryDirectory() as td:
                target_local = Path(td) / "target.fasta"
                entrap_local = Path(td) / "entrap.fasta"
                s3_download_file(fx.inputs.fasta, str(target_local))
                s3_download_file(fx.inputs.entrapment_fasta, str(entrap_local))
                results_parquet = res_dir / raw_stem / "results.parquet"
                out_parquet = (
                    out_root / "parquets" / f"{fx.name}-{ts}-classified.parquet"
                )
                out_plot = out_root / "plots" / f"{fx.name}-fdr_curve-{ts}.png"
                scalars = analyse(
                    results_parquet=results_parquet,
                    target_fasta=target_local,
                    entrapment_fasta=entrap_local,
                    out_parquet=out_parquet,
                    out_plot=out_plot,
                    title=f"{fx.name} entrapment FDR",
                )
            run.log(scalars)
            run.log({"entrap/fdr_curve": wandb.Image(str(out_plot))})
            artifact = wandb.Artifact(f"{fx.name}-classified", type="dataset")
            artifact.add_file(str(out_parquet))
            run.log_artifact(artifact)
    finally:
        run.finish()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fixtures", nargs="*", help="Fixture names to run")
    p.add_argument("--all", dest="all_", action="store_true")
    p.add_argument("--match", help="Glob pattern over fixture names")
    p.add_argument("--notes", help="Free-form note added to wandb run")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    fixtures = select_fixtures(args.fixtures, args.all_, args.match)
    for fx in fixtures:
        run_one(fx, notes=args.notes, dry_run=args.dry_run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
