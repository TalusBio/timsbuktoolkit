"""Build a fixture and push its inputs to S3.

Resolves polymorphic --db specs into concatenated FASTAs, uploads them and
the raw .d directory, builds the speclib via speclib_build_cli, and writes
the fixture TOML to bench/fixtures/<name>.toml.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

from loguru import logger

from bench._db_resolver import resolve_dbs
from bench._s3 import s3_upload_dir, s3_upload_file


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--name", required=True, help="Fixture name (used as filename + S3 subdir)"
    )
    p.add_argument("--bucket", required=True, help="S3 bucket")
    p.add_argument("--prefix", required=True, help="S3 prefix under the bucket")
    p.add_argument(
        "--db",
        action="append",
        default=[],
        required=True,
        metavar="SPEC",
        help="Target FASTA source (repeatable)",
    )
    p.add_argument("--raw", required=True, help="Raw .d / .idx (local dir or s3://...)")
    p.add_argument(
        "--config", required=True, help="Local timsseek config TOML to embed"
    )
    p.add_argument("--entrap-db", action="append", default=[], metavar="SPEC")
    p.add_argument("--calib-db", action="append", default=[], metavar="SPEC")
    p.add_argument(
        "--speclib",
        dest="speclib_uri",
        help="Skip main speclib build, reference this URI",
    )
    p.add_argument(
        "--calibration-speclib",
        dest="calibration_speclib_uri",
        help="Skip calib speclib build, reference this URI",
    )
    p.add_argument("--koina-url")
    p.add_argument(
        "--request-delay-ms",
        type=int,
        default=500,
        help="Per-request delay passed to speclib_build_cli (ms; default 500)",
    )
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument(
        "--force",
        action="store_true",
        help="Re-upload S3 objects even if they already exist",
    )
    return p.parse_args(argv)


def run_speclib_build(
    fasta_s3: str,
    speclib_s3: str,
    koina_url: str | None,
    request_delay_ms: int = 500,
) -> None:
    cmd = [
        "cargo",
        "run",
        "--release",
        "-p",
        "speclib_build_cli",
        "--",
        "--fasta",
        fasta_s3,
        "--fixed-mod",
        "C[U:4]",
        "--max-ions",
        "10",
        "-o",
        speclib_s3,
    ]
    if koina_url:
        cmd.extend(["--koina-url", koina_url])
    cmd.extend(["--request-delay-ms", str(request_delay_ms)])
    logger.info("$ {}", " ".join(cmd))
    subprocess.run(cmd, check=True)


def build_fixture_toml(
    name: str,
    description: str,
    config_path: Path,
    fasta_uri: str,
    speclib_uri: str,
    raw_uri: str,
    entrapment_fasta_uri: str | None,
    calibration_speclib_uri: str | None,
) -> str:
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
    config_text = config_path.read_text().strip()
    for raw_line in config_text.splitlines():
        is_section_header = (
            raw_line.startswith("[")
            and raw_line.rstrip().endswith("]")
            and not raw_line.startswith("[[")
        )
        if is_section_header:
            inner = raw_line.strip()[1:-1]
            lines.append(f"[config.{inner}]")
        else:
            lines.append(raw_line)
    lines.append("")
    return "\n".join(lines)


def _resolve_and_upload_fasta(
    specs: list[str],
    s3_dest: str,
    label: str,
    workdir: Path,
    skip_if_exists: bool = False,
) -> None:
    local = workdir / f"{label}.fasta"
    resolve_dbs(specs, local)
    s3_upload_file(str(local), s3_dest, skip_if_exists=skip_if_exists)


def run_pipeline(
    *,
    name: str,
    bucket: str,
    prefix: str,
    db: list[str],
    raw: str,
    config: str,
    entrap_db: list[str],
    calib_db: list[str],
    speclib_uri: str | None,
    calibration_speclib_uri: str | None,
    koina_url: str | None,
    fixture_target: Path,
    overwrite: bool,
    dry_run: bool,
    force: bool = False,
    request_delay_ms: int = 500,
) -> None:
    """Execute the full upload + build + write-toml flow."""
    dest_prefix = f"s3://{bucket}/{prefix.rstrip('/')}/{name}"
    target_fasta_uri = f"{dest_prefix}/proteome.fasta"
    entrap_fasta_uri = f"{dest_prefix}/entrap.fasta" if entrap_db else None
    calib_fasta_uri = f"{dest_prefix}/calib.fasta" if calib_db else None

    main_speclib_uri = speclib_uri or f"{dest_prefix}/lib.msgpack.zst"
    # When entrap_db is present, the speclib must cover both target+entrap so
    # the search can score entrapment peptides. We upload a concatenated fasta
    # to a separate URI and point speclib_build_cli at it. The per-fasta
    # target/entrap files are still uploaded separately so analyse() can
    # classify hits by source.
    speclib_input_fasta_uri = (
        f"{dest_prefix}/speclib_input.fasta" if entrap_db else None
    )
    final_calib_speclib_uri: str | None = calibration_speclib_uri
    if final_calib_speclib_uri is None and calib_db:
        final_calib_speclib_uri = f"{dest_prefix}/calib_lib.msgpack.zst"

    # Raw is either a local dir we upload or an existing s3 URI we just reference
    if raw.startswith("s3://"):
        raw_uri = raw
    else:
        raw_uri = f"{dest_prefix}/sample.d"

    if fixture_target.exists() and not overwrite and not dry_run:
        raise FileExistsError(
            f"fixture TOML already exists: {fixture_target}"
            " (pass --overwrite to replace)"
        )

    plan = {
        "name": name,
        "dest_prefix": dest_prefix,
        "target_fasta_uri": target_fasta_uri,
        "entrap_fasta_uri": entrap_fasta_uri,
        "calib_fasta_uri": calib_fasta_uri,
        "raw_uri": raw_uri,
        "main_speclib_uri": main_speclib_uri,
        "calib_speclib_uri": final_calib_speclib_uri,
        "build_main_speclib": speclib_uri is None,
        "build_calib_speclib": (calib_db != [] and calibration_speclib_uri is None),
    }
    logger.info("plan: {}", plan)
    if dry_run:
        logger.info("--dry-run: stopping before any side effects")
        return

    with tempfile.TemporaryDirectory() as td:
        workdir = Path(td)

        # 1. Resolve and upload target FASTA
        _resolve_and_upload_fasta(
            db, target_fasta_uri, "proteome", workdir, skip_if_exists=not force
        )

        # 2. Optional entrapment FASTA
        if entrap_db:
            assert entrap_fasta_uri is not None
            _resolve_and_upload_fasta(
                entrap_db, entrap_fasta_uri, "entrap", workdir, skip_if_exists=not force
            )

        # 3. Optional calibration FASTA
        if calib_db:
            assert calib_fasta_uri is not None
            _resolve_and_upload_fasta(
                calib_db, calib_fasta_uri, "calib", workdir, skip_if_exists=not force
            )

        # 4. Upload raw dir if local
        if not raw.startswith("s3://"):
            s3_upload_dir(raw, raw_uri, idempotent=not force)

        # 4b. If entrap_db present, build a merged target+entrap fasta locally
        # and upload it as the speclib build input.
        if speclib_input_fasta_uri is not None:
            merged_local = workdir / "speclib_input.fasta"
            target_local = workdir / "proteome.fasta"
            entrap_local = workdir / "entrap.fasta"
            with merged_local.open("wb") as out:
                for src in (target_local, entrap_local):
                    out.write(src.read_bytes())
                    if not src.read_bytes().endswith(b"\n"):
                        out.write(b"\n")
            s3_upload_file(
                str(merged_local),
                speclib_input_fasta_uri,
                skip_if_exists=not force,
            )

        # 5. Build speclib(s) if not user-provided
        if speclib_uri is None:
            speclib_input_uri = speclib_input_fasta_uri or target_fasta_uri
            run_speclib_build(
                speclib_input_uri, main_speclib_uri, koina_url, request_delay_ms
            )
        if calib_db and calibration_speclib_uri is None:
            assert calib_fasta_uri is not None
            assert final_calib_speclib_uri is not None
            run_speclib_build(
                calib_fasta_uri,
                final_calib_speclib_uri,
                koina_url,
                request_delay_ms,
            )

    # 6. Emit fixture TOML
    body = build_fixture_toml(
        name=name,
        description="",
        config_path=Path(config),
        fasta_uri=target_fasta_uri,
        speclib_uri=main_speclib_uri,
        raw_uri=raw_uri,
        entrapment_fasta_uri=entrap_fasta_uri,
        calibration_speclib_uri=final_calib_speclib_uri,
    )
    fixture_target.parent.mkdir(parents=True, exist_ok=True)
    fixture_target.write_text(body)
    logger.info("wrote fixture: {}", fixture_target)
    logger.info("remember to: git add {}", fixture_target)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    fixture_target = Path("bench/fixtures") / f"{args.name}.toml"
    run_pipeline(
        name=args.name,
        bucket=args.bucket,
        prefix=args.prefix,
        db=args.db,
        raw=args.raw,
        config=args.config,
        entrap_db=args.entrap_db,
        calib_db=args.calib_db,
        speclib_uri=args.speclib_uri,
        calibration_speclib_uri=args.calibration_speclib_uri,
        koina_url=args.koina_url,
        fixture_target=fixture_target,
        overwrite=args.overwrite,
        dry_run=args.dry_run,
        force=args.force,
        request_delay_ms=args.request_delay_ms,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
