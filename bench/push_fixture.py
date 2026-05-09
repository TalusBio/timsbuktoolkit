"""Build a fixture and push its inputs to S3.

Resolves polymorphic --db specs into peptide lists, uploads them and the raw
.d directory, builds the speclib via speclib_build_cli --peptides, and writes
the fixture TOML to bench/fixtures/<name>.toml.
"""

from __future__ import annotations

import argparse
import random
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

from loguru import logger

from bench._db_resolver import resolve_dbs
from bench._digest import digest_proteins, length_filter, parse_fasta
from bench._s3 import s3_upload_dir, s3_upload_file
from bench._shuffle import generate_shuffled_entrapment


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
        help="Target proteome source (repeatable)",
    )
    p.add_argument("--raw", required=True, help="Raw .d / .idx (local dir or s3://...)")
    p.add_argument(
        "--config", required=True, help="Local timsseek config TOML to embed"
    )
    p.add_argument(
        "--entrap-db",
        action="append",
        default=[],
        metavar="SPEC",
        help="Foreign entrapment db specs, OR exactly 'SHUFFLED' for Algorithm 1",
    )
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
    p.add_argument(
        "--entrap-ratio",
        type=float,
        default=1.0,
        help="Entrapment ratio r >= 1.0 (default 1.0)",
    )
    p.add_argument(
        "--peptide-min-len",
        type=int,
        default=7,
        help="Minimum peptide length after digestion (default 7)",
    )
    p.add_argument(
        "--peptide-max-len",
        type=int,
        default=30,
        help="Maximum peptide length after digestion (default 30)",
    )
    p.add_argument(
        "--missed-cleavages",
        type=int,
        default=1,
        help="Number of missed cleavages for trypsin digestion (default 1)",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=42,
        help="RNG seed for shuffle / subsample (default 42)",
    )
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument(
        "--force",
        action="store_true",
        help="Re-upload S3 objects even if they already exist",
    )

    args = p.parse_args(argv)

    if args.entrap_ratio < 1.0:
        p.error("--entrap-ratio must be >= 1.0")

    if len(args.entrap_db) > 1 and "SHUFFLED" in args.entrap_db:
        p.error("SHUFFLED cannot be mixed with other --entrap-db specs")

    return args


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _subsample_set(s: set[str], k: int, seed: int) -> set[str]:
    rng = random.Random(seed)
    if k >= len(s):
        return set(s)
    return set(rng.sample(sorted(s), k))


def _write_peptides(peptides: set[str], path: Path) -> None:
    path.write_text("\n".join(sorted(peptides)) + "\n")


def _write_pairing(pairs: list[tuple[str, str]], path: Path) -> None:
    lines = ["target_peptide\tentrap_peptide"]
    for t, s in pairs:
        lines.append(f"{t}\t{s}")
    path.write_text("\n".join(lines) + "\n")


def _digest_fasta(
    fasta_path: Path,
    missed_cleavages: int,
    min_len: int,
    max_len: int,
) -> set[str]:
    proteins = parse_fasta(fasta_path)
    raw = digest_proteins(proteins, missed_cleavages=missed_cleavages)
    return length_filter(raw, min_len=min_len, max_len=max_len)


# ---------------------------------------------------------------------------
# Speclib build
# ---------------------------------------------------------------------------


def run_speclib_build(
    peptides_uri: str,
    speclib_uri: str,
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
        "--peptides",
        peptides_uri,
        "--fixed-mod",
        "C[U:4]",
        "--max-ions",
        "10",
        "-o",
        speclib_uri,
    ]
    if koina_url:
        cmd.extend(["--koina-url", koina_url])
    cmd.extend(["--request-delay-ms", str(request_delay_ms)])
    logger.info("$ {}", " ".join(cmd))
    subprocess.run(cmd, check=True)


# ---------------------------------------------------------------------------
# Fixture TOML builder
# ---------------------------------------------------------------------------


def build_fixture_toml(
    name: str,
    description: str,
    config_path: Path,
    target_peptides_uri: str,
    speclib_uri: str,
    raw_uri: str,
    entrapment_peptides_uri: str | None = None,
    entrapment_ratio: float | None = None,
    entrapment_mode: str | None = None,
    pairing_uri: str | None = None,
    calibration_speclib_uri: str | None = None,
) -> str:
    lines: list[str] = []
    lines.append(f'name = "{name}"')
    desc = description.replace('"', '\\"')
    lines.append(f'description = "{desc}"')
    lines.append("")
    lines.append("[inputs]")
    lines.append(f'target_peptides = "{target_peptides_uri}"')
    lines.append(f'speclib = "{speclib_uri}"')
    lines.append(f'raw = "{raw_uri}"')
    if entrapment_peptides_uri is not None:
        lines.append(f'entrapment_peptides = "{entrapment_peptides_uri}"')
    if entrapment_ratio is not None:
        lines.append(f"entrapment_ratio = {entrapment_ratio}")
    if entrapment_mode is not None:
        lines.append(f'entrapment_mode = "{entrapment_mode}"')
    if pairing_uri is not None:
        lines.append(f'pairing = "{pairing_uri}"')
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


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------


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
    entrap_ratio: float = 1.0,
    peptide_min_len: int = 7,
    peptide_max_len: int = 30,
    missed_cleavages: int = 1,
    seed: int = 42,
) -> None:
    """Execute the full upload + build + write-toml flow."""
    # Validate SHUFFLED mixing up-front (defensive; parse_args also checks)
    if len(entrap_db) > 1 and "SHUFFLED" in entrap_db:
        raise ValueError("SHUFFLED cannot be mixed with other --entrap-db specs")
    if entrap_ratio < 1.0:
        raise ValueError("entrap_ratio must be >= 1.0")

    dest_prefix = f"s3://{bucket}/{prefix.rstrip('/')}/{name}"

    main_speclib_uri = speclib_uri or f"{dest_prefix}/lib.msgpack.zst"
    final_calib_speclib_uri: str | None = calibration_speclib_uri
    if final_calib_speclib_uri is None and calib_db:
        final_calib_speclib_uri = f"{dest_prefix}/calib_lib.msgpack.zst"

    raw_uri = raw if raw.startswith("s3://") else f"{dest_prefix}/sample.d"

    if fixture_target.exists() and not overwrite and not dry_run:
        raise FileExistsError(
            f"fixture TOML already exists: {fixture_target}"
            " (pass --overwrite to replace)"
        )

    use_shuffled = entrap_db == ["SHUFFLED"]
    use_foreign = bool(entrap_db) and not use_shuffled

    plan = {
        "name": name,
        "dest_prefix": dest_prefix,
        "use_shuffled": use_shuffled,
        "use_foreign": use_foreign,
        "entrap_ratio": entrap_ratio,
        "raw_uri": raw_uri,
        "main_speclib_uri": main_speclib_uri,
        "calib_speclib_uri": final_calib_speclib_uri,
        "build_main_speclib": speclib_uri is None,
        "build_calib_speclib": (bool(calib_db) and calibration_speclib_uri is None),
    }
    logger.info("plan: {}", plan)
    if dry_run:
        logger.info("--dry-run: stopping before any side effects")
        return

    with tempfile.TemporaryDirectory() as td:
        workdir = Path(td)

        # 1. Resolve target proteome → target peptides
        target_fasta = workdir / "target.fasta"
        resolve_dbs(db, target_fasta)
        p_target = _digest_fasta(
            target_fasta,
            missed_cleavages=missed_cleavages,
            min_len=peptide_min_len,
            max_len=peptide_max_len,
        )

        # 2. Entrapment handling
        p_foreign: set[str] = set()
        pairs: list[tuple[str, str]] = []
        actual_r: float | None = None
        mode: str | None = None
        emit_pairing = False

        if use_shuffled:
            # Algorithm 1: shuffled entrapment
            r_int = int(entrap_ratio)
            pairs = generate_shuffled_entrapment(p_target, r=r_int, seed=seed)
            kept_targets = {t for (t, _) in pairs}
            p_target = kept_targets
            p_foreign = {s for (_, s) in pairs}
            actual_r = float(r_int)
            mode = "shuffled"
            emit_pairing = r_int == 1

        elif use_foreign:
            # Algorithm 2: foreign entrapment
            entrap_fasta = workdir / "entrap.fasta"
            resolve_dbs(entrap_db, entrap_fasta)
            all_foreign = _digest_fasta(
                entrap_fasta,
                missed_cleavages=missed_cleavages,
                min_len=peptide_min_len,
                max_len=peptide_max_len,
            )
            # Remove any peptides that appear in the target
            all_foreign -= p_target
            n_needed = int(entrap_ratio * len(p_target))
            if n_needed <= len(all_foreign):
                p_foreign = _subsample_set(all_foreign, n_needed, seed)
                actual_r = entrap_ratio
            else:
                p_foreign = all_foreign
                actual_r = len(p_foreign) / len(p_target) if p_target else 0.0
                warnings.warn(
                    f"Not enough foreign peptides: needed {n_needed}, "
                    f"got {len(p_foreign)}. Actual r = {actual_r:.4f}",
                    stacklevel=2,
                )
            pairs = []
            mode = "foreign"
            emit_pairing = False

        # 3. Build database peptide list (target ∪ entrap)
        p_database = p_target | p_foreign

        # 4. Write local peptide files
        target_pep_local = workdir / "target.peptides.txt"
        database_pep_local = workdir / "database.peptides.txt"
        _write_peptides(p_target, target_pep_local)
        _write_peptides(p_database, database_pep_local)

        entrap_pep_local: Path | None = None
        pairing_local: Path | None = None
        if entrap_db:
            entrap_pep_local = workdir / "entrap.peptides.txt"
            _write_peptides(p_foreign, entrap_pep_local)
        if emit_pairing:
            pairing_local = workdir / "pairing.tsv"
            _write_pairing(pairs, pairing_local)

        # 5. Upload peptide files
        target_pep_uri = f"{dest_prefix}/target.peptides.txt"
        database_pep_uri = f"{dest_prefix}/database.peptides.txt"
        s3_upload_file(str(target_pep_local), target_pep_uri, skip_if_exists=not force)
        s3_upload_file(
            str(database_pep_local), database_pep_uri, skip_if_exists=not force
        )

        entrap_pep_uri: str | None = None
        pairing_uri: str | None = None
        if entrap_pep_local is not None:
            entrap_pep_uri = f"{dest_prefix}/entrap.peptides.txt"
            s3_upload_file(
                str(entrap_pep_local), entrap_pep_uri, skip_if_exists=not force
            )
        if pairing_local is not None:
            pairing_uri = f"{dest_prefix}/pairing.tsv"
            s3_upload_file(str(pairing_local), pairing_uri, skip_if_exists=not force)

        # 6. Upload raw if local
        if not raw.startswith("s3://"):
            s3_upload_dir(raw, raw_uri, idempotent=not force)

        # 7. Calibration db peptides (Algorithm 2 path, no entrapment subtract)
        calib_pep_uri: str | None = None
        if calib_db:
            calib_fasta = workdir / "calib.fasta"
            resolve_dbs(calib_db, calib_fasta)
            p_calib = _digest_fasta(
                calib_fasta,
                missed_cleavages=missed_cleavages,
                min_len=peptide_min_len,
                max_len=peptide_max_len,
            )
            calib_pep_local = workdir / "calib.peptides.txt"
            _write_peptides(p_calib, calib_pep_local)
            calib_pep_uri = f"{dest_prefix}/calib.peptides.txt"
            s3_upload_file(
                str(calib_pep_local), calib_pep_uri, skip_if_exists=not force
            )

        # 8. Build speclib(s)
        if speclib_uri is None:
            run_speclib_build(
                database_pep_uri, main_speclib_uri, koina_url, request_delay_ms
            )
        if calib_db and calibration_speclib_uri is None:
            assert calib_pep_uri is not None
            assert final_calib_speclib_uri is not None
            run_speclib_build(
                calib_pep_uri,
                final_calib_speclib_uri,
                koina_url,
                request_delay_ms,
            )

    # 9. Emit fixture TOML
    body = build_fixture_toml(
        name=name,
        description="",
        config_path=Path(config),
        target_peptides_uri=target_pep_uri,
        speclib_uri=main_speclib_uri,
        raw_uri=raw_uri,
        entrapment_peptides_uri=entrap_pep_uri,
        entrapment_ratio=actual_r,
        entrapment_mode=mode,
        pairing_uri=pairing_uri,
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
        entrap_ratio=args.entrap_ratio,
        peptide_min_len=args.peptide_min_len,
        peptide_max_len=args.peptide_max_len,
        missed_cleavages=args.missed_cleavages,
        seed=args.seed,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
