"""Build a fixture and push its inputs to S3.

Resolves polymorphic --db specs into concatenated FASTAs, uploads them and
the raw .d directory, builds the speclib via speclib_build_cli, and writes
the fixture TOML to bench/fixtures/<name>.toml.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from loguru import logger


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
        help=(
            "Target FASTA source (repeatable). Accepts: path/to/*.fasta(.gz),"
            " path/to/ids.txt, s3://..., UPxxxxxxxxx, accession"
        ),
    )
    p.add_argument(
        "--raw", required=True, help="Raw .d / .idx (local dir or s3://...)"
    )
    p.add_argument(
        "--config", required=True, help="Local timsseek config TOML to embed"
    )

    p.add_argument(
        "--entrap-db",
        action="append",
        default=[],
        metavar="SPEC",
        help="Entrapment FASTA source (repeatable)",
    )
    p.add_argument(
        "--calib-db",
        action="append",
        default=[],
        metavar="SPEC",
        help="Calibration FASTA source (repeatable)",
    )

    p.add_argument(
        "--speclib",
        help="If set, skip main speclib build and reference this URI",
    )
    p.add_argument(
        "--calibration-speclib",
        help="If set, skip calib speclib build and reference this URI",
    )

    p.add_argument("--koina-url", help="Koina URL passed to speclib_build_cli")

    p.add_argument(
        "--dry-run", action="store_true", help="Print the resolved plan and exit"
    )
    p.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing fixture TOML"
    )
    return p.parse_args(argv)


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
    """Emit a fixture TOML body as a string. The body is valid against
    bench._fixture_schema.Fixture."""
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
    config_text = config_path.read_text().strip()
    lines.append("# === embedded timsseek config ===")
    # Re-emit each non-empty top-level section under `[config.<section>]`. The
    # simplest approach: rewrite section headers `[X]` → `[config.X]` in the
    # source. Sub-sections `[X.Y]` → `[config.X.Y]`. This preserves comments
    # and ordering of the user's config.
    for raw_line in config_text.splitlines():
        is_section = (
            raw_line.startswith("[")
            and raw_line.rstrip().endswith("]")
            and not raw_line.startswith("[[")
        )
        if is_section:
            inner = raw_line.strip()[1:-1]
            lines.append(f"[config.{inner}]")
        else:
            lines.append(raw_line)
    lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    logger.info("push_fixture.py invoked: {}", vars(args))
    raise NotImplementedError("upload + build pipeline lands in Task 10")


if __name__ == "__main__":
    sys.exit(main())
