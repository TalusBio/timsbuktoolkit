"""Bench runner: load named fixtures, run timsseek, log to wandb."""

from __future__ import annotations

import argparse
import fnmatch
import sys
from pathlib import Path

from loguru import logger

from bench._fixture_schema import Fixture, load_fixture

DEFAULT_FIXTURES_DIR = Path("bench/fixtures")
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
    """Resolve CLI selection flags into a list of loaded fixtures."""
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


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fixtures", nargs="*", help="Fixture names to run")
    p.add_argument("--all", dest="all_", action="store_true")
    p.add_argument("--match", help="Glob pattern over fixture names")
    p.add_argument("--notes", help="Free-form note added to wandb run")
    p.add_argument("--dry-run", action="store_true", help="Print plan and stop")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    fixtures = select_fixtures(args.fixtures, args.all_, args.match)
    for fx in fixtures:
        logger.info("would run fixture: {}", fx.name)
    raise NotImplementedError("timsseek invocation lands in Task 12")


if __name__ == "__main__":
    sys.exit(main())
