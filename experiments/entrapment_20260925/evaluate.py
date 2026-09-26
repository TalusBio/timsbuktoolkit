"""Score an entrapment run by original-target IDs at empirical paired FDP <= 1%."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path


def evaluate(rows, fdp_limit: float = 0.01, reported_q: float = 0.01):
    best = None
    at_reported_q = None
    last_q = -1.0
    n_rows = 0
    for row in rows:
        q = float(row["q_value"])
        fdp = float(row["paired_fdp"])
        targets = int(row["n_targets"])
        entrapments = int(row["n_entrapments"])
        if not math.isfinite(q) or not math.isfinite(fdp) or q < last_q:
            raise ValueError("report must have finite FDP and ascending q-values")
        last_q = q
        n_rows += 1
        point = {
            "q_value": q,
            "paired_fdp": fdp,
            "n_targets": targets,
            "n_entrapments": entrapments,
        }
        if q <= reported_q:
            at_reported_q = point
        if fdp <= fdp_limit and (
            best is None
            or targets > best["n_targets"]
            or (targets == best["n_targets"] and fdp < best["paired_fdp"])
        ):
            best = point
    if n_rows == 0:
        raise ValueError("report has no hits")
    return {
        "fdp_limit": fdp_limit,
        "reported_q": reported_q,
        "n_rows": n_rows,
        "best_at_empirical_fdp": best,
        "at_reported_q": at_reported_q,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    with args.report.open(newline="") as stream:
        result = evaluate(csv.DictReader(stream, delimiter="\t"))
    result["report"] = str(args.report)
    rendered = json.dumps(result, indent=2) + "\n"
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(rendered)
    sys.stdout.write(rendered)


if __name__ == "__main__":
    main()
