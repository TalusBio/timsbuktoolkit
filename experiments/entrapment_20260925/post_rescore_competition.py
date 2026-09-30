"""Compete saved no-competition scores, then recalculate q-values and paired FDP."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import polars as pl

from bench.entrapment import analyze

from .evaluate import evaluate


def assign_qvalues(scores: np.ndarray, is_target: np.ndarray) -> np.ndarray:
    """Mirror timsseek's sorted-score q-value calculation, including score ties."""
    if len(scores) != len(is_target) or not np.isfinite(scores).all():
        raise ValueError("scores and labels must be aligned and finite")
    if np.any(scores[:-1] < scores[1:]):
        raise ValueError("results must be sorted by descending discriminant score")

    target_count = np.cumsum(is_target, dtype=np.int64)
    decoy_count = np.cumsum(~is_target, dtype=np.int64) + 1
    ratio = np.full(len(scores), np.inf, dtype=np.float32)
    np.divide(
        decoy_count.astype(np.float32),
        target_count.astype(np.float32),
        out=ratio,
        where=target_count != 0,
    )
    qvalues = np.minimum.accumulate(ratio[::-1])[::-1].copy()
    np.minimum(qvalues, np.float32(1), out=qvalues)

    # The engine gives every score tie the first member's q-value.
    starts = np.concatenate(([0], np.flatnonzero(scores[1:] != scores[:-1]) + 1))
    return np.repeat(qvalues[starts], np.diff(np.append(starts, len(scores))))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    if args.out.exists():
        raise FileExistsError(args.out)
    columns = [
        "library_id",
        "decoy_group_id",
        "sequence",
        "precursor_charge",
        "is_target",
        "discriminant_score",
        "qvalue",
    ]
    frame = pl.read_parquet(args.input, columns=columns)
    scores = frame["discriminant_score"].to_numpy()
    labels = frame["is_target"].to_numpy()
    replayed = assign_qvalues(scores, labels)
    if not np.array_equal(replayed, frame["qvalue"].to_numpy()):
        raise ValueError(
            "saved q-values differ from replay; input is not the E2 full output"
        )

    # E2's Parquet rows are already in descending score order. Keep the first
    # member of each (group, charge), so tied scores retain the engine's stable
    # preexisting rank order.
    winners = frame.unique(
        subset=["decoy_group_id", "precursor_charge"],
        keep="first",
        maintain_order=True,
    )
    winners = winners.with_columns(
        pl.Series(
            "qvalue",
            assign_qvalues(
                winners["discriminant_score"].to_numpy(),
                winners["is_target"].to_numpy(),
            ),
        )
    )
    args.out.mkdir(parents=True)
    results = args.out / "results.parquet"
    winners.write_parquet(results)
    report = analyze(results, args.pairs, args.out / "fdp")
    with report.open(newline="") as stream:
        metric = evaluate(csv.DictReader(stream, delimiter="\t"))
    metric.update({
        "input": str(args.input),
        "input_rows": frame.height,
        "winners": winners.height,
        "output": str(results),
        "report": str(report),
        "method": (
            "best discriminant per (decoy_group_id, precursor_charge), then recompute q"
        ),
    })
    rendered = json.dumps(metric, indent=2) + "\n"
    (args.out / "metric.json").write_text(rendered)
    sys.stdout.write(rendered)


if __name__ == "__main__":
    main()
