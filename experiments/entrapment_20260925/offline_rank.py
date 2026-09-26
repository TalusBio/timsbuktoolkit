"""Try score blends on saved competed candidates without rerunning extraction."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import polars as pl

from bench.entrapment import load_pairs, stripped_sequence


def assign_qvalues(scores: np.ndarray, is_target: np.ndarray) -> np.ndarray:
    """Match timsseek's pseudocount, reverse cumulative minimum, and score ties."""
    order = np.argsort(-scores, kind="stable")
    ranked_scores = scores[order]
    targets = np.cumsum(is_target[order], dtype=np.int64)
    decoys = np.cumsum(~is_target[order], dtype=np.int64) + 1
    raw = np.divide(
        decoys.astype(np.float32),
        targets.astype(np.float32),
        out=np.full(len(scores), np.inf, dtype=np.float32),
        where=targets > 0,
    )
    ranked_q = np.minimum.accumulate(np.minimum(raw, 1.0)[::-1])[::-1]
    starts = np.empty(len(scores), dtype=bool)
    starts[0] = True
    starts[1:] = ranked_scores[1:] != ranked_scores[:-1]
    first = np.maximum.accumulate(np.where(starts, np.arange(len(scores)), 0))
    ranked_q = ranked_q[first]
    result = np.empty_like(ranked_q)
    result[order] = ranked_q
    return result


def paired_metric(hits, labels, pairing, limit=0.01):
    """Evaluate FDRBench's one-fold paired estimator at complete q-value ties."""
    hits.sort(key=lambda h: (h[2], -h[3], h[0], h[1]))
    target_rank = {
        (peptide, charge): rank
        for rank, (peptide, charge, _, _) in enumerate(hits)
        if labels[peptide] == "target"
    }
    target_to_entrapment = {target: entrap for entrap, target in pairing.items()}
    seen_entrapments = set()
    targets = entrapments = orphaned = outranked = 0
    best = None
    at_reported_q = None
    start = 0
    while start < len(hits):
        q_value = hits[start][2]
        end = start + 1
        while end < len(hits) and hits[end][2] == q_value:
            end += 1
        for rank in range(start, end):
            peptide, charge, _, _ = hits[rank]
            if labels[peptide] == "target":
                targets += 1
                mate = target_to_entrapment.get(peptide)
                if (mate, charge) in seen_entrapments:
                    orphaned -= 1
                    outranked += 1
            else:
                entrapments += 1
                seen_entrapments.add((peptide, charge))
                t_rank = target_rank.get((pairing[peptide], charge))
                if t_rank is None or rank < t_rank:
                    orphaned += 1
        fdp = (entrapments + orphaned + 2 * outranked) / (targets + entrapments)
        point = {
            "q_value": float(q_value),
            "paired_fdp": fdp,
            "n_targets": targets,
            "n_entrapments": entrapments,
        }
        if q_value <= 0.01:
            at_reported_q = point
        if fdp <= limit and (best is None or targets > best["n_targets"]):
            best = point
        start = end
    return {"best_at_empirical_fdp": best, "at_reported_q": at_reported_q}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--alpha", type=float, nargs="+", required=True)
    args = parser.parse_args()

    labels, pairing = load_pairs(args.pairs)
    frame = pl.read_parquet(
        args.results,
        columns=[
            "sequence",
            "precursor_charge",
            "is_target",
            "discriminant_score",
            "main_score",
            "qvalue",
        ],
    )
    is_target = frame["is_target"].to_numpy()
    mlp = frame["discriminant_score"].to_numpy()
    main = frame["main_score"].to_numpy()
    main_log = np.log1p(main)
    target_rows = np.flatnonzero(is_target)
    peptides = [stripped_sequence(seq) for seq in frame["sequence"][target_rows]]
    if any(peptide not in labels for peptide in peptides):
        raise ValueError("target result is absent from the pair file")
    base = pl.DataFrame({
        "peptide": peptides,
        "charge": frame["precursor_charge"][target_rows],
    })

    for alpha in args.alpha:
        scores = (mlp + np.float32(alpha) * main_log).astype(np.float32)
        q_values = assign_qvalues(scores, is_target)
        if alpha == 0:
            stored = frame["qvalue"].to_numpy()
            mismatch = np.count_nonzero(q_values != stored)
            sys.stderr.write(f"baseline q-value mismatches: {mismatch}\n")
            if mismatch:
                raise ValueError("offline q-value calculation disagrees with timsseek")
        unique = (
            base
            .with_columns(
                pl.Series("q_value", q_values[target_rows]),
                pl.Series("score", scores[target_rows]),
            )
            .sort(
                ["peptide", "charge", "q_value", "score"],
                descending=[False, False, False, True],
            )
            .unique(subset=["peptide", "charge"], keep="first", maintain_order=True)
        )
        hits = unique.select("peptide", "charge", "q_value", "score").iter_rows()
        result = paired_metric(list(hits), labels, pairing)
        result["alpha"] = alpha
        sys.stdout.write(json.dumps(result) + "\n")


if __name__ == "__main__":
    main()
