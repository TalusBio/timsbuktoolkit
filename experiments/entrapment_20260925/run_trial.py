"""Run one fixed HeLa entrapment search and record its empirical FDP metric."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

from bench.entrapment import analyze

from .evaluate import evaluate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--name", required=True)
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--root", type=Path, default=Path("shitshit/hela_entrapment"))
    parser.add_argument(
        "--timsseek", type=Path, default=Path("target/release/timsseek")
    )
    parser.add_argument(
        "--rescore-model", choices=("mlp", "lda", "gbm", "hybrid"), default="mlp"
    )
    args = parser.parse_args()

    trial = args.root / "trials" / args.name
    if trial.exists():
        raise FileExistsError(trial)
    trial.mkdir(parents=True)
    search = trial / "search"
    fdp = trial / "fdp"
    binary_hash = hashlib.sha256(args.timsseek.read_bytes()).hexdigest()
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], check=True, capture_output=True, text=True
    ).stdout.strip()
    command = [
        str(args.timsseek),
        "--speclib-uri",
        str(args.root / "entrapment.mzspeclib.txt.gz"),
        "--raw-inputs",
        str(args.raw),
        "--output-uri",
        str(search),
        "--max-qvalue",
        "1",
        "--decoy-strategy",
        "never",
        "--rescore-model",
        args.rescore_model,
    ]
    subprocess.run(command, check=True)
    report_path = analyze(
        search / args.raw.name / "results.parquet",
        args.root / "peptide_pairs.tsv",
        fdp,
    )
    with report_path.open(newline="") as stream:
        metric = evaluate(csv.DictReader(stream, delimiter="\t"))
    metric.update({
        "trial": args.name,
        "git_commit": commit,
        "binary_sha256": binary_hash,
        "rescore_model": args.rescore_model,
        "report": str(report_path),
    })
    rendered = json.dumps(metric, indent=2) + "\n"
    result = Path(__file__).parent / "results" / f"{args.name}.json"
    result.write_text(rendered)
    sys.stdout.write(rendered)


if __name__ == "__main__":
    main()
