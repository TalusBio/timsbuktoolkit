import argparse
import logging
from pathlib import Path

import polars as pl
from rich.pretty import pprint

from .backends.mlp import mlp_stuff
from .backends.mokapot import mokapot_stuff
from .backends.xgb import xgboost_stuff
from .feateng import read_files, to_mokapot_df
from .plotting import main_score_hist


# TODO: Rename
def main(args):
    paths = [Path(p) for p in args.results_dir]
    for p in paths:
        if not p.exists():
            raise FileNotFoundError(f"Path {p} does not exist")
    data = read_files(paths)
    data = data.filter(pl.col("obs_mobility").is_not_nan())

    outdir = Path(args.output_dir)
    if not outdir.exists():
        outdir.mkdir(parents=True)

    main_score_hist(data, outdir)
    xgboost_stuff(data, outdir)
    data, cols = to_mokapot_df(data)
    score = mlp_stuff(data, cols=cols, output_dir=outdir)
    mokapot_stuff(data, outdir)

    exit()
    # Leaving the code I used for hparam tuning of the MLP
    # (0.2, 4, 512, 0)
    # (0.05, 3, 64, 0.1)
    scores = {}
    best_score = score
    best_params = None
    combs = []
    for pw in [0.05, 0.2, 0.5, 0.8]:
        # for do in [0, 0.05, 0.1, 0.2, 0.5]:
        for do in [0, 0.05, 0.1]:
            for hl in [3, 4, 5]:
                for hd in [64, 128]:
                    combs.append((hl, hd, pw, do))
            for hl in [1, 2, 3, 4, 5]:
                for hd in [256, 512]:
                    combs.append((hl, hd, pw, do))

    # shuffle combs
    import random

    random.shuffle(combs)
    for hl, hd, pw, do in combs:
        tmp_outdir = outdir / f"hl-{hl}_hd-{hd}-pw-{pw}-do-{do}"
        tmp_outdir.mkdir(parents=True, exist_ok=True)
        score = mlp_stuff(
            data,
            cols=cols,
            output_dir=tmp_outdir,
            pos_weight=pw,
            nhidden_layers=hl,
            hidden_dims=hd,
            dropout=do,
        )
        scores[(hl, hd, do)] = score
        if score > best_score:
            pprint(f"New best score: {score} at {pw} {hl} {hd} {do}")
            best_score = score
            best_params = (pw, hl, hd, do)
        else:
            pprint("Keeping old score")
            pprint(f"Current score {score}")
            pprint(
                f"Old best score: {best_score} at {best_params if best_params is not None else 'BASELINE'}"
            )

    sorted_scores = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    pprint(sorted_scores)
    pprint(best_params)


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--results_dir",
        type=str,
        nargs="+",
        help="Path to the directories containing the results.csv files",
    )
    parser.add_argument(
        "-o", "--output_dir", type=str, default=".", help="Output directory"
    )
    return parser


def cli_main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    parser = build_parser()
    args, unkargs = parser.parse_known_args()
    if unkargs:
        raise ValueError(f"Unknown arguments: {unkargs}")

    main(args)


if __name__ == "__main__":
    cli_main()
