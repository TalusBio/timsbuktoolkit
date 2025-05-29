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
    # mokapot_stuff(data, outdir)
    # xgboost_stuff(data, outdir)
    data, cols = to_mokapot_df(data)
    score = mlp_stuff(data, cols=cols, output_dir=outdir)


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
