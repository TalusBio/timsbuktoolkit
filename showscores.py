# /// script
# dependencies = [
#   "polars",
#   "rich",
#   "matplotlib",
#   "numpy",
#   "tqdm",
#   "mokapot @ git+https://github.com/jspaezp/mokapot.git@feat/re_add_compound_index_spec",
#   "xgboost",
# ]
# ///

import argparse
import logging
from dataclasses import dataclass
from pathlib import Path
from xgboost import XGBClassifier
from sklearn.model_selection import GridSearchCV

import matplotlib.pyplot as plt
import mokapot
import numpy as np
import pandas as pd
import polars as pl
import xgboost as xgb
from mokapot.column_defs import ColumnGroups, OptionalColumns
from rich.pretty import pprint
from tqdm.auto import tqdm

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
)


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


def read_files(results_dirs: list[Path]) -> pl.LazyFrame:
    files = set()
    for results_dir in results_dirs:
        files.update(results_dir.glob("*.parquet"))

    files = list(files)
    pprint(f"Scanning {len(files)} files")
    data = pl.scan_parquet(files)
    pprint("Done scanning")
    return data


def lazy_abs_and_maxfill(df: pl.LazyFrame, columns: list[str]) -> pl.LazyFrame:
    exprs_first = []
    exprs_later = []
    for column in columns:
        exprs_first.append(pl.col(column).abs())
        exprs_later.append(pl.col(column).fill_nan(pl.col(column).max()))

    return df.with_columns(exprs_first).with_columns(exprs_later)


def lazy_zero_fill(df: pl.LazyFrame, columns: list[str]) -> pl.LazyFrame:
    exprs = []
    for column in columns:
        exprs.append(pl.col(column).fill_nan(0))
    return df.with_columns(exprs)


def log_cols(df: pl.LazyFrame, columns: list[str]) -> pl.LazyFrame:
    exprs = []
    for col in columns:
        exprs.append(pl.col(col).log1p().fill_nan(0))
    return df.with_columns(exprs)


def add_id(df: pl.LazyFrame) -> pl.LazyFrame:
    # Add an id column ... pretty simple incremental integer
    df = df.with_row_index(name="id")
    return df


def check_noninf(df: pl.DataFrame, columns: list[str]):
    pprint("Checking for infinite values")
    any_inf = False
    df_inf = df.select(columns).filter(pl.any_horizontal(pl.all().is_infinite()))
    pprint(df_inf)
    for col in columns:
        if df_inf[col].is_infinite().any():
            pprint(f"Column {col} has infinite values")
            any_inf = True
    if any_inf:
        raise ValueError("Data contains infinite values")


def check_nonnan(df: pl.LazyFrame, columns: list[str]):
    pprint("Checking for NaN values")
    nan_cols = []
    df_nan = df.select(columns).filter(pl.any_horizontal(pl.all().is_nan()))
    pprint(df_nan)
    for col in columns:
        if df_nan[col].is_nan().any():
            pprint(f"Column {col} has NaN values")
            nan_cols.append(col)
    if nan_cols:
        raise ValueError(f"Data contains NaN values: {nan_cols}")


def cast_f32(df: pl.LazyFrame, cols: list[str]) -> pl.LazyFrame:
    exprs = []
    for col in cols:
        exprs.append(pl.col(col).cast(pl.Float32))
    return df.with_columns(exprs)


def calculate_delta_rt(df: pl.LazyFrame) -> pl.LazyFrame:
    return df.with_columns(
        delta_rt=(pl.col("obs_rt_seconds") - pl.col("precursor_rt_query_seconds")),
    ).with_columns(
        abs_delta_rt=(
            pl.col("obs_rt_seconds") - pl.col("precursor_rt_query_seconds")
        ).abs(),
    )


def to_mokapot_df(df: pl.LazyFrame) -> tuple[pd.DataFrame, ColumnGroups]:
    pprint("Starting to_mokapot")
    loggable_cols = (
        # Log
        "npeaks",
        "lazyerscore",
        "ms1_summed_precursor_intensity",
        "main_score",
        "delta_next",
        "lazyerscore_vs_baseline",
        "norm_lazyerscore_vs_baseline",
    )
    imputable_cols = (
        # Abs impute
        "ms2_mz_error_0",
        "ms2_mz_error_1",
        "ms2_mz_error_2",
        "ms2_mz_error_3",
        "ms2_mz_error_4",
        "ms2_mz_error_5",
        "ms2_mz_error_6",
        "ms2_mobility_error_0",
        "ms2_mobility_error_1",
        "ms2_mobility_error_2",
        "ms2_mobility_error_3",
        "ms2_mobility_error_4",
        "ms2_mobility_error_5",
        "ms2_mobility_error_6",
        "ms1_mz_error_0",
        "ms1_mz_error_1",
        "ms1_mz_error_2",
        "ms1_mobility_error_0",
        "ms1_mobility_error_1",
        "ms1_mobility_error_2",
        "sq_delta_ms1_ms2_mobility",
        "delta_ms1_ms2_mobility",
    )
    # zero_imputable_cols = ("ms1_ms2_correlation",)
    zero_imputable_cols = ()
    feat_cols = (
        (
            "precursor_charge",
            "precursor_mz",
            "precursor_mobility_query",
            "precursor_rt_query_seconds",
            "obs_rt_seconds",
            "obs_mobility",
            "ms2_cosine_ref_similarity",
            "ms2_coelution_score",
            "ms2_summed_transition_intensity",
            "ms1_cosine_ref_similarity",
            "ms1_coelution_score",
            "nqueries",
            "delta_theo_rt",
            "sq_delta_theo_rt",
            "ms1_inten_ratio_2",
            "ms2_inten_ratio_4",
            "ms2_inten_ratio_6",
            "ms1_inten_ratio_1",
            "ms2_inten_ratio_2",
            "ms2_inten_ratio_1",
            "ms2_inten_ratio_3",
            "ms1_inten_ratio_0",
            "ms2_inten_ratio_5",
            "ms2_inten_ratio_0",
        )
        + loggable_cols
        + imputable_cols
        + zero_imputable_cols
    )
    df_use = log_cols(df, loggable_cols)
    df_use = lazy_abs_and_maxfill(df_use, imputable_cols)
    df_use = lazy_zero_fill(df_use, zero_imputable_cols)
    df_use = cast_f32(df_use, feat_cols)
    pprint("Collecting")
    df_use = add_id(df_use).collect(streaming=True)
    check_noninf(df_use, feat_cols)
    check_nonnan(df_use, feat_cols)
    pprint("Stripping sequences")
    stripped_seqs = (
        df_use["sequence"].str.replace_all("\/\d+", "").str.replace_all("\[.*?\]", "")
    )
    df_use = df_use.with_columns(
        td_id=pl.Series(
            derive_td_pair(
                stripped_seqs.to_list(), charges=df_use["precursor_charge"].to_list()
            )
        )
    )
    pprint("Initial T/D competition")
    init_nrow = len(df_use)
    df_use = df_use.filter(
        pl.col("main_score")
        == pl.col("main_score").max().over(["td_id", "precursor_charge"])
    )
    pprint(
        f"T/D competition Removed {init_nrow - len(df_use)} rows (kept {len(df_use)})"
    )
    pprint("Converting to pandas")
    df_use = df_use.to_pandas()
    cols = ColumnGroups(
        columns=df_use.columns,
        target_column="is_target",
        peptide_column="sequence",
        # spectrum_columns=("id", "td_id"),
        spectrum_columns=("td_id",),
        feature_columns=feat_cols,
        extra_confidence_level_columns=(),
        optional_columns=OptionalColumns(
            id=None,
            filename=None,
            scan=None,
            calcmass="precursor_mz",
            # charge="precursor_charge",
            charge=None,
            expmass=None,
            rt="obs_rt_seconds",
            protein=None,
        ),
    )
    nonfeat_cols = set(cols.columns) - set(cols.feature_columns)
    pprint(f"Non-feature columns: {nonfeat_cols}")
    return df_use, cols


def derive_td_pair(sequences: list[str], charges: list[int]) -> list[int]:
    pprint("Deriving TD pairs")
    td_pairs = {}
    max_id = 0

    out = []
    for seq, charge in tqdm(zip(sequences, charges, strict=True), total=len(sequences)):
        if (seq, charge) in td_pairs:
            out.append(td_pairs[(seq, charge)])
            continue

        dec = seq[0] + seq[1:-1][::-1] + seq[-1]
        td_pairs[(seq, charge)] = td_pairs[(dec, charge)] = max_id
        out.append(max_id)
        max_id += 1

    pprint(f"Found {max_id} TD pairs in {len(sequences)} sequences/charge pairs")
    return out


def to_mokapot(df: pl.LazyFrame) -> mokapot.LinearPsmDataset:
    df_use, cols = to_mokapot_df(df)
    return mokapot.LinearPsmDataset(
        df_use,
        column_groups=cols,
    )


# def to_xgb(df: pl.LazyFrame) -> xgb.DMatrix:
#     df_use, cols = to_mokapot_df(df)
#     X = df_use.loc[:, cols.feature_columns].to_numpy()
#     y = df_use.loc[:, cols.target_column].to_numpy()
#     return xgb.DMatrix(X, y, feature_names=cols.feature_columns)


def to_folds_xgb(
    *,
    shuffled_df: pd.DataFrame,
    cols: ColumnGroups,
    num_folds: int,
) -> list[xgb.DMatrix]:
    out = np.array_split(shuffled_df, num_folds)
    out = [
        xgb.DMatrix(
            x.loc[:, cols.feature_columns].to_numpy(),
            np.where(x.loc[:, cols.target_column].to_numpy(), 1, 0),
            feature_names=cols.feature_columns,
        )
        for x in out
    ]
    return out


def plot_importances(importances: dict[str, list[float]]):
    # Lollipop plot showing the importance of each feature
    fig, ax = plt.subplots(figsize=(6, 8))
    rev_imps = [(k, v) for k, v in importances.items()]
    rev_imps = rev_imps[::-1]
    for feature, importance in rev_imps:
        expanded_feat = [feature] * len(importance)
        # Stem of the lollipop
        ax.hlines(expanded_feat, 0, importance)

        # Dot at the top of the lollipop
        ax.scatter(importance, expanded_feat, c="k", marker="o")

    # Rotate the labels
    ax.tick_params(axis="x", rotation=45)

    ax.set_ylabel("Importance ('gain' as defined by xgboost)")
    ax.set_xlabel("Feature")
    # square root scale the x axis
    ax.set_xscale("log")
    fig.tight_layout()
    return fig


def xgboost_stuff(df: pl.LazyFrame, output_dir: Path):
    df_use, cols = to_mokapot_df(df)
    pprint("Shuffling")
    df_use = df_use.sample(frac=1).reset_index(drop=True, inplace=False)
    folds = to_folds_xgb(shuffled_df=df_use, num_folds=5, cols=cols)
    fold_model = KFoldModel.from_folds(folds)
    fold_model.train()
    fold_model.score()
    importances = fold_model.get_importances()
    pprint(importances)
    fig = plot_importances(importances)
    outfile = output_dir / "importances.png"
    fig.savefig(outfile)
    pprint(f"Wrote {outfile}")
    plt.close()

    ctargs = fold_model.concat_targets()
    cscores = fold_model.concat_scores()
    df_use["rescore_score"] = cscores
    df_use["qvalue"] = mokapot.qvalues.qvalues_from_scores(cscores, ctargs == 1)
    df_use = df_use.sort_values("rescore_score", ascending=False)
    outfile = output_dir / "rescored_values.parquet"
    df_use.to_parquet(outfile, index=False)
    pprint(f"Wrote {outfile}")
    order = np.argsort(-cscores)
    ctargs = ctargs[order]
    cscores = cscores[order]
    cumtargs = np.cumsum(ctargs)

    target_preds = cscores[ctargs == 1]
    decoy_preds = cscores[ctargs == 0]
    qvals = mokapot.qvalues.qvalues_from_scores(cscores, ctargs == 1)
    for ct in [0.01, 0.05, 0.1, 0.5, 1.0]:
        num_at_thresh = int(np.sum(ctargs[qvals < ct]))
        ssc = cscores[qvals < ct]
        if len(ssc) == 0:
            pprint(f"No scores at {ct}")
            continue
        score_at_thresh = np.min(ssc)
        pprint(f"Score at {ct}: {score_at_thresh}")
        pprint(f"Number of targets at {ct}: {num_at_thresh}")

    report = Report(
        targets_at_1=np.sum(ctargs[qvals < 0.01]).item(),
        targets_at_5=np.sum(ctargs[qvals < 0.05]).item(),
        targets_at_10=np.sum(ctargs[qvals < 0.1]).item(),
    )
    pprint(report)
    outfile = output_dir / "report.toml"
    report.save_to_toml(outfile)
    pprint(f"Wrote {outfile}")

    # plt.plot(qvals, cumtargs, label="All Scores")

    mask = qvals < 0.1
    plt.plot(qvals[mask], cumtargs[mask], label="QValues < 0.1")
    plt.legend(loc="upper right")
    plt.xlim(0, 0.1)
    plt.title("Cumulative number of accepted peptides.")
    outfile = output_dir / "plot_qvalues.png"
    plt.savefig(outfile)
    pprint(f"Wrote {outfile}")
    plt.close()

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].hist(target_preds, alpha=0.5, bins=100, label="Targets")
    ax[0].hist(decoy_preds, alpha=0.5, bins=100, label="Decoys")
    # ax[0].axvline(x=score_at_onepct, color="k", linestyle="--", alpha=0.5)
    ax[0].legend()

    ax[1].hist(target_preds, alpha=0.5, bins=100, label="Targets")
    ax[1].hist(decoy_preds, alpha=0.5, bins=100, label="Decoys")
    # ax[1].axvline(x=score_at_onepct, color="k", linestyle="--", alpha=0.5)
    ax[1].set_yscale("log")
    ax[1].legend()
    plt.title("Histogram of 'rescoring' score.")
    outfile = output_dir / "plot_hist.png"
    plt.savefig(outfile)
    pprint(f"Wrote {outfile}")
    plt.close()


def main_score_hist(df: pl.LazyFrame, output_dir: Path):
    pprint("Plotting main scores")
    scores_df = df.select(["is_target", "main_score"]).collect()
    target_scores = scores_df.filter(pl.col("is_target") == "true")["main_score"]
    decoy_scores = scores_df.filter(pl.col("is_target") == "false")["main_score"]

    pprint("Number of targets: {}".format(len(target_scores)))
    pprint("Number of decoys: {}".format(len(decoy_scores)))
    if (len(target_scores) + len(decoy_scores)) != len(scores_df):
        raise ValueError("Error filtering targets and decoys")

    scores = np.log1p(np.concatenate((target_scores, decoy_scores)))
    bins = np.histogram_bin_edges(scores, bins=50)

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10, 8), sharex=True)

    ax[0].hist(scores, bins, alpha=0.5, label="All Scores")
    ax[0].hist(
        np.log1p(target_scores.to_numpy()), bins, alpha=0.5, label="Target Scores"
    )
    ax[0].hist(np.log1p(decoy_scores.to_numpy()), bins, alpha=0.5, label="Decoy Scores")
    ax[0].set_xlabel("Main Score (log1p)")
    ax[0].set_ylabel("Count")
    ax[0].legend(loc="upper right")
    ax[0].set_yscale("log")

    ax[1].hist(scores, bins, alpha=0.5, label="All Scores")
    ax[1].hist(
        np.log1p(target_scores.to_numpy()), bins, alpha=0.5, label="Target Scores"
    )
    ax[1].hist(np.log1p(decoy_scores.to_numpy()), bins, alpha=0.5, label="Decoy Scores")
    ax[1].set_xlabel("Main Score (log1p)")
    ax[1].set_ylabel("Count")
    ax[1].legend(loc="upper right")

    target_file = output_dir / "plot_mainscores.png"
    pprint(f"Saving plot to {target_file}")
    plt.title("Histogram of 'main_score' scores.")
    plt.savefig(target_file)
    plt.close()


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

    ds = to_mokapot(data)
    pprint("Brewing Mokapot")
    models, scores = mokapot.brew([ds])
    qvals = mokapot.qvalues.qvalues_from_scores(scores[0], ds.targets)
    for ct in [0.01, 0.05, 0.1, 0.5, 1.0]:
        num_at_thresh = np.sum(ds.targets[qvals < ct])
        ssc = scores[0][qvals < ct]
        if len(ssc) == 0:
            pprint(f"No scores at {ct}")
            continue
        score_at_thresh = np.min(ssc)
        pprint(
            f"Mokapot Number of targets at {ct}: {num_at_thresh}; Score: {score_at_thresh}"
        )


@dataclass
class Report:
    targets_at_1: int
    targets_at_5: int
    targets_at_10: int

    def save_to_toml(self, path: Path):
        with open(path, "w") as f:
            f.write("[report]\n")
            f.write(f"targets_at_1 = {self.targets_at_1}\n")
            f.write(f"targets_at_5 = {self.targets_at_5}\n")
            f.write(f"targets_at_10 = {self.targets_at_10}\n")


@dataclass
class KFoldModel:
    """
    The idea here is to have N folds that will be used in this way:
    1. Train on 1 fold.
    2. Use 1 fold for early stopping.
    3. For prediction, use all other folds (models where the fold was not used for training or early stopping).
    """

    folds: list[xgb.DMatrix]
    models: list[xgb.Booster | None]
    scores: list[np.ndarray | None]

    @staticmethod
    def from_folds(folds: list[xgb.DMatrix]):
        return KFoldModel(folds, [None] * len(folds), [None] * len(folds))

    def train(self):
        for i in range(len(self.folds)):
            pprint(f"Training model {i}/{len(self.folds)}")
            train = self.folds[i]
            early_stop = self.folds[(i + 1) % len(self.folds)]
            model = xgb.train(
                {"objective": "binary:logistic", "scale_pos_weight": 0.2},
                train,
                num_boost_round=200,
                evals=[
                    (early_stop, "validation"),
                ],
                early_stopping_rounds=5,
                verbose_eval=1,
            )
            self.models[i] = model

    def score(self):
        for i in range(len(self.folds)):
            pprint(f"Scoring fold {i}/{len(self.folds)}")
            scores = []
            for j, model in enumerate(self.models):
                if j == i:
                    continue
                if j == (i + 1) % len(self.folds):
                    continue
                scores.append(model.predict(self.folds[i]))
            self.scores[i] = np.mean(scores, axis=0)

    def get_importances(self):
        imps = [model.get_score(importance_type="gain") for model in self.models]
        imps_order = sorted(imps[0].items(), key=lambda x: x[1], reverse=True)
        out = {k[0]: [w.get(k[0], 0) for w in imps] for k in imps_order}
        return out

    def concat_scores(self):
        if self.scores[0] is None:
            raise ValueError("Scores not computed")
        return np.concatenate(self.scores)

    def concat_targets(self):
        targets = [x.get_label() for x in self.folds]
        return np.concatenate(targets)


if __name__ == "__main__":
    parser = build_parser()
    args, unkargs = parser.parse_known_args()
    if unkargs:
        raise ValueError(f"Unknown arguments: {unkargs}")

    main(args)
