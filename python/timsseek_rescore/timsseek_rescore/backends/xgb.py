from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import mokapot
import numpy as np
import polars as pl
import xgboost as xgb
from rich.pretty import pprint
from uniplot import histogram

from ..datamodels import Report
from ..feateng import to_mokapot_df
from ..folding import to_folds_xgb
from ..plotting import plot_importances


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
            model = self._train_fold(i)
            self.models[i] = model

    def _train_fold(self, fold_index: int) -> xgb.Booster:
        train = self.folds[fold_index]
        early_stop = self.folds[(fold_index + 1) % len(self.folds)]

        model = self._train_model(train, early_stop)
        return model

    def _train_model(self, train: xgb.DMatrix, early_stop: xgb.DMatrix) -> xgb.Booster:
        model = xgb.train(
            {"objective": "binary:logistic", "scale_pos_weight": 0.5},
            # {"objective": "binary:logistic"},
            train,
            num_boost_round=200,
            evals=[
                (early_stop, "validation"),
            ],
            early_stopping_rounds=5,
            verbose_eval=1,
        )
        return model

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


def tdc(scores: np.ndarray[float], target: np.ndarray[bool], desc: bool = True):
    """Estimate q-values using target decoy competition.

    Estimates q-values using the simple target decoy competition method.
    For set of target and decoy PSMs meeting a specified score threshold,
    the false discovery rate (FDR) is estimated as:

    ...math:
        FDR = \frac{Decoys + 1}{Targets}

    More formally, let the scores of target and decoy PSMs be indicated as
    :math:`f_1, f_2, ..., f_{m_f}` and :math:`d_1, d_2, ..., d_{m_d}`,
    respectively. For a score threshold :math:`t`, the false discovery
    rate is estimated as:

    ...math:
        E\\{FDR(t)\\} = \frac{|\\{d_i > t; i=1, ..., m_d\\}| + 1}
        {\\{|f_i > t; i=1, ..., m_f|\\}}

    The reported q-value for each PSM is the minimum FDR at which that
    PSM would be accepted.

    With one exception, the lowest score will always have a q-value of 1.0.

    Parameters
    ----------
    scores : numpy.ndarray of float
        A 1D array containing the score to rank by
    target : numpy.ndarray of bool
        A 1D array indicating if the entry is from a target or decoy
        hit. This should be boolean, where `True` indicates a target
        and `False` indicates a decoy. `target[i]` is the label for
        `metric[i]`; thus `target` and `metric` should be of
        equal length.
    desc : bool
        Are higher scores better? `True` indicates that they are,
        `False` indicates that they are not.

    Returns
    -------
    numpy.ndarray
        A 1D array with the estimated q-value for each entry. The
        array is the same length as the `scores` and `target` arrays.
    """
    # Since numpy 2.x relying in attribute errors is not viable here
    # https://numpy.org/neps/nep-0050-scalar-promotion.html#impact-on-can-cast
    # So I am manually checking the constraints.
    if (
        np.issubdtype(target.dtype, np.integer)
        and target.max() <= 1
        and target.min() >= 0
    ):
        target = target.astype(bool)

    if np.issubdtype(target.dtype, np.floating):
        like_one = target == np.ones_like(target)
        like_zero = target == np.zeros_like(target)
        if np.all(like_one | like_zero):
            target = target.astype(bool)

    if not np.issubdtype(target.dtype, bool):
        err = ValueError(
            f"'target' should be boolean. passed type: {target.dtype}"
            f" with value: {target}"
        )
        raise err

    if scores.shape[0] != target.shape[0]:
        raise ValueError("'scores' and 'target' must be the same length")

    # Unsigned integers can cause weird things to happen.
    # Convert all scores to floats to for safety.
    scores = scores.astype(np.float32)

    # Sort and estimate FDR
    # Sort order is first by score and then ties are
    # sorted by target
    if desc:
        srt_idx = np.lexsort((target, -scores))
    else:
        srt_idx = np.lexsort((target, scores))

    scores = scores[srt_idx]
    target = target[srt_idx]

    cum_targets = target.cumsum()
    cum_decoys = (~target).cumsum()

    # Handles zeros in denominator
    fdr = np.divide(
        (cum_decoys + 1),
        cum_targets,
        out=np.ones_like(cum_targets, dtype=np.float32),
        where=(cum_targets != 0),
    )
    # Clamp the FDR to 1.0
    fdr = np.minimum(fdr, 1.0)

    # Sort by scores and resolve ties with fdr
    # This implementation is 1.5x slower with small data
    # from 0.05 to 0.07 miliseconds.
    # But up to 100x faster with large data
    # and does not need compilation (or numba as a dependency)
    # For the former implementation check the git history
    if desc:
        sorting = np.lexsort((fdr, scores))
    else:
        sorting = np.lexsort((fdr, -scores))
    tmp = np.minimum.accumulate(fdr[sorting])
    # Set the FDR to 1 for the lowest score
    # This prevent prevents a bug where features like the charge
    # would seem to be very good, because ... since they tie a lot
    # of PSMs, and we report the 'best' FDR for all ties, it would
    # artifially yield very low q-values.
    tmp[tmp == tmp[0]] = 1.0
    np_qval = np.flip(tmp)[np.argsort(srt_idx)]

    return np_qval


def xgboost_stuff(shuffled_df: pl.DataFrame, cols: list[str], output_dir: Path):
    folds = to_folds_xgb(shuffled_df=shuffled_df, num_folds=5, cols=cols)
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
    shuffled_df["rescore_score"] = cscores
    shuffled_df["qvalue"] = mokapot.qvalues.qvalues_from_scores(cscores, ctargs == 1)
    shuffled_df = shuffled_df.sort_values("rescore_score", ascending=False)
    outfile = output_dir / "rescored_values.parquet"
    shuffled_df.to_parquet(outfile, index=False)

    one_pct_df = shuffled_df[shuffled_df["qvalue"] < 0.01]
    if len(one_pct_df) > 0:

        def plot_hexbin(
            target_df,
            decoy_df,
            x_col_label,
            y_col_label,
            outfile,
            title,
            add_one_to_one: bool = False,
        ):
            fig, ax = plt.subplots(1, 2, figsize=(10, 5))
            x_col, xlabel = x_col_label
            y_col, ylabel = y_col_label
            for i, (df, sub_title) in enumerate(
                zip([target_df, decoy_df], ["Targets", "Decoys"])
            ):
                ax[i].hexbin(
                    df[x_col],
                    df[y_col],
                    gridsize=50,
                    cmap="viridis",
                    bins="log",
                )
                ax[i].set_title(f"{title} ({sub_title})")
                ax[i].set_xlabel(xlabel)
                ax[i].set_ylabel(ylabel)
                if add_one_to_one:
                    min_val = min(
                        df[x_col].min(),
                        df[y_col].min(),
                    )
                    max_val = max(
                        df[x_col].max(),
                        df[y_col].max(),
                    )
                    ax[i].plot([min_val, max_val], [min_val, max_val], "r--", lw=1)
            fig.tight_layout()
            fig.savefig(outfile)
            plt.close()
            pprint(f"Wrote {outfile}")

        target_df = one_pct_df[one_pct_df["is_target"]]
        decoy_df = one_pct_df[~one_pct_df["is_target"]]

        plot_hexbin(
            target_df,
            decoy_df,
            ("ms2_mz_error_0", "Mass error (m/z)"),
            ("obs_rt_seconds", "Observed RT (s)"),
            output_dir / "mass_error_rt_1pct.png",
            "1% FDR",
        )

        plot_hexbin(
            target_df,
            decoy_df,
            ("ms2_mz_error_0", "Mass error (m/z)"),
            ("precursor_mz", "Precursor m/z"),
            output_dir / "mass_error_mz_1pct.png",
            "1% FDR",
        )

        plot_hexbin(
            target_df,
            decoy_df,
            ("obs_mobility", "Observed Mobility"),
            ("obs_rt_seconds", "Observed RT (s)"),
            output_dir / "mobility_rt_1pct.png",
            "1% FDR",
        )

        # For mobility error, subtract precursor_mobility_query from obs_mobility
        one_pct_df = one_pct_df.copy()
        one_pct_df["mobility_error"] = (
            one_pct_df["obs_mobility"] - one_pct_df["precursor_mobility_query"]
        )
        target_df = one_pct_df[one_pct_df["is_target"]]
        decoy_df = one_pct_df[~one_pct_df["is_target"]]
        plot_hexbin(
            target_df,
            decoy_df,
            ("mobility_error", "Observed Mobility Error"),
            ("obs_rt_seconds", "Observed RT (s)"),
            output_dir / "mobility_error_rt_1pct.png",
            "1% FDR",
        )

        plot_hexbin(
            target_df,
            decoy_df,
            ("precursor_mz", "Precursor m/z"),
            ("obs_mobility", "Observed Mobility"),
            output_dir / "mz_mobility_1pct.png",
            "1% FDR",
        )

        plot_hexbin(
            target_df,
            decoy_df,
            ("precursor_rt_query_seconds", "Predicted RT (s)"),
            ("obs_rt_seconds", "Observed RT (s)"),
            output_dir / "predicted_rt_obs_rt_1pct.png",
            "1% FDR",
            add_one_to_one=True,
        )
    else:
        pprint("No values at 1% FDR")

    pprint(f"Wrote {outfile}")
    order = np.argsort(-cscores)
    ctargs = ctargs[order]
    cscores = cscores[order]
    cumtargs = np.cumsum(ctargs)

    target_preds = cscores[ctargs == 1]
    decoy_preds = cscores[ctargs == 0]

    histogram(target_preds, title="Target scores")
    histogram(decoy_preds, title="Decoy scores")

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
