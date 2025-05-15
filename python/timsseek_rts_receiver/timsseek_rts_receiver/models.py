from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from pydantic import BaseModel, Field


class ArrayResponse(BaseModel):
    arr: List[List[float]]
    rts_ms: List[int]
    mz_order: Union[List[Tuple[int, float]], list[Tuple[str, float]]]

    def plot_transition(self, min_rt, max_rt):
        # In the same figure/axes show every intensity across retention time
        rt_use = np.array(self.rts_ms)
        ranges = np.searchsorted(rt_use, [min_rt, max_rt])
        rt_plot = (rt_use[ranges[0] : ranges[1]] / 1000) / 60

        fig, ax = plt.subplots()
        for k, v in zip(self.mz_order, self.arr, strict=True):
            ax.plot(rt_plot, v[ranges[0] : ranges[1]], label=k)

        ax.set_xlabel("Retention Time (min)")
        ax.set_ylabel("Intensity")
        ax.legend()
        fig.tight_layout()
        return fig


class Extractions(BaseModel):
    # eg: ElutionGroup
    fragments: ArrayResponse
    precursors: ArrayResponse


class MainScoreElements(BaseModel):
    ms1_coelution_score: List[float]
    ms1_cosine_ref_sim: List[float]
    ms2_coelution_score: List[float]
    ms2_cosine_ref_sim: List[float]
    ms2_lazyscore: List[float]
    # Should Nones be allowed??
    ms2_lazyscore_vs_baseline: List[float | None]
    ms2_corr_v_gauss: List[float]
    # TODO: REMOVE
    # hyperscore: List[float]
    split_lazyscore: List[float]
    # END

    ref_time_ms: List[int]
    ms2_lazyscore_vs_baseline_std: float

    def plot(self, min_rt_ms, max_rt_ms, vlines_ms: Optional[List[int]] = None):
        # Make a plot grid, where each row is a different score element
        # but all share the same retention time axis
        ncol = 3
        nrow = 3
        fig, ax = plt.subplots(nrows=nrow, ncols=ncol, figsize=(10, 12))

        rt_use = np.array(self.ref_time_ms)
        ranges = np.searchsorted(rt_use, [min_rt_ms, max_rt_ms])
        rt_plot = (rt_use[ranges[0] : ranges[1]] / 1000) / 60

        mask_label = np.zeros(len(self.ref_time_ms), dtype=bool)
        if vlines_ms is not None:
            for vline in vlines_ms:
                local_range = np.searchsorted(rt_use, [vline - 5_000, vline + 5_000])
                mask_label[local_range[0] : local_range[1]] = True

        score_name_pairs = [
            ("MS1 Coelution Score", self.ms1_coelution_score),
            ("MS1 Cosine Ref Sim", self.ms1_cosine_ref_sim),
            ("MS2 Coelution Score", self.ms2_coelution_score),
            ("MS2 Cosine Ref Sim", self.ms2_cosine_ref_sim),
            ("MS2 LazyScore", self.ms2_lazyscore),
            ("MS2 LazyScore Baseline", self.ms2_lazyscore_vs_baseline),
            # ("Hyperscore", self.hyperscore),
            ("Split LazyScore", self.split_lazyscore),
            ("MS2 Corr v Gauss", self.ms2_corr_v_gauss),
        ]

        for i, (name, score) in enumerate(score_name_pairs):
            local_score = np.array(score).astype(float)
            score_ir = local_score[mask_label]
            score_or = local_score[~mask_label]

            # Signal to noise ratio would be (max(ir) - max(or))/std(or)
            max_ir = np.nanmax(score_ir) if len(score_ir) > 0 else 0
            max_or = np.nanmax(score_or) if len(score_or) > 0 else 0
            std_or = np.nanstd(score_or) if len(score_or) > 0 else 1e-3
            snr = (max_ir - max_or) / std_or

            ax[i // ncol, i % ncol].plot(rt_plot, score[ranges[0] : ranges[1]])
            ax[i // ncol, i % ncol].set_title(name)

            ax[i // ncol, i % ncol].text(
                0.95,
                0.95,
                f"SNR: {np.round(snr, 2)}",
                transform=ax[i // ncol, i % ncol].transAxes,
                ha="right",
                va="top",
            )

        if vlines_ms is not None:
            vlines_minutes = np.array(vlines_ms) / 1000 / 60
            for li in range(len(vlines_minutes)):
                for r in range(nrow):
                    for c in range(ncol):
                        ax[r, c].axvline(
                            x=vlines_minutes[li], color="k", linestyle="--", alpha=0.5
                        )

        for i in range(nrow):
            for j in range(ncol):
                ax[i, j].set_xlabel("Retention Time (min)")

        fig.tight_layout()
        return fig

    def min_rt(self):
        return min(self.ref_time_ms)

    def max_rt(self):
        return max(self.ref_time_ms)


class SearchResults(BaseModel):
    is_target: bool
    lazyerscore: float
    norm_lazyerscore_vs_baseline: float
    lazyerscore_vs_baseline: float
    ms2_corr_v_gauss: float
    main_score: float
    delta_next: float
    delta_second_next: float
    nqueries: int
    ms1_coelution_score: float
    ms1_cosine_ref_similarity: float
    ms1_mobility_error_0: float
    ms1_mobility_error_1: float
    ms1_mobility_error_2: float
    ms1_mz_error_0: float
    ms1_mz_error_1: float
    ms1_mz_error_2: float
    ms1_summed_precursor_intensity: float
    ms2_coelution_score: float
    ms2_cosine_ref_similarity: float
    ms2_mobility_error_0: None | float
    ms2_mobility_error_1: None | float
    ms2_mobility_error_2: None | float
    ms2_mobility_error_3: None | float
    ms2_mobility_error_4: None | float
    ms2_mobility_error_5: None | float
    ms2_mobility_error_6: None | float
    ms2_mz_error_0: None | float
    ms2_mz_error_1: None | float
    ms2_mz_error_2: None | float
    ms2_mz_error_3: None | float
    ms2_mz_error_4: None | float
    ms2_mz_error_5: None | float
    ms2_mz_error_6: None | float

    ms1_inten_ratio_0: float
    ms1_inten_ratio_1: float
    ms1_inten_ratio_2: float
    ms2_inten_ratio_0: float
    ms2_inten_ratio_1: float
    ms2_inten_ratio_2: float
    ms2_inten_ratio_3: float
    ms2_inten_ratio_4: float
    ms2_inten_ratio_5: float
    ms2_inten_ratio_6: float

    delta_ms1_ms2_mobility: None | float
    sq_delta_ms1_ms2_mobility: None | float

    delta_theo_rt: float
    sq_delta_theo_rt: float

    ms2_summed_transition_intensity: float
    npeaks: int
    obs_mobility: float
    obs_rt_seconds: float
    precursor_charge: int
    precursor_mobility_query: float
    precursor_mz: float
    precursor_rt_query_seconds: float
    sequence: str

    def as_table(self):
        return pd.DataFrame({
            "key": self.model_dump().keys(),
            "value": [str(x) for x in self.model_dump().values()],
        })


class ResponseData(BaseModel):
    extractions: Extractions
    main_score_elements: MainScoreElements
    longitudinal_main_score: List[float]
    search_results: SearchResults

    def plot_main_score(self, min_rt, max_rt):
        fig, ax = plt.subplots()
        peptide = self.search_results.sequence
        charge = self.search_results.precursor_charge

        rts = np.array(self.main_score_elements.ref_time_ms)
        ranges = np.searchsorted(rts, [min_rt, max_rt])
        rt_plot = (rts[ranges[0] : ranges[1]] / 1000) / 60
        ax.plot(rt_plot, self.longitudinal_main_score[ranges[0] : ranges[1]])
        ax.set_xlabel("Retention Time (min)")
        ax.set_ylabel("Main Score")
        # Title
        ax.set_title(
            f"Main Score for {peptide} (z={charge})\n"
            f"Best RT: {self.search_results.obs_rt_seconds / 60:.2f} min"
        )
        fig.tight_layout()
        return fig


class Response(BaseModel):
    status: str
    data: ResponseData
