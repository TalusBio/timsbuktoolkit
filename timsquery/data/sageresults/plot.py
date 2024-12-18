# /// script
# dependencies = [
#   "polars",
#   "matplotlib",
# ]
# ///

import json
import numpy as np
import matplotlib.pyplot as plt
import argparse


def ms_to_mins(x):
    return (x / 1_000) / 60


def infinite_color_cycle():
    colors = ["#308AAD", "#C8102E", "#96D8D8", "#B2EE79"]
    while True:
        for c in colors:
            yield c


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--elution-groups", type=str, default="ubb_elution_groups.json")
    parser.add_argument(
        "--query-results", type=str, default="query_results/results.json"
    )
    parser.add_argument("--output", type=str, default="ubb_peptide_plot.png")
    return parser


def main(args):
    data = json.load(open(args.query_results))
    labels = json.load(open(args.elution_groups))
    id_to_label = {x["id"]: x["human_id"] for x in labels}

    num_elution_groups = len(data)
    fig, ax = plt.subplots(ncols=num_elution_groups, nrows=2, figsize=(10, 7))

    for i, x in enumerate(data):
        ms1_rts = np.array(x["result"]["ms1_stats"]["retention_time_miliseconds"])
        ms1_rts = ms_to_mins(ms1_rts)
        ms2_rts = np.array(x["result"]["ms2_stats"]["retention_time_miliseconds"])
        ms2_rts = ms_to_mins(ms2_rts)

        colors = infinite_color_cycle()
        sorted_ms1_keys = sorted(
            x["result"]["ms1_stats"]["transition_intensities"].keys()
        )
        sorted_ms2_keys = sorted(
            x["result"]["ms2_stats"]["transition_intensities"].keys()
        )

        rt_plot_center = x["elution_group"]["rt_seconds"] / 60
        rt_plot_right = rt_plot_center + 0.5
        rt_plot_left = rt_plot_center - 0.5
        max_ms1 = 0
        max_ms2 = 0

        for j in sorted_ms1_keys:
            col = next(colors)
            j = x["result"]["ms1_stats"]["transition_intensities"][j]
            max_ms1 = max(max_ms1, max(j))
            ax[0, i].plot(ms1_rts, j, color=col, alpha=0.4)
            # Draw the points
            ax[0, i].scatter(ms1_rts, j, color=col, s=5)

        for j in sorted_ms2_keys:
            col = next(colors)
            j = x["result"]["ms2_stats"]["transition_intensities"][j]
            max_ms2 = max(max_ms2, max(j))
            ax[1, i].plot(ms2_rts, j, color=col, alpha=0.4)
            # Draw the points
            ax[1, i].scatter(ms2_rts, j, color=col, s=5)

        ax[0, i].set_xlim(rt_plot_left, rt_plot_right)
        ax[1, i].set_xlim(rt_plot_left, rt_plot_right)

        ax[0, i].vlines(
            rt_plot_center,
            0,
            max_ms1,
            colors="red",
            linestyles="dashed",
            alpha=0.3,
        )
        ax[1, i].vlines(
            rt_plot_center,
            0,
            max_ms2,
            colors="red",
            linestyles="dashed",
            alpha=0.3,
        )

        # Label axes ...
        ax[0, i].set_xlabel("Retention Time (min)")
        ax[1, i].set_xlabel("Retention Time (min)")
        ax[0, i].set_ylabel("Intensity")
        ax[1, i].set_ylabel("Intensity")

        # Label titles
        human_lab = id_to_label[x["elution_group"]["id"]]
        ax[0, i].set_title(human_lab + " MS1")
        ax[1, i].set_title(human_lab + " MS2")

    plt.tight_layout()
    plt.savefig(args.output)


if __name__ == "__main__":
    parser = build_parser()
    args = parser.parse_args()
    main(args)
