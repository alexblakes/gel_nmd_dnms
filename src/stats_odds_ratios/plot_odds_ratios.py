"""Plot odds ratios."""

import logging
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
import pandas as pd
import seaborn as sns

import src
from src import constants as C
from src import visualisation as vis

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FILE_IN = "data/statistics/case_solved_odds_ratios.tsv"
_SVG = "data/plots/case_solved_odds_ratios.svg"
_PNG = "data/plots/case_solved_odds_ratios.png"

logger = logging.getLogger(__name__)


def read_data(path=_FILE_IN):
    return pd.read_csv(path, sep="\t")


def h_point_range(x, y, ax=None, scatter_kwargs={}, errorbar_kwargs={}):
    ax = ax or plt.gca()

    ax.scatter(x=x, y=y, **scatter_kwargs)
    ax.errorbar(x=x, y=y, **errorbar_kwargs)

    return ax


def customise_axes(ax=None, ylabel=None, **set_kwargs):
    ax = ax or plt.gca()

    ax.set_ylim(-1, 2)
    ax.set_yticks(ticks=[np.mean(ax.get_ylim())], labels=ylabel)
    ax.tick_params(axis="both", length=0)
    ax.tick_params(axis="y", pad=5)

    ax.xaxis.set_major_locator(ticker.MaxNLocator(min_n_ticks=10))

    ax.axvline(x=1, c="white", linestyle="--", zorder=-1)

    ax.set_xlabel("Odds ratio\n(case solved)")
    ax.label_outer()

    for spine in ax.spines.values():
        spine.set_edgecolor("white")
    ax.spines["left"].set_visible(False)

    ax.set(**set_kwargs)

    return ax


def plot_odds_ratios(
    df, ax, facecolor="white", label=None, xticks=False, legend=False, **kwargs
):

    if not xticks:
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(axis="x", bottom=False)
        ax.set_xlabel(None)

    pass


def add_legend(ax=None, line2d_kwargs={}, legend_kwargs={}):
    line2d_kwargs.setdefault("color", "black")
    line2d_kwargs.setdefault("marker", "o")
    line2d_kwargs.setdefault("markersize", 4)

    legend_kwargs.setdefault("markerfirst", False)
    legend_kwargs.setdefault("loc", "center right")

    ax = ax or plt.gca()

    legend_elements = [plt.Line2D([0], [0], **line2d_kwargs)]

    ax.legend(handles=legend_elements, **legend_kwargs)

    return ax


def plot(axs):

    df = read_data()

    data_subsets = [g for _, g in df.groupby("label", sort=False)]
    bg_colors = [vis.adjust_alpha(c, 0.7) for c in sns.color_palette()]

    for data, ax, color in zip(data_subsets, axs, bg_colors):
        x = data["odds_ratio"]
        y = np.arange(len(data))
        ylabel = data["label"].unique()
        xerr = data[["err_lo", "err_hi"]].T
        ps = data["bfr_sig"]
        geom_colors = ["black", "dimgrey"]

        h_point_range(
            x,
            y,
            ax=ax,
            scatter_kwargs=dict(color=geom_colors),
            errorbar_kwargs=dict(xerr=xerr, linestyle="", ecolor=geom_colors),
        )
        customise_axes(ax, ylabel, facecolor=color)
        vis.add_significance_asterisk(x, y, ps, ax, color="white", s=3)

    add_legend(axs[0], line2d_kwargs=dict(label="Unconstrained", color="dimgrey"))
    add_legend(axs[1], line2d_kwargs=dict(label="Constrained", color="black"))

    return axs


def main():
    """Run as script."""

    plt.style.use([C.STYLE_DEFAULT, C.COLOR_ENRICHMENT])
    fig, axs = plt.subplots(7, 1, figsize=(7 * C.CM, 5 * C.CM), sharex=True)
    plt.subplots_adjust(hspace=0)

    plot(axs)

    plt.savefig(_SVG)
    plt.savefig(_PNG, dpi=600)
    plt.close("all")

    return axs


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
