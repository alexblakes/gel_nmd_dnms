"""Plot the number of transcripts with recurrent dnPTVs per region."""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import src
from src import constants as C
from src import visualisation as vis

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FILE_IN = "data/statistics/recurrent_dnms.tsv"
_PNG = "data/plots/recurrent_dnms.png"
_SVG = "data/plots/recurrent_dnms.svg"

logger = logging.getLogger(__name__)


def read_data(path):
    return pd.read_csv(path, sep="\t", index_col="n")


def customise_axes(ax=None, **set_kwargs):
    if not ax:
        ax = plt.gca()

    ax.set(**set_kwargs)

    ax.set_yscale("log")
    ax.set_xticks(ticks=ax.get_xticks(), labels=["1", "2", "3", "4", "5+"])

    # Bar labels
    bars = ax.containers[0]
    ax.bar_label(bars)

    return ax


def same_ylims(axs):
    y_lims = [ax.get_ylim() for ax in axs]
    y_max = max([x for y in y_lims for x in y])
    for ax in axs:
        ax.set_ylim(top=y_max)

    return axs


def main():
    """Run as script."""

    df = read_data(_FILE_IN)

    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_REGIONS)

    fig, axs = plt.subplots(2, 3, figsize=(12 * C.CM, 8 * C.CM), layout="constrained")
    axs = axs.flatten()

    colors = [sns.color_palette()[x] for x in [0, 0, 0, 1, 3, 3]]

    ylabs = ["Regions", None, None, "Transcripts", None, None]

    titles = [
        "Any region",
        "Constrained region in any gene",
        "Constrained region\nin a non-morbid gene",
        "Constrained NMD-target region\nin a non-morbid gene",
        "Across any constrained\nNMD-escape region\nin a non-morbid gene",
        "Within one constrained\nNMD-escape region\nin a non-morbid gene",
    ]

    for ax, column, ylab, title, color in zip(axs, df.columns, ylabs, titles, colors):
        vis.vertical_bars(df.loc[:, column], ax, color=color)
        customise_axes(ax, ylabel=ylab, xlabel="$\it{dn}$PTVs", title=title)

    same_ylims(axs)

    plt.savefig(_PNG, dpi=600)
    plt.savefig(_SVG)
    plt.close("all")

    return df


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
