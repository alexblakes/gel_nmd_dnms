"""Plot DNV enrichment in different regions and gene sets."""

import logging

import matplotlib.pyplot as plt
from matplotlib import ticker
import pandas as pd
import numpy as np

import src
from src import constants as C
from src import visualisation as vis

_PATHS = [
    "data/statistics/dnms_enrichment_all_genes_tidy.tsv",
    "data/statistics/dnms_enrichment_ad_tidy.tsv",
    "data/statistics/dnms_enrichment_ar_tidy.tsv",
    "data/statistics/dnms_enrichment_non_morbid_tidy.tsv",
]
_TITLES = ["All genes", "OMIM morbid (AD)", "OMIM morbid (AR)", "Non-morbid"]
_PNG = "data/plots/dnm_enrichment.png"
_SVG = "data/plots/dnm_enrichment.svg"


logger = logging.getLogger(__name__)


def read_data(path):
    return pd.read_csv(path, sep="\t", index_col="csq")


def customise_axes(ax=None, title=None):
    if not ax:
        plt.gca()

    ax.axhline(1, linestyle="--", color="grey", alpha=0.5)
    ax.set_title(title, pad=7)
    ax.set_ylabel("Fold enrichment")
    ax.label_outer()
    ax.tick_params(axis="y", labelleft=True)
    ax.set_xticks(
        ticks=ax.get_xticks(),
        labels=ax.get_xticklabels(),
        rotation=45,
        rotation_mode="anchor",
        ha="right",
    )

    # If using a log scale
    ax.set_ylim(0.5)
    ax.set_yscale("log", base=2)
    ax.yaxis.set_major_locator(ticker.SymmetricalLogLocator(base=2, linthresh=0.1))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))

    return ax


def plot(path, title, ax=None):
    ax = ax or plt.gca()

    df = read_data(path)
    vis.vertical_bars(df.fc, ax, yerr=df[["fc_ci_lo", "fc_ci_hi"]].T)
    vis.add_significance_asterisk(
        xs=np.arange(len(df)),
        ys=np.array(df["fc"] + df["fc_ci_hi"]),
        ax=ax,
        ps=df["p"] < (0.05 / 28),
        y_adj=3,
    )
    customise_axes(ax, title)

    return ax


def plots(axs, paths=_PATHS, titles=_TITLES):
    for ax, path, title in zip(axs, paths, titles):
        plot(path, title, ax)

    return axs


def main():
    """Run as script."""

    plt.style.use([C.STYLE_DEFAULT, C.COLOR_ENRICHMENT])
    fig, axs = plt.subplots(1, 4, figsize=(18 * C.CM, 5 * C.CM), layout="constrained")

    plots(axs)
    vis.same_lims(axs, y=True)

    plt.savefig(_PNG, dpi=600)
    plt.savefig(_SVG)
    plt.close("all")

    return fig


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
