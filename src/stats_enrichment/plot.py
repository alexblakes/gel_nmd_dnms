"""Plot DNV enrichment in different regions and gene sets."""

import logging

import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
import pandas as pd
import seaborn as sns

import src
from src import constants as C
from src import visualisation as vis

_PNG = "data/plots/dnm_enrichment.png"
_SVG = "data/plots/dnm_enrichment.svg"
_PATHS = [
    "data/statistics/dnms_enrichment_all_genes_tidy.tsv",
    "data/statistics/dnms_enrichment_ad_tidy.tsv",
    "data/statistics/dnms_enrichment_ar_tidy.tsv",
    "data/statistics/dnms_enrichment_non_morbid_tidy.tsv",
]
_TITLES = ["All genes", "OMIM morbid (AD)", "OMIM morbid (AR)", "Non-morbid"]


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
    ax.set_xticks(ticks=ax.get_xticks(), labels=ax.get_xticklabels(), rotation=45, rotation_mode="anchor", ha="right")

    # If using a log scale
    ax.set_ylim(bottom=0.5)
    ax.set_yscale("log", base=2)
    ax.yaxis.set_major_locator(ticker.SymmetricalLogLocator(base=2, linthresh=0.1))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))

    return ax


def add_significance_label(ax, df, alpha=0.05, n_tests=28):
    if not ax:
        plt.gca()

    xs = ax.get_xticks()
    ys = df.fc + df.fc_ci_hi
    ps = df.p < (alpha / n_tests)

    for x, y, p in zip(xs, ys, ps):
        if p:
            ax.text(x, y, r"$\star$", ha="center", va="bottom")

    return ax


def plot(path, title, ax=None):
    ax = ax or plt.gca()

    df = read_data(path)
    vis.vertical_bars(df.fc, ax, yerr=df[["fc_ci_lo", "fc_ci_hi"]].T)
    customise_axes(ax, title)
    add_significance_label(ax, df)

    return ax

def main():
    """Run as script."""

    plt.style.use([C.STYLE_DEFAULT, C.COLOR_ENRICHMENT])
    fig, axs = plt.subplots(1, 4, figsize=(18 * C.CM, 5 * C.CM), layout="constrained")

    for path, title, ax in zip(_PATHS, _TITLES, axs):
        plot(path, title, ax)

    vis.same_lims(axs, y=True)

    plt.savefig(_PNG, dpi=600)
    plt.savefig(_SVG)
    plt.close("all")

    return fig


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
