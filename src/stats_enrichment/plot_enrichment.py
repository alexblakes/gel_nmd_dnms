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
    "data/statistics/dnms_enrichment_all_genes.tsv",
    "data/statistics/dnms_enrichment_ad.tsv",
    "data/statistics/dnms_enrichment_ar.tsv",
    "data/statistics/dnms_enrichment_non_morbid.tsv",
]
_TITLES = ["All genes", "OMIM morbid (AD)", "OMIM morbid (AR)", "Non-morbid"]


logger = logging.getLogger(__name__)


def read_data(path):
    return pd.read_csv(path, sep="\t", index_col="csq")


def reverse_data(df):
    return df.iloc[::-1]


def normalise_error(df):
    return df.assign(
        fc_ci_lo=lambda x: abs(x.fc_ci_lo - x.fc), fc_ci_hi=lambda x: x.fc_ci_hi - x.fc, fc_ci_null=0
    )


def add_labels(df):
    sigfig = lambda x: x.round(1).astype(str)

    fc = lambda x: sigfig(x.fc)
    ci_lo = lambda x: sigfig(x.fc - x.fc_ci_lo)
    ci_hi = lambda x: sigfig(x.fc + x.fc_ci_hi)

    return df.assign(labels= lambda x: fc(x) + " (" + ci_lo(x) + "-" + ci_hi(x) + ")")


def log_transform(df):
    return df.transform({"fc": np.log2, "fc_ci_lo": np.log2, "fc_ci_hi": np.log2})


def customise_axes(ax=None, title=None):
    if not ax:
        plt.gca()

    ax.axvline(1, linestyle="--", color="grey", alpha=0.5)
    ax.set_title(title)
    ax.set_xlabel("Fold enrichment")
    ax.label_outer()
    ax.xaxis.set_major_locator(ticker.MaxNLocator(5))

    # If using a log scale on x
    # ax.set_xlim(left=0.5)
    # ax.set_xscale("log", base=2)
    # ax.xaxis.set_major_locator(ticker.SymmetricalLogLocator(base=2, linthresh=0.1))
    # ax.xaxis.set_major_formatter(ticker.LogFormatterExponent(base=2))

    return ax


def add_bar_labels(df, ax=None):
    if not ax:
        ax = plt.gca()
    
    labels = df.labels
    bars = ax.containers[1]
    ax.bar_label(bars, labels, padding=7, fontsize=5.5)

    return ax


def add_significance_label(ax, df):
    if not ax:
        plt.gca()

    alpha = 0.05
    n_tests = 28

    xs = df.fc + df.fc_ci_null
    ys = np.arange(len(df))
    ps = df.p < (alpha / n_tests)

    for x, y, p in zip(xs, ys, ps):
        if p:
            ax.text(x, y, "$\\star$", ha="left", va="center", fontsize=10)

    return ax


def same_xlims(axs):
    x_lims = [ax.get_xlim() for ax in axs]
    x_max = max([x for y in x_lims for x in y])
    for ax in axs:
        ax.set_xlim(right=x_max)

    return axs


def main():
    """Run as script."""

    # Instantiate figure
    plt.style.use([C.STYLE_DEFAULT, C.COLOR_ENRICHMENT])
    fig, axs = plt.subplots(1, 4, figsize=(18 * C.CM, 4 * C.CM), layout="constrained")

    # Create plots
    for path, title, ax in zip(_PATHS, _TITLES, axs):
        df = read_data(path).pipe(reverse_data).pipe(normalise_error).pipe(add_labels)
        vis.horizontal_bars(df.fc, ax, xerr=df[["fc_ci_lo", "fc_ci_null"]].T)
        customise_axes(ax, title)
        add_significance_label(ax, df)
        add_bar_labels(df, ax)

    same_xlims(axs)

    # Save figure
    plt.savefig(_PNG, dpi=600)
    plt.savefig(_SVG)
    plt.close("all")

    return fig


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
