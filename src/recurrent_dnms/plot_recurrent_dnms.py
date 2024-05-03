"""Module docstring."""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import src
from src import constants as C
from src.visualisation import plots

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FIG_OUT = "data/plots/recurrent_dnms"

logger = logging.getLogger(__name__)


def read_dnm_ptvs(path):
    return pd.read_csv(
        path,
        sep="\t",
        usecols=[
            "enst",
            "chr",
            "pos",
            "ref",
            "alt",
            "symbol",
            "omim_inheritance_simple",
            "region",
            "constraint",
            "loeuf",
            "cohort",
            "id",
            "csq",
        ],
    ).query("csq == 'frameshift_variant' | csq == 'stop_gained'")


def count_dnms_per_region(df, name):
    counts = (
        df.groupby(["enst", "region"])["chr"]
        .count()
        .groupby("enst")
        .max()
        .value_counts()
        .reset_index()
        .set_axis(["n", "count"], axis=1)
    )
    counts["n"] = counts["n"].clip(upper=5)
    counts = counts.groupby("n")["count"].sum().rename(name)

    return counts


def combine_dnm_counts(ptvs):
    m1 = ptvs.columns
    m2 = ptvs.constraint == "constrained"
    m3 = ~ptvs.omim_inheritance_simple.str.contains("AD|XL").fillna(False)
    m4 = m2 & m3

    masks = [m1, m2, m3, m4]
    names = [
        "All regions",
        "Constrained regions",
        "Not OMIM morbid (AD or XL)",
        "Constrained and not OMIM morbid (AD or XL)",
    ]
    zipped = zip(masks, names)

    combined_counts = (
        pd.concat([count_dnms_per_region(ptvs[m], n) for m, n in zipped], axis=1)
        .melt(ignore_index=False, var_name="mask", value_name="count")
        .set_index("mask", append=True)
        .swaplevel(0, 1)
    )

    return combined_counts


def plot_bars(data, ax=None, legend=False):
    if not ax:
        ax = plt.gca()

        return None


def main():
    """Run as script."""

    df = read_dnm_ptvs(C.DNMS_ANNOTATED_CLINICAL).pipe(combine_dnm_counts)

    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_VIBRANT)

    fig, ax = plt.subplots(1, 1, figsize=(12 * C.CM, 4 * C.CM), layout="constrained")

    plots.grouped_bars(df)
    plots.annotate_grouped_bars(
        legend=list(df.index.get_level_values(0).unique()),
        legend_kwargs=dict(bbox_to_anchor=(1, 1.05), borderpad=0),
        set_kwargs=dict(
            xticklabels=[1, 2, 3, 4, "5+"],
            xlabel="$\it{De\ novo}$ PTV count",
            ylabel="Number of transcripts",
            yscale="log",
        ),
    )

    plt.savefig(f"{_FIG_OUT}.png", dpi=600)
    plt.savefig(f"{_FIG_OUT}.svg")
    plt.close("all")

    return df


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
