"""Plot relative enrichment of DNMs in constrained regions."""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import src
from src import constants as C
from src import visualisation as vis

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_ALL_GENES = "data/statistics/dnms_enrichment_all_genes.tsv"
_MORBID = "data/statistics/dnms_enrichment_morbid_genes.tsv"
_NON_MORBID = "data/statistics/dnms_enrichment_non_morbid_genes.tsv"
_SVG = "data/plots/dnm_enrichment_constrained_regions.svg"
_PNG = "data/plots/dnm_enrichment_constrained_regions.png"

logger = logging.getLogger(__name__)


def read_enrichment_data(path):
    df = pd.read_csv(
        path,
        sep="\t",
        index_col="csq",
    )
    return df


def adjust_ci(df):
    df = df.assign(ci_l=lambda x: abs(x["relative_enrichment"] - x["ci_l"]))
    df = df.assign(ci_r=lambda x: abs(x["relative_enrichment"] - x["ci_r"]))
    return df


def parse_data(path):
    return read_enrichment_data(path).pipe(adjust_ci)


def horizontal_bars(
    values,
    ax=None,
    **kwargs,
):

    kwargs.setdefault("tick_label", values.index)
    kwargs.setdefault("color", PALETTE)
    kwargs.setdefault("ecolor", [vis.adjust_lightness(c, 0.8) for c in PALETTE])

    if not ax:
        ax = plt.gca()

    n = len(values)  # Number of bars
    height = 1 - (1/n)
    y = np.arange(n)

    ax.barh(y=y, width=values, height=height, **kwargs)

    ax.axvline(x=1, linestyle="--", color="grey", alpha=0.5)

    return None


def main():
    """Run as script."""

    # Load data
    all_genes = parse_data(_ALL_GENES)
    morbid = parse_data(_MORBID)
    non_morbid = parse_data(_NON_MORBID)

    # # Subset the morbid genes by mode of inheritance
    morbid_ad = morbid[morbid["inheritance_simple"] == "AD"]
    morbid_ar = morbid[morbid["inheritance_simple"] == "AR"]

    fig, axs = plt.subplots(
        1, 4, figsize=(18 * C.CM, 4 * C.CM), layout="constrained", sharey=True,
    )

    horizontal_bars(
        all_genes["relative_enrichment"],
        axs[0],
        xerr=[all_genes["ci_l"], all_genes["ci_r"]],
    )

    horizontal_bars(
        morbid_ad["relative_enrichment"],
        axs[1],
        xerr=[morbid_ad["ci_l"], morbid_ad["ci_r"]],
    )

    horizontal_bars(
        morbid_ar["relative_enrichment"],
        axs[2],
        xerr=[morbid_ar["ci_l"], morbid_ar["ci_r"]],
    )

    horizontal_bars(
        non_morbid["relative_enrichment"],
        axs[3],
        xerr=[non_morbid["ci_l"], non_morbid["ci_r"]],
    )

    axs[0].set_title("All genes")
    axs[1].set_title("Morbid genes\n(dominant)")
    axs[2].set_title("Morbid genes\n(recessive)")
    axs[3].set_title("Non-morbid genes")
  
    for ax in axs:
        ax.set_xlim(0,11)
        ax.set_xlabel("Fold enrichment")
        
    plt.savefig(_SVG)
    plt.savefig(_PNG, dpi=1000)
    
    pass


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)

    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_ENRICHMENT)
    PALETTE = sns.color_palette()

    main()