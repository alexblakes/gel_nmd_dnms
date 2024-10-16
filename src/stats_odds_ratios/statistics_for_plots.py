"""Statistics for OR plots."""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import src
from src import constants as C
from src import visualisation as vis

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FILE_IN = "data/statistics/case_solved_odds_ratios.tsv"
_SVG = "data/plots/case_solved_odds_ratios.svg"
_PNG = "data/plots/case_solved_odds_ratios.png"

logger = logging.getLogger(__name__)


def plot_odds_ratios(
    df,
    ax,
    facecolor="white",
    label=None,
    xticks=False,
    legend=False,
    **kwargs
):

    y = np.arange(len(df))[::-1]
    colors = {
        "unconstrained": "grey",
        "constrained": "black",
    }

    for i, constraint in enumerate(df["constraint"]):
        ax.errorbar(
            x=[df.iloc[i]["odds_ratio"]],
            y=[(len(df) - 1) - i],
            fmt="o",
            markersize=4,
            color=colors[constraint],
            xerr=(
                [abs(df.iloc[i]["ci_lo"] - df.iloc[i]["odds_ratio"])],
                [abs(df.iloc[i]["ci_hi"] - df.iloc[i]["odds_ratio"])],
            ),
            label=constraint.capitalize()
        )
        
    if legend:
        ax.legend(loc="center right",)

    ax.axvline(x=1, c="grey", linestyle="--")

    ax.set_xlabel("Odds ratio\n(case solved)")
    ax.set_ylim(y.min() - 1, y.max() + 1)
    ax.set_xlim(-1, 14)
    ax.set_yticks(ticks=[y.max()/2], labels=[label])
    ax.tick_params(axis="y", length=0)
    # ax.label_outer(remove_inner_ticks=True)

    ax.set_facecolor(facecolor)
    ax.patch.set_alpha(0.3)

    ax.spines["left"].set_visible(False)
    ax.spines["left"].set_position(("outward", 7))

    if not xticks:
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(axis="x", bottom=False)
        ax.set_xlabel(None)
    
    pass

def main():
    """Run as script."""

    or_stats = pd.read_csv(_FILE_IN, sep="\t")
    
    df0 = or_stats[or_stats["region"] == "Whole transcript"]
    df1 = or_stats[or_stats["region"] == "NMD target"]
    df2 = or_stats[or_stats["region"] == "Start proximal"]
    df3 = or_stats[or_stats["region"] == "Long exon"]
    df4 = or_stats[or_stats["region"] == "Distal"]

    pass


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
