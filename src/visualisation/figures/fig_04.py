"""Plot Figure 4"""

import logging

from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns

import src
from src import constants as C
from src import visualisation as vis
from src.stats_odds_ratios import plot_odds_ratios as por
from src.stats_enrichment import plot as pse
from src.stats_recurrent_dnms import plot as prd

PNG = "data/plots/figure_04.png"
SVG = "data/plots/figure_04.svg"

logger = logging.getLogger(__name__)


def main():
    """Run as script."""

    plt.style.use([C.STYLE_DEFAULT, C.COLOR_ENRICHMENT])
    fig = plt.figure(figsize=(18 * C.CM, 10 * C.CM), layout="constrained")
    subfigs = fig.subfigures(2, 1, height_ratios=[4, 3], hspace=0.1)
    subsubfigs = subfigs[1].subfigures(1, 2, width_ratios=[3, 4])
    axs_top = subfigs[0].subplots(1, 3).flatten()
    axs_bottom_left = subsubfigs[0].subplots(7, 1, sharex=True).flatten()
    axs_bottom_right = subsubfigs[1].subplots(1, 2).flatten()

    # Enrichment plots
    paths = [
        "data/statistics/dnms_enrichment_all_genes_tidy.tsv",
        "data/statistics/dnms_enrichment_ad_tidy.tsv",
        "data/statistics/dnms_enrichment_non_morbid_tidy.tsv",
    ]
    titles = ["All genes", "OMIM morbid (AD)", "Non-morbid"]

    for path, title, ax in zip(paths, titles, axs_top):
        pse.plot(path, title, ax)

    vis.same_lims(axs_top, y=True)

    # Odds ratio plots
    por.plot(axs_bottom_left)
    subsubfigs[0].get_layout_engine().set(h_pad=0)

    # Recurrent dnPTV plots
    recurrent_dnms = prd.read_data()
    columns = ["nmd_target", "one_nmd_escape"]
    colors = [sns.color_palette()[x] for x in [3, 5]]
    xlabels = [
        "$\it{dn}$PTVs in a constrained\nNMD target region\nin a non-morbid gene",
        "$\it{dn}$PTVs in a constrained\nNMD escape region\nin a non-morbid gene",
    ]

    for ax, column, xlabel, color in zip(axs_bottom_right, columns, xlabels, colors):
        vis.vertical_bars(recurrent_dnms[column], ax=ax, color=color)
        prd.customise_axes(ax, ylabel="Transcripts", xlabel=xlabel)

    vis.same_lims(axs_bottom_right, y=True)

    # Add panel labels
    labels = list("abcdef")
    axs = axs_top[0], axs_top[1], axs_top[2], axs_bottom_left[0], axs_bottom_right[0], axs_bottom_right[1]
    for ax, l in zip(axs, labels):
        vis.panel_label(ax, l)

    plt.savefig(PNG, dpi=600)
    plt.savefig(SVG)
    plt.close("all")

    return fig


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
