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
    fig = plt.figure(figsize=(10 * C.CM, 18 * C.CM), layout="constrained")
    subfigs = fig.subfigures(3, 1, height_ratios=[5,8,4], hspace=0.05)
    # subsubfigs = subfigs[0].subfigures(1, 2)
    axs_middle = subfigs[1].subplots(2, 2).flatten()
    axs_top_left = subfigs[0].subplots(7, 1, sharex=True).flatten()
    axs_bottom = subfigs[2].subplots(1, 2).flatten()

    # Odds ratio plots
    por.plot(axs_top_left)
    subfigs[0].get_layout_engine().set(h_pad=0)

    # Enrichment plots
    paths = [
        "data/statistics/dnms_enrichment_all_genes_tidy.tsv",
        "data/statistics/dnms_enrichment_ad_tidy.tsv",
        "data/statistics/dnms_enrichment_ar_tidy.tsv",
        "data/statistics/dnms_enrichment_non_morbid_tidy.tsv",
    ]
    titles = ["All genes", "OMIM morbid (AD)", "OMIM morbid (AR)", "Non-morbid"]

    for path, title, ax in zip(paths, titles, axs_middle):
        pse.plot(path, title, ax)

    vis.same_lims(axs_middle, y=True)

    # Recurrent dnPTV plots
    # axs_bottom[2].axis("off")
    # axs_bottom[3].axis("off")
    recurrent_dnms = prd.read_data()
    columns = ["nmd_target", "one_nmd_escape"]
    colors = [sns.color_palette()[x] for x in [3, 5]]
    xlabels = [
        "$\it{De}$ $\it{novo}$ nonsense & FS variants\nin a constrained NMD target\nregion in a non-morbid gene",
        "$\it{De}$ $\it{novo}$ nonsense & FS variants\nin a constrained NMD escape\nregion in a non-morbid gene",
    ]

    for ax, column, xlabel, color in zip(axs_bottom, columns, xlabels, colors):
        vis.vertical_bars(recurrent_dnms[column], ax=ax, color=color)
        prd.customise_axes(ax, ylabel="Transcripts", xlabel=xlabel)

    vis.same_lims(axs_bottom, y=True)

    # Add panel labels
    labels = list("abcdefg")
    axs = axs_top_left[0], axs_middle[0], axs_middle[1], axs_middle[2], axs_middle[3], axs_bottom[0], axs_bottom[1]
    for ax, l in zip(axs, labels):
        vis.panel_label(ax, l)

    plt.savefig(PNG, dpi=600)
    plt.savefig(SVG)
    plt.close("all")

    return fig


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
