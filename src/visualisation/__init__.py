"""Visualisation tools for regional nonsense constraint paper."""

# Imports
from collections import namedtuple
from pathlib import Path

import colorsys
import matplotlib.colors as mc
import matplotlib.pyplot as plt
import seaborn as sns

from src import constants as C
from src import setup_logger


# Functions
def color_palette(style="default"):
    """Choose a color palette."""

    if not style in ["default","regions", "enrichment"]:
        raise ValueError("style must be one of 'default', 'regions', or 'enrichment'.")
    
    if style == "default":
        plt.style.use(C.COLOR_VIBRANT)

    if style == "regions":
        plt.style.use(C.COLOR_REGIONS)

    if style == "enrichment":
        plt.style.use(C.COLOR_ENRICHMENT)

    return sns.color_palette().as_hex()


def adjust_lightness(color, amount=0.5):
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))

    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


def panel_label(ax, s, x=-0.05, y=1.05, **kwargs):
    ax.text(
        x,
        y,
        s,
        transform=ax.transAxes,
        va="bottom",
        ha="right",
        fontsize=8,
        fontweight="bold",
        **kwargs,
    )

