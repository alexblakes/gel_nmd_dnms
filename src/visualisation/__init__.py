"""Visualisation tools for regional nonsense constraint paper."""

import colorsys
import matplotlib.colors as mc
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def adjust_lightness(color, amount=0.5):
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))

    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


def adjust_alpha(color, alpha):
    return mc.to_rgba(color, alpha)


def panel_label(ax, s, x=-0.05, y=1.05, **kwargs):
    kwargs.setdefault("fontsize", 8)
    kwargs.setdefault("fontweight", "bold")
    kwargs.setdefault("transform", ax.transAxes)
    kwargs.setdefault("va", "bottom")
    kwargs.setdefault("ha", "right")

    ax.text(
        x,
        y,
        s,
        **kwargs
    )

    return ax


def vertical_bars(series, ax=None, **kwargs):
    """Vertical bar chart for values in a series.

    xticklabels are taken from the series' index.
    """

    kwargs.setdefault("tick_label", series.index)
    kwargs.setdefault("color", sns.color_palette())
    kwargs.setdefault("ecolor", [adjust_lightness(c, 0.8) for c in sns.color_palette()])

    if not ax:
        ax = plt.gca()

    n_bars = len(series)
    x = np.arange(n_bars)

    ax.bar(x=x, height=series, **kwargs)

    return ax


def horizontal_bars(
    values,
    ax=None,
    **kwargs
):

    kwargs.setdefault("tick_label", values.index)
    kwargs.setdefault("color", sns.color_palette())
    kwargs.setdefault("ecolor", [adjust_lightness(c, 0.8) for c in sns.color_palette()])

    if not ax:
        ax = plt.gca()

    n = len(values)  # Number of bars
    height = 1 - (1 / n)
    y = np.arange(n)

    ax.barh(y=y, width=values, height=height, **kwargs)

    return ax