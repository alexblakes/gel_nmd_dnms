"""Module docstring."""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


_LOGFILE = f"data/logs/{Path(__file__).stem}.log"

logger = logging.getLogger(__name__)


def grouped_bars(data, ax=None):
    """Plot a grouped bar chart.
    
    Parameters:
        data : Series or single-column data frame of values for plotting. Must have a 
            two-level multi-index. Level 0 is the bars, level 1 is the xticks.
        ax : The Axes on which to plot. Defaults to None.                    
    """    
    
    if not ax:
        ax = plt.gca()
    
    outer_level_values = data.index.get_level_values(0).unique()

    for i, level_value in enumerate(outer_level_values):
        multiplier = i
        data_subset = data.xs(level_value, level=0).squeeze()
        n_groups = len(data_subset)
        x = np.arange(n_groups)
        n_bars = len(outer_level_values)
        bar_width = 1 / (n_bars + 1)
        bar_offset = bar_width * multiplier
        tick_offset = bar_width * ((n_bars / 2) - 0.5)

        ax.bar(x + bar_offset, data_subset, bar_width)
        ax.set_xticks(x + tick_offset, labels=data_subset.index)

    return None

def annotate_grouped_bars(ax=None, bar_label=True, legend=False, legend_kwargs={}, set_kwargs={}):
    """Customise grouped bar plot.
    
    Parameters:
        ax : The Axes on which to plot. Defaults to none.
        bar_label : Bool. Whether to add bar labels.
        legend : List of strings for bar legend. Defaults to False.
        **set_kwargs : kwargs passed to ax.set()
    """

    if not ax:
        ax = plt.gca()
    
    bars = [c for c in ax.containers]

    if bar_label:
        for b in bars:
            ax.bar_label(b)

    if legend:
        legend_bars = bars[0:len(bars)]
        assert len(legend) == len(legend_bars),\
            "The number of items in 'legend' does not equal the number of bars"
        ax.legend(legend_bars, legend, **legend_kwargs)

    ax.set(**set_kwargs)

    return None