"""Calculate relative enrichment of de novo variants."""

import pandas as pd

from src import constants as C

_REGION_LABELS = {
    "transcript": "Full CDS",
    "nmd_target": "NMD target",
    "start_proximal": "Start proximal",
    "long_exon": "Long exon",
    "distal_nmd": "Distal",
}


def sort_column(df, column="region", labels=_REGION_LABELS, **kwargs):
    """Rename and reorder the column of a dataframe."""

    kwargs.setdefault("ordered", True)

    assert all(
        i in labels for i in df[column]
    ), "Some values in the column are missing from the labels."

    df[column] = pd.Categorical(
        df[column], categories=labels, **kwargs
    ).rename_categories(labels)

    return df.sort_values(column)


def sort_index(df, labels=_REGION_LABELS, **kwargs):
    """Rename and reorder the index of a series or dataframe."""

    kwargs.setdefault("name", "region")

    assert all(
        i in labels.keys() for i in df.index
    ), "Some values in the index are missing from the labels."

    index = pd.Index(labels, **kwargs)

    return df.reindex(index).rename(index=labels)


def categorical_regions_column(series, categories=C.REGIONS, labels=C.REGION_LABELS):
    return pd.Categorical(
        series,
        categories=categories,
        ordered=True,
    ).rename_categories(labels)


def sort_region_column(df, column="region", **kwargs):
    df[column] = categorical_regions_column(df[column], **kwargs)
    return df.sort_values(column)


def categorical_regions_index(index, categories, labels, name="region"):
    return pd.CategoricalIndex(
        index,
        categories=categories,
        ordered=True,
        name=name,
    ).rename_categories(labels)


def sort_region_index(df, **kwargs):
    df.index = categorical_regions_index(df.index, **kwargs)
    return df.sort_index()
