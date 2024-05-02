"""Module docstring."""

import logging
from pathlib import Path

import pandas as pd

import src
from src import constants as C

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"

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


def main():
    """Run as script."""
    ptvs = read_dnm_ptvs(C.DNMS_ANNOTATED_CLINICAL)

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

    df = pd.concat([count_dnms_per_region(ptvs[m], n) for m, n in zipped], axis=1).melt(
        ignore_index=False, var_name="mask"
    )

    # counts = count_dnms_per_region(ptvs)
    return df


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
