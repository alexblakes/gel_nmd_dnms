"""Plot the number of transcripts with recurrent dnPTVs per region."""

import logging
from pathlib import Path

import pandas as pd

import src

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FILE_IN = "data/interim/dnms_annotated_clinical.tsv"
_FILE_OUT = "data/statistics/recurrent_dnms.tsv"

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
    )


def filter_for_ptvs(df):
    return df.query("csq == 'frameshift_variant' | csq == 'stop_gained'")


def filter_regions(df, masks=None):
    if not masks:
        return df

    else:
        return df.loc[masks, :]


def group_nmd_escape_regions(df, region=None):
    if not region:
        return df
    else:
        return df.assign(region=region)


def count_dnms_per_region(df, max_per_transcript=False):
    """Calculate the number of dnPTVs per region per transcript."""

    count = df.groupby(["enst", "region"])["chr"].count()

    count_of_counts = (
        lambda x: x.value_counts().reset_index().set_axis(["n", "count"], axis=1)
    )

    if max_per_transcript:
        return count.groupby("enst").max().pipe(count_of_counts)

    return count.pipe(count_of_counts)


def clip_dnm_counts(df, clip=5):
    return df.assign(n=lambda x: x.n.clip(upper=clip))


def recount_after_clipping(df):
    return df.groupby("n")["count"].sum()


def combine_dnm_counts(ptvs):
    constrained = lambda x: x.constraint == "constrained"
    non_morbid = lambda x: x.omim_inheritance_simple.isna()
    nmd_target = lambda x: x.region == "nmd_target"
    nmd_escape = lambda x: x.region.isin(["start_proximal","long_exon","distal_nmd"])

    filters = [
        None,
        lambda x: constrained(x),
        lambda x: constrained(x) & non_morbid(x),
        lambda x: constrained(x) & non_morbid(x) & nmd_target(x),
        lambda x: constrained(x) & non_morbid(x) & nmd_escape(x),
        lambda x: constrained(x) & non_morbid(x) & nmd_escape(x),
    ]

    rename_region = [None, None, None, None, "nmd_escape", None]

    max_per_transcript = [False, False, False, False, False, True]

    titles = [
        "any_region",
        "constrained",
        "non_morbid",
        "nmd_target",
        "any_nmd_escape",
        "one_nmd_escape",
    ]

    return pd.concat(
        [
            ptvs.pipe(filter_regions, filter)
            .pipe(group_nmd_escape_regions, rename)
            .pipe(count_dnms_per_region, _max)
            .pipe(clip_dnm_counts)
            .pipe(recount_after_clipping)
            .rename(title)
            for filter, rename, _max, title in zip(
                filters, rename_region, max_per_transcript, titles
            )
        ],
        axis=1,
    )

def write_out(df, path):
    df.to_csv(path, sep="\t")
    return df

def main():
    """Run as script."""

    return read_dnm_ptvs(_FILE_IN).pipe(filter_for_ptvs).pipe(combine_dnm_counts).pipe(write_out, _FILE_OUT)


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
