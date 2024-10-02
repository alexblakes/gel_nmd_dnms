"""Drop duplicate variants across DDD and GEL.

Where a DNV appears both in the DDD and GEL cohorts, drop it from the DDD cohort."""

import logging
from pathlib import Path

import pandas as pd

import src

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FILE_IN = "data/interim/dnms_38_combined_af_vep_tidy.tsv"
_FILE_OUT = "data/interim/dnms_38_combined_af_vep_tidy_dedup.tsv"

logger = logging.getLogger(__name__)


def read_data(path):
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chr", "pos", "ref", "alt", "csq", "enst", "cohort", "id"],
    )

    logger.info(f"Starting DNMs: {len(df)}")
    logger.info(f"DNMs per cohort:\n{df.cohort.value_counts()}")

    return df


def get_cohort_dnms(df: pd.DataFrame, cohort: str):
    return df.loc[df["cohort"] == cohort, ["chr", "pos", "ref", "alt"]]


def find_ddd_dnms_in_gel(df):
    gel = get_cohort_dnms(df, "GEL")
    ddd = get_cohort_dnms(df, "DDD")

    # Find DDD DNMs present in GEL. Keep the original indices of the DDD DNMs.
    dups = ddd.reset_index().merge(gel, how="inner").set_index("index")

    # Drop DDD DNMs with duplicate indices (i.e. with multiple matches in GEL)
    dups = dups.loc[~dups.index.duplicated(keep="first")]

    logger.info(f"DDD DNMs present in GEL: {len(dups)}")
    logger.info(f"DDD DNMs present in GEL (unique): {len(dups.drop_duplicates())}")

    return dups


def write_out(df, path):
    logger.info("Writing to output.")
    df.to_csv(path, sep="\t", header=False, index=False)
    return df


def drop_ddd_dnms_in_gel(df):
    drop_index = find_ddd_dnms_in_gel(df).index

    df = df.loc[df.index.difference(drop_index)]

    logger.info(f"Remaining DNMs in total: {len(df)}")
    logger.info(f"Remaining DNMs per cohort:\n{df.cohort.value_counts()}")

    return df


def main():
    """Run as script."""

    df = read_data(_FILE_IN).pipe(drop_ddd_dnms_in_gel).pipe(write_out, _FILE_OUT)

    return df


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
