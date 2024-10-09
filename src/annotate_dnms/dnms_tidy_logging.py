"""Log summary information for tidied DNMs."""

import logging
from pathlib import Path

import pandas as pd

import src

_FILE_IN = "data/interim/dnms_38_combined_af_vep_tidy_dedup.tsv"
_NAMES = "chr pos ref alt csq enst cohort id".split()

logger = logging.getLogger(__name__)


def read_tidy_dnms(path=_FILE_IN):
    """Read tidied DNMs."""

    return pd.read_csv(path, sep="\t", header=None, names=_NAMES)


def main():
    """Run as script."""

    df = read_tidy_dnms(_FILE_IN)

    logger.info(f"DNMs after tidying: {len(df)}")
    logger.info(f"Unique variants {len(df.drop_duplicates(['chr','pos','ref','alt']))}")
    logger.info(
        f"Unique by variant / ID: {len(df.drop_duplicates(['chr','pos','ref','alt','id']))}"
    )
    logger.info(f"Consequence value counts:\n{df.csq.value_counts()}")
    logger.info(f"Cohort value counts:\n{df.cohort.value_counts()}")
    logger.info(f"Unique identifiers: {len(df[['cohort','id']].drop_duplicates())}")
    return df


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
