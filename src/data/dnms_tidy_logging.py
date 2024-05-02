"""Log summary information for tidied DNMs."""

# Imports
from pathlib import Path

import numpy as np
import pandas as pd

from src import setup_logger
from src import constants as C


# Logging
logger = setup_logger(Path(__file__).stem)


# Module constants
_NAMES = "chr pos ref alt csq enst cohort id".split()


# Functions
def read_tidy_dnms(path):
    """Read tidied DNMs."""

    return pd.read_csv(path, sep="\t", header=None, names=_NAMES)


def main():
    """Run as script."""
    
    df = read_tidy_dnms(C.DNMS_VEP_TIDY)

    logger.info(f"DNMs after tidying: {len(df)}")
    logger.info(f"Unique variants {len(df.drop_duplicates(['chr','pos','ref','alt']))}")
    logger.info(
        f"Unique by variant / ID: {len(df.drop_duplicates(['chr','pos','ref','alt','id']))}"
    )
    logger.info(f"Consequence value counts:\n{df.csq.value_counts()}")
    logger.info(f"Cohort value counts:\n{df.cohort.value_counts()}")
    return df


if __name__ == "__main__":
    main()
