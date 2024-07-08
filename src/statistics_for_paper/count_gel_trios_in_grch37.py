# Count the trio offspring in GRCh37 in GEL.

import logging
from pathlib import Path

import pandas as pd

import src

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FILE_IN="data/interim/gel_dnm_offspring_clean.tsv"

logger = logging.getLogger(__name__)

def main():
    df = pd.read_csv(_FILE_IN, sep="\t")

    logger.info(f"Participants with at least one DNM: {df['Base Filter Total'].count()}")
    logger.info(f"Assembly value counts:\n{df['Assembly'].value_counts()}")

    return df

if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()