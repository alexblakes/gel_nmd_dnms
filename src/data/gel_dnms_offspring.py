"""Clean the GEL DNM cohort data.

Find unique offspring in the cohort.
For individuals with both GRCh37 and GRCh38 data, choose GRCh38. 
"""

# Import modules
from pathlib import Path

import pandas as pd

from src import setup_logger
from src import constants as C


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def get_offspring(cohort):
    """Filter for offspring."""

    offspring = cohort[cohort["Member"] == "Offspring"]

    logger.info(f"Number of offspring: {len(offspring)}")

    return offspring


def filter_unique_offspring(offspring):
    """Find unique individuals amongst the offspring.

    Each individual should be present in only one trio. Where an individual has data for
    both GRCH37 and GRCh38, drop the GRCh37 entry.
    """
    dup = offspring["Participant Id"][offspring.duplicated("Participant Id")]

    logger.info(f"Duplicated participant IDs: {len(dup)}")

    m1 = offspring["Participant Id"].isin(dup)
    m2 = offspring["Assembly"] == "GRCh37"

    logger.info(f"Duplicated IDs in GRCh37: {(m1 & m2).sum()}")

    unique_offspring = offspring[~(m1 & m2)]

    logger.info(
        f"Duplicated IDs after filtering: {unique_offspring.duplicated('Participant Id').sum()}"
    )
    logger.info(f"Remaining offspring: {len(unique_offspring)}")

    return unique_offspring


def main():
    """Run as script."""

    offspring = (
        pd.read_csv(C.GEL_DE_NOVO_COHORT, sep="\t")
        .pipe(get_offspring)
        .pipe(filter_unique_offspring)
    )

    logger.info("Writing to output.")
    offspring.to_csv(C.DE_NOVO_OFFSPRING, sep="\t", index=False)


if __name__ == "__main__":
    main()
