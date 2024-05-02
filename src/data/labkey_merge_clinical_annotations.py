"""Docstring."""

# Imports
from pathlib import Path

import pandas as pd

from src import setup_logger
from src import constants as C


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def read_participant_data(path):
    """Read participant data."""

    logger.info("Reading participant data.")

    df = pd.read_csv(
        path,
        sep="\t",
        low_memory=False,
        usecols=[
            "Participant Id",
            "Handling Gmc Trust",
            "Year Of Birth",
            "Participant Type",
        ],
    ).set_axis(["participant_id", "gmc", "birth_year", "participant_type"], axis=1)

    # Check there are no duplicate participants
    logger.info(f"Unique participant IDs: {df['participant_id'].nunique()}")
    logger.info(f"Duplicated participant IDs: {df.duplicated('participant_id').sum()}")
    logger.info(f"GMC value counts:\n{df.gmc.value_counts()}")
    logger.info(
        f"Participant type value counts:\n{df['participant_type'].value_counts()}"
    )

    return df


def main():
    """Run as script."""

    case_solved = pd.read_csv(C.LABKEY_EQ_CASE_SOLVED, sep="\t")
    max_tiers = pd.read_csv(C.LABKEY_TIERS_HIGH, sep="\t")
    hpo = pd.read_csv(C.LABKEY_PHENOTYPES_CLEAN, sep="\t")
    participant = read_participant_data(C.LABKEY_PARTICIPANTS)

    logger.info("Merging clinical annotations.")

    clinical = (
        case_solved.merge(max_tiers, how="outer", on="participant_id")
        .merge(hpo, how="outer", on="participant_id")
        .merge(participant, how="left")
    )

    logger.info(
        f"Duplicated participant IDs: {clinical.participant_id.duplicated().sum()}"
    )
    logger.info(f"Non-null values (all participants):\n{clinical.notna().sum()}")
    logger.info(
        f"Non-null values (probands):\n{clinical[clinical.participant_type == 'Proband'].notna().sum()}"
    )

    logger.info("Writing to output.")
    clinical.to_csv(C.LABKEY_CLINICAL, sep="\t", index=False)

    pass


if __name__ == "__main__":
    main()
