"""Merge clinical annotations from Labkey."""

import logging
from pathlib import Path

import pandas as pd

import src

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_CASE_SOLVED = "data/interim/labkey_exit_questionnaires_case_solved.tsv"
_HIGH_TIERS = "data/interim/labkey_tiering_highest_tiers.tsv"
_PHENOTYPES_CLEAN = "data/interim/labkey_phenotypes_clean.tsv"
_PARTICIPANTS = "data/raw/participant_2023-08-09_12-18-36.tsv"
_FILE_OUT = "data/interim/labkey_participant_clinical.tsv"

logger = logging.getLogger(__name__)


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

    case_solved = pd.read_csv(_CASE_SOLVED, sep="\t")
    max_tiers = pd.read_csv(_HIGH_TIERS, sep="\t")
    hpo = pd.read_csv(_PHENOTYPES_CLEAN, sep="\t")
    participant = read_participant_data(_PARTICIPANTS)

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
    clinical.to_csv(_FILE_OUT, sep="\t", index=False)

    pass


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
