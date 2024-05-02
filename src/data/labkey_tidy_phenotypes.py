"""Aggregate HPO terms per participant."""

# Imports
from pathlib import Path

import pandas as pd

from src import setup_logger
from src import constants as C


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def read_hpo(path):
    """Get HPO terms, which are present in each participant, from LabKey."""

    logger.info("Reading participant phenotype (HPO) data.")

    hpo = (
        pd.read_csv(
            path,
            sep="\t",
            usecols=["Participant Id", "Hpo Term", "Hpo Present"],
        )
        .set_axis(["participant_id", "hpo", "hpo_present"], axis=1)
        .query("hpo_present == 'Yes'")
        .drop("hpo_present", axis=1)
    )

    logger.info(f"Participants with >=1 HPO term: {hpo.participant_id.nunique()}")
    logger.info(f"Participants with duplicated HPO terms: {hpo.duplicated().sum()}")
    logger.info(f"HPO terms per participant:\n{hpo.groupby('participant_id').hpo.count().agg(['mean', 'median', 'min', 'max'])}")
    logger.info(f"The ten most common HPO terms in the cohort:\n{hpo.hpo.value_counts().head(10)}")

    return hpo


def aggregate_hpo_terms_per_participant(hpo):
    """Join all HPO terms per participant."""

    hpo_joined = hpo.groupby("participant_id")["hpo"].apply(",".join).reset_index()
    
    return hpo_joined


def main():
    """Run as script."""

    hpo = read_hpo(C.LABKEY_PHENOTYPES).pipe(aggregate_hpo_terms_per_participant)

    logger.info("Writing to output.")
    hpo.to_csv(C.LABKEY_PHENOTYPES_CLEAN, sep="\t", index=False)


if __name__ == "__main__":
    main()
