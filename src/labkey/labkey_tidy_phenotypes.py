"""Aggregate HPO terms per participant."""

import logging
from pathlib import Path

import pandas as pd

import src
from src import constants as C

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FILE_IN = "data/raw/rare_diseases_participant_phen_2023-08-01_09-40-58.tsv"
_FILE_OUT = "data/interim/labkey_phenotypes_clean.tsv"

logger = logging.getLogger(__name__)


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

    hpo = read_hpo(_FILE_IN).pipe(aggregate_hpo_terms_per_participant)

    logger.info("Writing to output.")
    hpo.to_csv(_FILE_OUT, sep="\t", index=False)


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
