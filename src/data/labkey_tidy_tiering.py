"""Tidy the tiering data, and find the highest tier for each participant."""

# Imports
from pathlib import Path

import pandas as pd

from src import setup_logger
from src import constants as C
from src.data.labkey_tidy_exit_questionnaires import tidy_gel_chr_names

# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def read_tiers(path="../data/tiering_data_2023-08-01_09-44-03.tsv"):
    """Get tiering data."""

    logger.info("Reading tiering data.")

    tiers = (
        pd.read_csv(
            path,
            sep="\t",
            usecols=[
                "Participant Id",
                "Assembly",
                "Chromosome",
                "Position",
                "Reference",
                "Alternate",
                "Tier",
            ],
            dtype={
                "Assembly": "category",
                "Reference": "category",
                "Alternate": "category",
                "Tier": "category",
            },
        )
        .set_axis(
            ["participant_id", "assembly", "chr", "pos", "ref", "alt", "tier"], axis=1
        )
        .drop_duplicates()
    ).pipe(tidy_gel_chr_names)

    logger.info(f"Participants with tiering data: {tiers.participant_id.nunique()}")
    logger.info(
        f"Tiered variants value counts:\n{tiers.groupby('assembly').tier.value_counts()}"
    )

    return tiers


def tidy_conflicting_tiers(df):
    """For a variant with mulitple tiers, assign it to the highest tier."""

    logger.info("Tidying conflicting tiers.")

    # Tiers as integer
    df["tier"] = df["tier"].str.slice(4).astype(int)

    # When a variant belongs to different tiers, assign it to the highest tier
    df["tier"] = (
        df.groupby(["participant_id", "assembly", "chr", "pos", "ref", "alt"])
        .tier.transform(min)
    )

    # Drop the resulting duplicate rows
    df = df.drop_duplicates()

    # Print summary information
    logger.info(
        f"Tiered variants value counts:\n{df.groupby('assembly').tier.value_counts()}"
    )

    return df


def get_max_tiers(df):
    """Get the highest tiered variant per participant.

    This is across both assemblies (GRCh37 and GRCh38).
    """

    logger.info("Getting highest tiers per participant.")

    max_tier = df.groupby("participant_id").tier.min().rename("max_tier").reset_index()

    logger.info(f"Highest tier value counts:\n{max_tier.max_tier.value_counts()}")

    return max_tier


def main():
    """Run as script."""

    tiers = read_tiers(C.LABKEY_TIERING).pipe(tidy_conflicting_tiers)
    max_tiers = get_max_tiers(tiers)

    logger.info("Writing to output.")
    tiers.to_csv(C.LABKEY_TIERS_CLEAN, sep="\t", index=False)
    max_tiers.to_csv(C.LABKEY_TIERS_HIGH, sep="\t", index=False)


if __name__ == "__main__":
    main()
