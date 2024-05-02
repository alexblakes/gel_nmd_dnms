"""Extract 'case solved' and ACMG data from the exit questionnaires."""

# Imports
from pathlib import Path

import numpy as np
import pandas as pd

from src import setup_logger
from src import constants as C


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def read_exit_questionnaire(path):
    """Read exit questionnaire data"""
    
    logger.info("Reading exit questionnaire data.")

    eq = (
        pd.read_csv(
            path,
            sep="\t",
            header=0,
            usecols=[
                "Participant Id",
                "Event Date",
                "Case Solved Family",
                "Acmg Classification",
                "Assembly",
                "Chromosome",
                "Position",
                "Reference",
                "Alternate",
            ],
            dtype={"Position": "Int64"},
        )
        .set_axis(
            [
                "participant_id",
                "date",
                "case_solved",
                "acmg",
                "assembly",
                "chr",
                "pos",
                "ref",
                "alt",
            ],
            axis=1,
        )
        .drop_duplicates()
    )

    logger.info(f"Exit questionnaires: {len(eq)}")

    return eq


def get_case_solved(df):
    """Get the "case solved" status for each participant.

    Where there are conflicting outcomes, use the most recent one.
    """

    logger.info("Getting 'case solved' annotations.")

    case_solved = (
        df[["participant_id", "case_solved", "date"]]
        .sort_values("date")
        .drop_duplicates(["participant_id"], keep="last")
        .drop("date", axis=1)
    )

    logger.info(f"Participants with an exit questionnaire: {len(case_solved)}")
    logger.info(f"Case solved value counts:\n{case_solved.case_solved.value_counts(dropna=False)}")

    return case_solved


def tidy_gel_chr_names(df):
    """Tidy the chr names for GRCh37 and 38 variants in GEL.

    For many GEL datasets, GRCh37 variants lack a "chr" prefix. Some GRCh38
    variants lack a "chr" prefix. Chroms "M" and "MT" are used interchangeably.
    """

    # Fix mitochondrial nomenclature
    df["chr"] = df.chr.replace("MT", "M")

    # Add chr prefix for all GRCh38 chromosomes
    m1 = df["assembly"] == "GRCh38"
    m2 = df.chr.str.startswith("chr")
    df.loc[m1 & ~m2, "chr"] = "chr" + df.chr

    return df


def get_exit_questionnaire_acmg(df):
    """Get ACMG annotations for variants in exit questionnaires."""

    logger.info(f"Getting ACMG annotations.")

    acmg = (
        df.drop("case_solved", axis=1)
        .replace("not_assessed", np.nan)
        .dropna()
        .pipe(tidy_gel_chr_names)
        .drop_duplicates()
    )

    logger.info(f"ACMG annotations: {len(acmg)}")
    logger.info(f"ACMG annotations by assembly:\n{acmg.groupby('assembly').chr.count()}")
    logger.info(f"ACMG value counts:\n{acmg.groupby('assembly').acmg.value_counts()}")

    id_37 = set(acmg.loc[acmg.assembly == "GRCh37", "participant_id"])
    id_38 = set(acmg.loc[acmg.assembly == "GRCh38", "participant_id"])
    both = id_37.intersection(id_38)
    logger.info(f"Participants with both GRCh37 & 38 assemblies: {len(both)}")

    return acmg


def main():
    """Run as script"""
    
    eq = read_exit_questionnaire(C.LABKEY_EXIT_QUESTIONNAIRE)

    case_solved = get_case_solved(eq)
    acmg = get_exit_questionnaire_acmg(eq)

    logger.info("Writing to output.")
    case_solved.to_csv(C.LABKEY_EQ_CASE_SOLVED, sep="\t", index=False)
    acmg.to_csv(C.LABKEY_EQ_ACMG, sep="\t", index=False)

    pass

if __name__ == "__main__":
    main()
