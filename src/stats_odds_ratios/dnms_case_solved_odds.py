"""Case solved odds ratio statistics."""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import contingency

import src
from src import stats_enrichment

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FILE_IN = "data/interim/dnms_annotated_clinical.tsv"
_FILE_OUT = "data/statistics/case_solved_odds_ratios.tsv"

logger = logging.getLogger(__name__)


def get_dnms(path):
    df = (
        pd.read_csv(
            path,
            sep="\t",
            usecols=[
                "participant_id",
                "omim_inheritance_simple",
                "region",
                "constraint",
                "cohort",
                "csq",
                "case_solved",
            ],
        )
        .query("cohort=='GEL'")
        .drop("cohort", axis=1)
        .fillna({"omim_inheritance_simple": ""})
        .assign(
            case_solved=lambda x: np.where(
                x.case_solved.isin(["yes", "partially"]), True, False
            )
        )
    )

    return df


def get_solved_count(df):
    solved = df["case_solved"].sum()
    unsolved = df["case_solved"].count() - solved
    return solved, unsolved


def get_odds_ratio(df, conditions, region, constraint):
    # Get cases
    case = df[conditions].drop_duplicates("participant_id")

    # Get controls
    case_ids = case["participant_id"]
    control = df[~df["participant_id"].isin(case_ids)].drop_duplicates("participant_id")

    # Create a contingency table
    table = pd.DataFrame(index=["solved", "unsolved"], columns=["case", "control"])
    table["case"] = get_solved_count(case)
    table["control"] = get_solved_count(control)

    # Get odds ratio
    result = contingency.odds_ratio(table)
    _or = result.statistic
    ci_lo, ci_hi = result.confidence_interval(0.95)

    stats = pd.DataFrame(
        index=[[region], [constraint]],
        columns=["odds_ratio", "ci_lo", "ci_hi"],
        data=[[_or, ci_lo, ci_hi]],
    )
    stats.index.set_names(["region", "constraint"], inplace=True)

    return stats


def main():
    df = get_dnms(_FILE_IN)

    # Masks
    ## Consequences
    m1 = df["csq"].isin(["frameshift_variant", "stop_gained"])

    ## Constraint
    m2 = df["constraint"] == "unconstrained"
    m3 = df["constraint"] == "constrained"

    ## NMD region
    m4 = df["region"] == "nmd_target"
    m5 = df["region"] == "long_exon"
    m6 = df["region"] == "distal_nmd"
    m7 = df["region"] == "start_proximal"

    # Get the odds of case solved for each group
    or_stats = (
        pd.concat(
            [
                # get_odds_ratio(df, (m1), "whole_transcript", "all"),
                get_odds_ratio(df, (m1 & m2), "whole_transcript", "unconstrained"),
                get_odds_ratio(df, (m1 & m3), "whole_transcript", "constrained"),
                # get_odds_ratio(df, (m1 & m4), "nmd_target", "all"),
                get_odds_ratio(df, (m1 & m4 & m2), "nmd_target", "unconstrained"),
                get_odds_ratio(df, (m1 & m4 & m3), "nmd_target", "constrained"),
                # get_odds_ratio(df, (m1 & m5), "long_exon", "all"),
                get_odds_ratio(df, (m1 & m5 & m2), "long_exon", "unconstrained"),
                get_odds_ratio(df, (m1 & m5 & m3), "long_exon", "constrained"),
                # get_odds_ratio(df, (m1 & m6), "distal", "all"),
                get_odds_ratio(df, (m1 & m6 & m2), "distal", "unconstrained"),
                get_odds_ratio(df, (m1 & m6 & m3), "distal", "constrained"),
                # get_odds_ratio(df, (m1 & m7),"start_proximal", "all"),
                get_odds_ratio(df, (m1 & m7 & m2), "start_proximal", "unconstrained"),
                get_odds_ratio(df, (m1 & m7 & m3), "start_proximal", "constrained"),
            ]
        )
        .dropna()
        .reset_index()
        .pipe(stats_enrichment.sort_region_column)
    )

    # Write to output
    or_stats.to_csv(_FILE_OUT, sep="\t", index=False)

    return or_stats  #! Testing


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
