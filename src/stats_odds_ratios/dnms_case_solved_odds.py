"""Case solved odds ratio statistics."""

import logging

import numpy as np
import pandas as pd
from statsmodels.stats import contingency_tables as ct

import src
from src import stats_enrichment as sen

_FILE_IN = "data/interim/dnms_annotated_clinical.tsv"
_FILE_OUT = "data/statistics/case_solved_odds_ratios.tsv"
_OR_LABELS = [
    "Synonymous (Full CDS)",
    "Missense (Full CDS)",
    "PTV (Full CDS)",
    "PTV (NMD target)",
    "PTV (Start proximal)",
    "PTV (Long exon)",
    "PTV (Distal)",
]
_CONSTRAINT_LABEL = {x: x.capitalize() for x in ["constrained", "unconstrained"]}

logger = logging.getLogger(__name__)


def read_data(path):
    return pd.read_csv(
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


def tidy_data(df):
    return (
        df.query("cohort=='GEL'")
        .drop("cohort", axis=1)
        .fillna({"omim_inheritance_simple": ""})
        .assign(
            case_solved=lambda x: np.where(
                x.case_solved.isin(["yes", "partially"]), True, False
            )
        )
    )


def get_cases_and_controls(df, masks):
    cases = df[masks].drop_duplicates("participant_id").assign(cohort="case")
    cases_ids = cases["participant_id"]
    controls = (
        df[~df["participant_id"].isin(cases_ids)]
        .drop_duplicates("participant_id")
        .assign(cohort="control")
    )

    return pd.concat([cases, controls])


def count_solved_unsolved(df):
    return (
        df.groupby("cohort")
        .agg(
            solved=("case_solved", "sum"),
            unsolved=("case_solved", lambda x: x.count() - x.sum()),
        )
        .T
    )


def get_odds_ratio(df):
    table = ct.Table2x2(np.array(df))
    _or = table.oddsratio
    ci_lo, ci_hi = table.oddsratio_confint()
    err_lo = _or - ci_lo
    err_hi = ci_hi - _or
    p = table.oddsratio_pvalue()

    return pd.DataFrame(
        columns=["odds_ratio", "err_lo", "err_hi", "p"], data=[[_or, err_lo, err_hi, p]]
    )


def label(df, **labels):
    return df.assign(**labels).set_index(list(labels.keys()))


def or_pipeline(df, masks, **labels):
    return (
        df.pipe(get_cases_and_controls, (masks))
        .pipe(count_solved_unsolved)
        .pipe(get_odds_ratio)
        .pipe(label, **labels)
    )


def bonferroni_correction(df, alpha=0.05):
    logger.info(f"Bonferroni P value: {alpha / len(df)}")
    return df.assign(bfr_sig = lambda x: x["p"] < alpha / len(x))


def sort_data(df):
    return (
        df.reset_index()
        .pipe(sen.sort_column, "constraint", _CONSTRAINT_LABEL)
        .pipe(sen.sort_column, "label", _OR_LABELS)
        .sort_values(["label", "constraint"])
        .set_index(["csq", "region", "constraint", "label"])
    )


def write_out(df, path):
    df.to_csv(path, sep="\t")
    return df


def main():

    df = read_data(_FILE_IN).pipe(tidy_data)

    # Masks
    synonymous = df["csq"] == "synonymous_variant"
    missense = df["csq"] == "missense_variant"
    ptv = df["csq"].isin(["frameshift_variant", "stop_gained"])
    unconstrained = df["constraint"] == "unconstrained"
    constrained = df["constraint"] == "constrained"
    nmd_target = df["region"] == "nmd_target"
    long_exon = df["region"] == "long_exon"
    distal_nmd = df["region"] == "distal_nmd"
    start_proximal = df["region"] == "start_proximal"

    constraint = [unconstrained, constrained]
    masks = [
        synonymous,
        missense,
        ptv,
        ptv & nmd_target,
        ptv & start_proximal,
        ptv & long_exon,
        ptv & distal_nmd,
    ]
    masks = [m & c for m in masks for c in constraint]

    # Labels
    csqs = ["synonymous"] * 2 + ["missense"] * 2 + ["ptv"] * 10
    regions = (
        ["transcript"] * 6
        + ["nmd_target"] * 2
        + ["start_proximal"] * 2
        + ["long_exon"] * 2
        + ["distal_nmd"] * 2
    )
    constraint_levels = ["unconstrained", "constrained"] * 7
    or_labels = [x for y in zip(_OR_LABELS, _OR_LABELS) for x in y]
    labels = zip(csqs, regions, constraint_levels, or_labels)
    header = ["csq", "region", "constraint", "label"]
    label_dicts = [{h: l for h, l in zip(header, _label)} for _label in labels]

    return (
        pd.concat([or_pipeline(df, m, **ld) for m, ld in zip(masks, label_dicts)])
        .pipe(bonferroni_correction)
        .pipe(sort_data)
        .pipe(write_out, _FILE_OUT)
    )


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
