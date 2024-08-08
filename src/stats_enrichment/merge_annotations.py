"""Add constraint and OMIM annotations to DNMs."""

import logging

import numpy as np
import pandas as pd

import src
from src.merge_annotations import dnms_annotate_constraint as dac

_DNMS_ANNOTATED = "data/interim/dnms_annotated.tsv"
_REGIONAL_NONSENSE_CONSTRAINT = "data/uom_csf/regional_nonsense_constraint.tsv"
_REGIONAL_CONSTRAINT_STATS = "data/uom_csf/regional_constraint_stats.tsv"
_GENE_IDS = "data/uom_csf/gene_ids.tsv"
_GENEMAP2_SIMPLE = "data/uom_csf/genemap2_simple.tsv"
_FILE_OUT = "data/interim/dnms_enrichment_counts.tsv"

logger = logging.getLogger(__name__)


def read_dnms(path):
    return pd.read_csv(path, sep="\t", usecols=["enst", "region", "csq"])


def add_transcript_as_region(df):
    return pd.concat([df, df.assign(region="transcript")])


def read_constraint_stats(path):
    return pd.read_csv(path, sep="\t", usecols=["enst", "region", "csq", "n_exp"])


def read_nonsense_constraint(path):
    return pd.read_csv(path, sep="\t", usecols=["enst", "region", "constraint"])


def tidy_nonsense_constraint(df):
    return df.loc[df["constraint"] != "indeterminate"].dropna(subset="constraint")


def unify_truncating_variants(df):
    """Convert all PTV annotations to "truncating"."""

    return df.assign(
        csq=lambda x: np.where(
            x["csq"].isin(["stop_gained", "frameshift", "frameshift_variant"]),
            "truncating",
            x["csq"],
        )
    )


def count_dnms(df):
    return (
        df.groupby(["enst", "region", "csq"]).agg(n_dnms=("csq", "count")).reset_index()
    )


def reorder_data(df):
    return df[
        ["enst", "region", "constraint", "csq", "inheritance_simple", "n_exp", "n_dnms"]
    ]


def write_out(df, path):
    df.to_csv(path, sep="\t", index=False)
    return df


def main():
    """Run as script."""

    logging.disable()

    dnms = (
        read_dnms(_DNMS_ANNOTATED)
        .pipe(unify_truncating_variants)
        .pipe(add_transcript_as_region)
        .pipe(count_dnms)
    )
    constraint_stats = read_constraint_stats(_REGIONAL_CONSTRAINT_STATS).pipe(
        unify_truncating_variants
    )
    nonsense_constraint = read_nonsense_constraint(_REGIONAL_NONSENSE_CONSTRAINT).pipe(
        tidy_nonsense_constraint
    )
    omim = dac.get_omim(_GENEMAP2_SIMPLE, _GENE_IDS).drop(
        ["phenotype", "inheritance"], axis=1
    )

    return (
        dnms.merge(constraint_stats, how="right", validate="m:1")
        .fillna({"n_dnms": 0})
        .merge(nonsense_constraint, how="left", validate="m:1")
        .merge(omim, how="left", validate="m:1")
        .fillna({"inheritance_simple": "non_morbid"})
        .pipe(reorder_data)
        .pipe(write_out, _FILE_OUT)
    )


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
