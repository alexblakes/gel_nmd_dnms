"""Quantify the enrichment of truncating DNMs in constrained transcripts and regions.
"""

# Imports
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chisquare
from scipy.stats import poisson

import src
from src.merge_annotations import dnms_annotate_constraint
from src import stats_enrichment

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_DNMS_ANNOTATED = "data/interim/dnms_annotated.tsv"
_REGIONAL_NONSENSE_CONSTRAINT = "data/uom_csf/regional_nonsense_constraint.tsv"
_REGIONAL_CONSTRAINT_STATS = "data/uom_csf/regional_constraint_stats.tsv"
_GENE_IDS = "data/uom_csf/gene_ids.tsv"
_GENEMAP2_SIMPLE = "data/uom_csf/genemap2_simple.tsv"
_DNM_ENRICHMENT_CATEGORIES = [
    "synonymous_variant",
    "missense_variant",
    "truncating",
    "nmd_target",
    "start_proximal",
    "long_exon",
    "distal",
][
    ::-1
]  # Reversed
_DNM_ENRICHMENT_LABELS = [
    "Synonymous",
    "Missense",
    "Truncating (Whole transcript)",
    "Truncating (NMD target)",
    "Truncating (Start proximal)",
    "Truncating (Long exon)",
    "Truncating (Distal)",
][
    ::-1
]  # Reversed
_GENE_SET_NAMES = ["all_genes", "morbid_genes", "non_morbid_genes"]
_GENE_SET_LABELS = ["All genes", "Morbid genes", "Non-morbid genes"]

logger = logging.getLogger(__name__)


def read_dnms(path):
    return pd.read_csv(
        path, sep="\t", usecols=["csq", "constraint", "inheritance_simple", "region"]
    )


def unify_truncating_variants(df):
    """Convert all PTV annotations to "truncating"."""

    df.loc[
        df["csq"].isin(["stop_gained", "frameshift", "frameshift_variant"]), "csq"
    ] = "truncating"

    return df


def read_nonsense_constraint(path):
    return pd.read_csv(path, sep="\t", usecols=["enst", "region", "constraint"])


def read_constraint_stats(path):
    return pd.read_csv(path, sep="\t", usecols=["enst", "region", "csq", "n_exp"])


def group_levels(obj):
    """Group a dataframe by levels of a multi-index. (All but the last level)."""

    n_levels = len(obj.index.levels)

    return obj.groupby(level=list(range(n_levels - 1)))


def get_expected_proportions(number_expected):
    """Find the proportion of variants expected in constrained / unconstrained regions."""

    # Absolute number of variants in constrained and unconstrained regions
    total = group_levels(number_expected).transform("sum")

    # Proportion of variants in constrained / unconstrained regions
    proportion_expected = (number_expected / total).rename("proportion_expected")

    assert np.allclose(
        group_levels(proportion_expected).sum(), 1
    ), "Proportions do not sum to 1."

    return proportion_expected


def get_expected_dnms(observed_dnms, number_expected):
    """Get the expected number of DNMs per group."""

    observed_dnms = observed_dnms.rename("dnms_observed")
    total_dnms = group_levels(observed_dnms).transform("sum")
    proportion_expected = get_expected_proportions(number_expected)
    expected_dnms = (total_dnms * proportion_expected).rename("dnms_expected")

    df = pd.concat([observed_dnms, proportion_expected, expected_dnms], axis=1)

    return df


def chi_square_test(df):
    """Get chi-square statistics for DNMs in constrained and unconstrained regions."""

    chi2 = group_levels(df).apply(
        lambda x: chisquare(x["dnms_observed"], x["dnms_expected"])
    )
    chi2 = pd.DataFrame(chi2.to_list(), columns=["chi2", "p"], index=chi2.index)

    return chi2


def relative_enrichment(df):
    oe = df["dnms_observed"] / df["dnms_expected"]

    return (oe.xs("constrained", level=-1) / oe.xs("unconstrained", level=-1)).rename(
        "relative_enrichment"
    )


def poisson_ci(df):
    # Get Poisson confidence interval (normal approximation)
    ci = poisson(df["dnms_observed"]).interval(0.95)
    ci = (
        pd.DataFrame.from_dict(ci)
        .T.set_index(df.index)
        .set_axis(["ci_l", "ci_r"], axis=1)
    )

    # Adjust for expected DNMs
    ci_adj = ci.div(df["dnms_expected"], axis=0)

    # Normalise relative to unconstrained regions
    ci_adj_constrained = ci_adj.xs("constrained", level=-1)

    oe = df["dnms_observed"] / df["dnms_expected"]
    oe_unconstrained = oe.xs("unconstrained", level=-1)

    ci_norm = ci_adj_constrained.div(oe_unconstrained, axis=0)

    return ci_norm


def unify_index(df, names=["csq"]):
    """Rename an index, if the number of index levels matches the number of names."""
    if df.index.nlevels == len(names):
        df.index = df.index.set_names(names)

    return df


def get_enrichment_statistics(df):
    chi2 = chi_square_test(df)
    enrichment = relative_enrichment(df)
    ci = poisson_ci(df)

    return pd.concat([chi2, enrichment, ci], axis=1).pipe(unify_index)


def tidy_data(df):
    df = df.reset_index()

    # Drop the needless "transcript" rows
    df = df[df["csq"] != "transcript"]

    # Rename the "distal_nmd" consequence label
    df = df.replace({"distal_nmd": "distal"})

    # Rename the "AD_AR" inheritance label
    df = df.replace({"AD_AR": "AD/AR"})

    # Capitalise constraint annotations (if they are present)
    try:
        df["constraint"] = df["constraint"].str.capitalize()
    except:
        pass

    # Convert consequences to a categorical index
    df = stats_enrichment.sort_region_column(
        df,
        column="csq",
        categories=_DNM_ENRICHMENT_CATEGORIES,
        labels=_DNM_ENRICHMENT_LABELS,
    )

    logger.debug(f"{df}")

    return df


def main():
    """Run as script."""

    # Load datasets
    dnms = read_dnms(_DNMS_ANNOTATED).pipe(unify_truncating_variants)
    nonsense_constraint = read_nonsense_constraint(_REGIONAL_NONSENSE_CONSTRAINT)
    constraint_stats = read_constraint_stats(_REGIONAL_CONSTRAINT_STATS).pipe(
        unify_truncating_variants
    )
    omim = dnms_annotate_constraint.get_omim(_GENEMAP2_SIMPLE, _GENE_IDS).drop(
        ["phenotype", "inheritance"], axis=1
    )

    # Merge constraint and OMIM annotations
    constraint = constraint_stats.merge(
        nonsense_constraint, validate="many_to_one"
    ).merge(omim, how="left", validate="many_to_one")

    """ 
    We want to find the number of DNMs (in our trio cohorts), and the expected
    number of rare SNVs (in gnomAD), in six subgroups:
    
     - Synonymous / missense / truncating variants in the whole transcript for:
       - All genes
       - OMIM morbid genes by inheritance pattern
       - Non-morbid genes
    
    - Truncating variants in NMD regions, for each of the above gene sets.

    The statistics in each subgroup will be stratified by constraint.
    """

    # Get subgroups of constrained regions
    m0 = constraint.region == "transcript"
    m1 = constraint.inheritance_simple.isna()
    m2 = constraint.csq == "truncating"

    constraint_transcripts = constraint[m0]
    constraint_transcripts_non_morbid = constraint[m0 & m1]

    constraint_truncating = constraint[m2].copy()
    constraint_truncating_non_morbid = constraint[m2 & m1].copy()

    # Get subgroups of DNMs
    m4 = dnms.csq == "truncating"
    m5 = dnms.inheritance_simple.isna()

    dnms_truncating = dnms[m4].copy()
    dnms_non_morbid = dnms[m5].copy()
    dnms_truncating_non_morbid = dnms[m4 & m5].copy()

    # 1. Transcript-level enrichment of synonymous / missense / truncating variants in all genes.
    group1 = lambda x: x.groupby(["csq", "constraint"])
    exp_transcript_all_genes = group1(constraint_transcripts)["n_exp"].sum()
    dnm_transcript_all_genes = group1(dnms)["csq"].count()

    # 2. Transcript-level enrichment of synonymous / missense / truncating variants in morbid genes
    group2 = lambda x: x.groupby(["csq", "inheritance_simple", "constraint"])
    exp_transcript_morbid = group2(constraint_transcripts)["n_exp"].sum()
    dnm_transcript_morbid = group2(dnms)["csq"].count()

    # 3. Transcript-level enrichment of synonymous / missense / truncating variants in non-morbid genes
    exp_transcript_non_morbid = group1(constraint_transcripts_non_morbid)["n_exp"].sum()
    dnm_transcript_non_morbid = group1(dnms_non_morbid)["csq"].count()

    # 4. Regional enrichment of truncating variants in all genes
    group3 = lambda x: x.groupby(["region", "constraint"])
    exp_regions_all_genes = group3(constraint_truncating)["n_exp"].sum()
    dnm_regions_all_genes = group3(dnms_truncating)["csq"].count()

    # 5. Regional enrichment of truncating variants in morbid genes
    group4 = lambda x: x.groupby(["region", "inheritance_simple", "constraint"])
    exp_regions_morbid = group4(constraint_truncating)["n_exp"].sum()
    dnm_regions_morbid = group4(dnms_truncating)["csq"].count()

    # 6. Regional enrichment of truncating variants in non-morbid genes
    exp_regions_non_morbid = group3(constraint_truncating_non_morbid)["n_exp"].sum()
    dnm_regions_non_morbid = group3(dnms_truncating_non_morbid)["csq"].count()

    # Organise into lists of dataframes, which can be zipped
    dnms_counts_transcripts = [
        dnm_transcript_all_genes,
        dnm_transcript_morbid,
        dnm_transcript_non_morbid,
    ]
    dnms_counts_regions = [
        dnm_regions_all_genes,
        dnm_regions_morbid,
        dnm_regions_non_morbid,
    ]
    expected_totals_transcripts = [
        exp_transcript_all_genes,
        exp_transcript_morbid,
        exp_transcript_non_morbid,
    ]
    expected_totals_regions = [
        exp_regions_all_genes,
        exp_regions_morbid,
        exp_regions_non_morbid,
    ]

    # Adjust DNM count by the number of SNVs expected
    exp_dnms_transcripts = [
        get_expected_dnms(dnms, exp)
        for dnms, exp in zip(dnms_counts_transcripts, expected_totals_transcripts)
    ]
    exp_dnms_regions = [
        get_expected_dnms(dnms, exp)
        for dnms, exp in zip(dnms_counts_regions, expected_totals_regions)
    ]

    # Get enrichment statistics
    stats_transcripts = [get_enrichment_statistics(x) for x in exp_dnms_transcripts]
    stats_regions = [get_enrichment_statistics(x) for x in exp_dnms_regions]

    # Organise the data into gene sets (all genes, morbid genes, non-morbid genes)
    counts = [pd.concat([a, b]) for a, b in zip(exp_dnms_transcripts, exp_dnms_regions)]
    stats = [pd.concat([a, b]) for a, b in zip(stats_transcripts, stats_regions)]

    # Tidy the data
    counts = [tidy_data(df) for df in counts]
    stats = [tidy_data(df) for df in stats]

    # Record count data in logs
    for df, label in zip(counts, _GENE_SET_LABELS):
        logger.info(f"DNM counts in {label}:\n{df}")

    # Write to output
    for df, name in zip(stats, _GENE_SET_NAMES):
        df.to_csv(f"data/statistics/dnms_enrichment_{name}.tsv", sep="\t", index=False)

    return stats


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
