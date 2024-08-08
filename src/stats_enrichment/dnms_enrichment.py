"""Quantify the enrichment of truncating DNMs in constrained transcripts and regions.
"""

import logging

import numpy as np
import pandas as pd
from scipy import stats

import src
from src import stats_enrichment

_FILE_IN = "data/interim/dnms_enrichment_counts.tsv"
_DNM_ENRICHMENT_NAMES = [
    "synonymous_variant",
    "missense_variant",
    "truncating",
    "nmd_target",
    "start_proximal",
    "long_exon",
    "distal",
]
_DNM_ENRICHMENT_LABELS = [
    "Synonymous",
    "Missense",
    "Truncating (Whole transcript)",
    "Truncating (NMD target)",
    "Truncating (Start proximal)",
    "Truncating (Long exon)",
    "Truncating (Distal)",
]
_GENE_SET_NAMES = ["all_genes", "morbid_genes", "non_morbid_genes"]
_GENE_SET_LABELS = ["All genes", "Morbid genes", "Non-morbid genes"]

logger = logging.getLogger(__name__)


def read_data(path):
    return pd.read_csv(path, sep="\t")


def count_variants(df, groupby):
    return df.groupby(groupby, dropna=False).agg(
        gnomad_exp=("n_exp", "sum"), dnms_obs=("n_dnms", "sum")
    )


def group_levels(obj, drop="constraint"):
    """Group a dataframe by all levels of a multi-index, except one."""

    return obj.groupby(level=[i for i in obj.index.names if i != drop])


def chi_square_test(df):
    return group_levels(df).apply(
        lambda x: stats.chisquare(x.dnms_obs, x.dnms_exp).pvalue
    )


def get_statistics(df):
    return df.assign(
        prop_exp=lambda x: x.gnomad_exp / group_levels(x)["gnomad_exp"].sum(),
        dnms_exp=lambda x: x.prop_exp * group_levels(x)["dnms_obs"].sum(),
        oe=lambda x: x.dnms_obs / x.dnms_exp,
        fc=lambda x: x.oe / x.xs("unconstrained", level="constraint").oe,
        log2_fc=lambda x: np.log2(x.fc),
        oe_ci=lambda x: (1.96 * np.sqrt(x.dnms_obs)) / x.dnms_exp,
        fc_ci=lambda x: x.oe_ci / x.xs("unconstrained", level="constraint").oe,
        chi2_p=lambda x: chi_square_test(x),
    )


def tidy_data(df):
    return df.xs("constrained", level="constraint").loc[:, ["fc", "fc_ci", "chi2_p"]]


def get_csqs_and_regions(df):
    slices = [
        ("synonymous_variant", "transcript"),
        ("missense_variant", "transcript"),
        ("truncating", "transcript"),
        ("truncating", "nmd_target"),
        ("truncating", "start_proximal"),
        ("truncating", "long_exon"),
        ("truncating", "distal_nmd"),
    ]

    return df.loc[slices, :]


def write_out(df, path):
    df.to_csv(path, sep="\t", index=False)
    return df


def main():
    """Run as script."""

    df = read_data(_FILE_IN)

    genes_by_omim_status = (
        count_variants(df, ["csq", "region", "inheritance_simple", "constraint"])
        .pipe(get_statistics)
        .pipe(tidy_data)
    )

    ad_genes = genes_by_omim_status.xs("AD", level="inheritance_simple").pipe(
        get_csqs_and_regions
    )
    ar_genes = genes_by_omim_status.xs("AR", level="inheritance_simple").pipe(
        get_csqs_and_regions
    )
    non_morbid_genes = genes_by_omim_status.xs(
        "non_morbid", level="inheritance_simple"
    ).pipe(get_csqs_and_regions)

    all_genes = count_variants(df, ["csq","region","constraint"]).pipe(get_statistics).pipe(tidy_data).pipe(get_csqs_and_regions)

    return all_genes

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
    logger = src.setup_logger(src.log_file(__file__))
    main()
