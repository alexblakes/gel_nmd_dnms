"""Quantify the enrichment of truncating DNMs in constrained transcripts and regions.
"""

import logging
from pathlib import Path

import pandas as pd

import src

_FILE_IN = "data/interim/dnms_enrichment_counts.tsv"
_SLICES = [
    ("synonymous_variant", "transcript"),
    ("missense_variant", "transcript"),
    ("truncating", "transcript"),
    ("truncating", "nmd_target"),
    ("truncating", "start_proximal"),
    ("truncating", "long_exon"),
    ("truncating", "distal_nmd"),
]
_DNM_ENRICHMENT_LABELS = [
    "Synonymous (Full CDS)",
    "Missense (Full CDS)",
    "Nonsense & FS (Full CDS)",
    "Nonsense & FS (NMD target)",
    "Nonsense & FS (Start proximal)",
    "Nonsense & FS (Long exon)",
    "Nonsense & FS (Distal)",
]
_FILE_OUT = "data/interim/dnms_enrichment.tsv"

logger = logging.getLogger(__name__)


def read_data(path):
    return pd.read_csv(path, sep="\t")


def count_variants(df, grouping):
    return df.groupby(grouping).agg(
        gnomad_exp=("n_exp", "sum"), dnms_obs=("n_dnms", "sum")
    )


def group_levels(obj, drop="constraint"):
    """Group a dataframe by all levels of a multi-index, except one."""

    return obj.groupby(level=[i for i in obj.index.names if i != drop])


<<<<<<< HEAD
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
    total_dnms = group_levels(observed_dnms).transform("sum") #! DEBUG
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
=======
def get_statistics(df):
    return df.assign(
        prop_exp=lambda x: x.gnomad_exp / group_levels(x)["gnomad_exp"].sum(),
        dnms_exp=lambda x: x.prop_exp * group_levels(x)["dnms_obs"].sum(),
        oe=lambda x: x.dnms_obs / x.dnms_exp,
        fc=lambda x: x.oe / x.xs("unconstrained", level="constraint").oe,
>>>>>>> enrichment
    )


def get_fc(df):
    return df.xs("constrained", level="constraint").loc[:, ["fc"]]


def resample(df, grouping):
    return df.groupby(grouping).sample(frac=1, replace=True)


def bootstrap_fc_confint(df, grouping, resamples):
    return pd.concat(
        [
            resample(df, grouping)
            .pipe(count_variants, grouping)
            .pipe(get_statistics)
            .pipe(get_fc)
            for i in range(resamples)
        ],
        axis=1,
    )


def get_confints(df, q_lo=0.025, q_hi=0.975):
    return df.quantile([q_lo, q_hi], axis=1).T.set_axis(
        ["fc_ci_lo", "fc_ci_hi"], axis=1
    )


def get_p_vals(df, h0=1):
    less = df.lt(h0, axis=0).sum(axis=1) / len(df.columns)
    greater = df.gt(h0, axis=0).sum(axis=1) / len(df.columns)

    return (pd.concat([less, greater], axis=1).min(axis=1) * 2).rename("p")


def bootstrap_statistics(df):
    return pd.concat([get_confints(df), get_p_vals(df)], axis=1)


def get_dnms_enrichment(df, grouping, bootstrap_samples=100):

    statistics = count_variants(df, grouping).pipe(get_statistics).pipe(get_fc)

    bootstrap_stats = bootstrap_fc_confint(df, grouping, bootstrap_samples).pipe(
        bootstrap_statistics
    )

    return pd.concat([statistics, bootstrap_stats], axis=1)


def tidy_index(df):
    return df.loc[_SLICES, :].set_axis(_DNM_ENRICHMENT_LABELS).rename_axis("csq")


def separate_omim_categories(df, moi):
    return df.xs(moi, level="inheritance_simple").pipe(tidy_index)


def get_path(label, path=_FILE_OUT):
    path = Path(path)
    return path.with_name(path.stem + label + path.suffix)


def write_out(df, path):
    df.to_csv(path, sep="\t")
    return df


def main():
    """Run as script."""

    df = read_data(_FILE_IN).dropna()

    # Get statistics for all genes
    group_all_genes = ["csq", "region", "constraint"]
    all_genes = get_dnms_enrichment(df, group_all_genes, bootstrap_samples=10000).pipe(
        tidy_index
    )

<<<<<<< HEAD
    # Merge constraint and OMIM annotations
    constraint = constraint_stats.merge(
        nonsense_constraint, how="inner", validate="many_to_one"
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

    constraint_transcripts = constraint[m0].copy()
    constraint_transcripts_non_morbid = constraint[m0 & m1].copy()

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
=======
    # Get statistics for genes per OMIM category
    group_omim_genes = ["csq", "region", "inheritance_simple", "constraint"]
    omim_genes = get_dnms_enrichment(df, group_omim_genes, bootstrap_samples=10000)
    moi = ["AD", "AR", "non_morbid"]
    ad, ar, non_morbid = [separate_omim_categories(omim_genes, m) for m in moi]
>>>>>>> enrichment

    # Write to output
    gene_sets = [all_genes, ad, ar, non_morbid]
    labels = ["_all_genes", "_ad", "_ar", "_non_morbid"]

    for df, label in zip(gene_sets, labels):
        write_out(df, get_path(label))

    return gene_sets


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
