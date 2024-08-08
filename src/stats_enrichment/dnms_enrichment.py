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
    "Synonymous",
    "Missense",
    "PTV (Whole transcript)",
    "PTV (NMD target)",
    "PTV (Start proximal)",
    "PTV (Long exon)",
    "PTV (Distal)",
]
_FILE_OUT = "data/statistics/dnms_enrichment.tsv"

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


def get_statistics(df):
    return df.assign(
        prop_exp=lambda x: x.gnomad_exp / group_levels(x)["gnomad_exp"].sum(),
        dnms_exp=lambda x: x.prop_exp * group_levels(x)["dnms_obs"].sum(),
        oe=lambda x: x.dnms_obs / x.dnms_exp,
        fc=lambda x: x.oe / x.xs("unconstrained", level="constraint").oe,
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
    all_genes = get_dnms_enrichment(df, group_all_genes, bootstrap_samples=1000).pipe(
        tidy_index
    )

    # Get statistics for genes per OMIM category
    group_omim_genes = ["csq", "region", "inheritance_simple", "constraint"]
    omim_genes = get_dnms_enrichment(df, group_omim_genes, bootstrap_samples=1000)
    moi = ["AD", "AR", "non_morbid"]
    ad, ar, non_morbid = [separate_omim_categories(omim_genes, m) for m in moi]

    # Write to output
    gene_sets = [all_genes, ad, ar, non_morbid]
    labels = ["_all_genes", "_ad", "_ar", "_non_morbid"]

    for df, label in zip(gene_sets, labels):
        write_out(df, get_path(label))

    return gene_sets


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
