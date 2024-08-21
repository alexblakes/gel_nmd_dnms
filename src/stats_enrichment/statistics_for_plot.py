"""Plot DNV enrichment in different regions and gene sets."""

import logging

import pandas as pd

import src

_PATHS = [
    "data/interim/dnms_enrichment_all_genes.tsv",
    "data/interim/dnms_enrichment_ad.tsv",
    "data/interim/dnms_enrichment_ar.tsv",
    "data/interim/dnms_enrichment_non_morbid.tsv",
]
_OUT_PATHS = [
    "data/statistics/dnms_enrichment_all_genes_tidy.tsv",
    "data/statistics/dnms_enrichment_ad_tidy.tsv",
    "data/statistics/dnms_enrichment_ar_tidy.tsv",
    "data/statistics/dnms_enrichment_non_morbid_tidy.tsv",
]


logger = logging.getLogger(__name__)


def read_data(path):
    return pd.read_csv(path, sep="\t", index_col="csq")


# def reverse_data(df):
#     return df.iloc[::-1]


def normalise_error(df):
    return df.assign(
        fc_ci_lo=lambda x: abs(x.fc_ci_lo - x.fc),
        fc_ci_hi=lambda x: x.fc_ci_hi - x.fc,
        fc_ci_null=0,
    )


def write_out(df, path):
    df.to_csv(path, sep="\t")
    return df


# def add_labels(df):
#     sigfig = lambda x: x.round(1).astype(str)

#     fc = lambda x: sigfig(x.fc)
#     ci_lo = lambda x: sigfig(x.fc - x.fc_ci_lo)
#     ci_hi = lambda x: sigfig(x.fc + x.fc_ci_hi)

#     return df.assign(labels= lambda x: fc(x) + " (" + ci_lo(x) + "-" + ci_hi(x) + ")")


# def log_transform(df):
#     return df.transform({"fc": np.log2, "fc_ci_lo": np.log2, "fc_ci_hi": np.log2})


def main():
    """Run as script."""

    for file_in, file_out in zip(_PATHS, _OUT_PATHS):
        read_data(file_in).pipe(normalise_error).pipe(write_out, file_out)


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
