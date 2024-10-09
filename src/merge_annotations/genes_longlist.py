"""Curate non-morbid genes in OMIM with >=3 dnPTVs."""

import logging

import pandas as pd

import src

FILE_IN = "data/interim/dnms_annotated_clinical.tsv"
FILE_OUT = "data/final/candidate_genes_longlist.tsv"
USECOLS = [
    "symbol",
    "enst",
    "region",
    "constraint",
    "n_trunc_region",
    "n_obs",
    "n_exp",
    "oe",
    "oe_ci_hi",
    "loeuf",
    "pli",
    "omim_inheritance_simple",
]

logger = logging.getLogger(__name__)


def write_out(df, path=FILE_OUT):
    df.to_csv(path, sep="\t", index=False)
    return df


def main():
    """Run as script."""
    
    return (
        pd.read_csv(FILE_IN, sep="\t", usecols=USECOLS)
        .query("omim_inheritance_simple not in ['AD', 'AD_AR', 'XL']")
        .query("n_trunc_region >= 3")
        .query("constraint == 'constrained'")
        .drop_duplicates()
        .loc[:, USECOLS]
        .sort_values("symbol")
        .rename(columns={"n_trunc_region":"de_novo_ptvs"})
        .pipe(write_out)
    )


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
