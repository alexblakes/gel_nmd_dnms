"""Combine DNMs after liftover."""

# Imports
from pathlib import Path

import pandas as pd

from src import constants as C
from src import setup_logger


# Module constants
_USECOLS = ["chr", "pos", "ref", "alt", "Trio Id"]
_NAMES = ["chr", "pos", "Trio Id", "ref", "alt", "qual", "filter", "info"]


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def read_vcf(path, cohort):
    """Read VCF data.

    Args:
        path (str): Path to VCF file.
        source (str): The source cohort. Either "gel" or "ddd".
    """

    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=_NAMES,
        usecols=_USECOLS,
    ).assign(**{"Trio Id": lambda x: cohort + "_" + x["Trio Id"]})

    logger.info(f"DNMs in {path}: {len(df)}")
    logger.info(f"Number of chromosomes: {df.chr.nunique()}")

    return df


def add_chr_prefix(df):
    """Add 'chr' prefix to GEL GRCh38 DNMs."""

    df["chr"] = "chr" + df["chr"].astype(str)

    return df


def concatenate_vcfs(*args):
    """Concatenate VCFs."""
    return (
        pd.concat(args)
        .assign(QUAL=".", FILTER=".", INFO=".")
        .loc[:, ["chr", "pos", "Trio Id", "ref", "alt", "QUAL", "FILTER", "INFO"]]
        .sort_values(["chr", "pos"])
    )


def main():
    """Run as script."""
    gel_38 = read_vcf(C.GEL_DNMS_VCF_38, "gel").pipe(add_chr_prefix)
    gel_lifted = read_vcf(C.GEL_37_LIFTED, "gel")
    ddd_lifted = read_vcf(C.KAPLANIS_37_LIFTED, "ddd")

    df = concatenate_vcfs(gel_38, gel_lifted, ddd_lifted)

    logger.info(f"Combined DNMs: {len(df)}")
    logger.info(f"Number of chromosomes: {df.chr.nunique()}")

    logger.info("Writing to output.")
    df.to_csv(C.DNMS_GRCH38_COMBINED, sep="\t", index=False, header=None)


if __name__ == "__main__":
    main()
