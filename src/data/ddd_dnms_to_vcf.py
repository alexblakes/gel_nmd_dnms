""" Script to convert the Kaplanis DNM data to VCF format."""

# Imports
from pathlib import Path

import pandas as pd

from src import constants as C
from src import setup_logger


# Module constants
_VCF_FORMAT = "##fileformat=VCFv4.1"
_VCF_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def reformat_to_vcf(df):
    logger.info(f"Variants in raw data: {len(df)}")
    
    df = (
        df.loc[:, ["chrom", "pos", "id", "ref", "alt"]]
        .assign(QUAL=".", FILTER=".", INFO=".")
        .dropna()
    )

    logger.info(f"DNMs after dropping NaNs: {len(df)}")

    return df


def write_vcf(df, path):
    """Write to VCF with a generic header."""

    with open(path, "w") as output:
        output.write(f"{_VCF_FORMAT}\n")
        output.write(f"{_VCF_HEADER}\n")
        df.to_csv(output, sep="\t", index=False, header=None)

    return None


def main():
    """Run the script."""
    pd.read_csv(C.KAPLANIS_DNMS, sep="\t").pipe(reformat_to_vcf).pipe(write_vcf, C.KAPLANIS_DNMS_VCF)


if __name__ == "__main__":
    main()
