"""Fill gaps in REF/ALT columns in Kaplanis DNM data."""

import logging
from pathlib import Path

import pandas as pd
import pysam

import src


_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FILE_IN = "data/raw/kaplanis_dnms.tsv"
_FILE_OUT = "data/interim/kaplanis_dnms_filled.tsv"
_FASTA = "/public_data_resources/reference/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"

logger = logging.getLogger(__name__)


def read_dnms(path):
    df = pd.read_csv(
        path,
        sep="\t",
        usecols=["chrom", "pos", "id", "ref", "alt", "study"],
    ).fillna("")

    logger.info(f"Gaps in REF: {(df.ref=='').sum()}")
    logger.info(f"Gaps in ALT: {(df.alt=='').sum()}")

    return df


def fill_gaps(row, fasta):
    """Fill gaps in REF and ALT alleles."""

    if (row["ref"] != "") & (row["alt"] != ""):
        pass

    if row["ref"] == "":
        pos = row["pos"] - 1
        region = f"{row.chrom}:{pos}-{pos}"
        ref = fasta.fetch(region=region)

        row["pos"] = pos
        row["ref"] = ref
        row["alt"] = ref + row["alt"]

    if row["alt"] == "":
        pos = row["pos"] - 1
        region = f"{row.chrom}:{pos}-{pos}"
        ref = fasta.fetch(region=region)

        row["pos"] = pos
        row["ref"] = ref + row["ref"]
        row["alt"] = ref

    return row


def main():
    (
        read_dnms(_FILE_IN)
        .apply(fill_gaps, fasta=pysam.FastaFile(_FASTA), axis=1)
        .to_csv(_FILE_OUT, sep="\t", header=False, index=False)
    )
    logger.info("Gaps filled.")


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
