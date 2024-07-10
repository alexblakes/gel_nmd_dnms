"""Boilerplate code for most modules."""

import logging
from pathlib import Path

import pandas as pd

import src

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_DNMS = "data/interim/dnms_38_combined.vcf.gz"
_VCF_HEADER = "chrom pos id ref alt qual filter info".split()

logger = logging.getLogger(__name__)


def read_dnms(path=_DNMS):
    df = 
    return df

def main():
    """Run as script."""

    dnms = pd.read_csv(_DNMS, sep="\t", comment="#", header=None, names=_VCF_HEADER)
    
    return dnms


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()