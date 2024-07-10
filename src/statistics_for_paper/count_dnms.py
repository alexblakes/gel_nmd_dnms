"""Count DNMs and trios in GEL and DDD."""

import logging
from pathlib import Path

import pandas as pd

import src

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_K37 = "data/raw/kaplanis_dnms.tsv"
_K38 = "data/interim/dnms_38_kaplanis.vcf.gz"
_COMBINED38 = "data/interim/dnms_38_combined.vcf.gz"
_VCF_HEADER = "chrom pos id ref alt qual filter info".split()

logger = src.setup_logger(_LOGFILE)

k37 = pd.read_csv(_K37, sep="\t")
k38 = pd.read_csv(_K38, sep="\t", comment="#", header=None, names=_VCF_HEADER)
combined_38 = pd.read_csv(_COMBINED38, sep="\t", comment="#", header=None, names=_VCF_HEADER)

logger.info(f"Unique IDs in Kaplanis data before liftover: {k37['id'].nunique()}")
logger.info(f"DNMs in Kaplanis data before liftover: {len(k37)}")
logger.info(f"Unique DNMs in Kaplanis data before liftover (by chrom, pos, ref, alt): {len(k37.drop_duplicates(['chrom','pos','ref','alt']))}")

logger.info(f"Unique IDs in Kaplanis data after liftover: {k38['id'].nunique()}")
logger.info(f"DNMs in Kaplanis data after liftover: {len(k38)}")
logger.info(f"Unique DNMs in Kaplanis data after liftover (by chrom, pos, ref, alt): {len(k38.drop_duplicates(['chrom','pos','ref','alt']))}")

logger.info(f"DNMs in combined data: {len(combined_38)}")
logger.info(f"Unique DNMs in combined data (by chrom, pos, ref, alt): {len(combined_38.drop_duplicates(['chrom','pos','ref','alt']))}")
