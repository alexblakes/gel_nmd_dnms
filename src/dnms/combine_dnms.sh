#! usr/bin/env bash

# Combine GEL and Kaplanis DNMs after liftover.

GEL_LIFTED="data/interim/dnms_38_gel_lifted.vcf.gz"
GEL_RAW="data/interim/dnms_38_gel_raw.vcf.gz"
KAPLANIS_LIFTED="data/interim/dnms_38_kaplanis.vcf.gz"
FILE_OUT="data/interim/dnms_38_combined.vcf.gz"

bcftools concat -a -Ou $GEL_LIFTED $GEL_RAW $KAPLANIS_LIFTED \
| bcftools sort --write-index=tbi -o "${FILE_OUT}"
