#! usr/bin/env bash

# Combine GEL DNMs after liftover.

LIFTED="data/interim/dnms_38_gel_lifted.vcf.gz"
RAW="data/interim/dnms_38_gel_raw.vcf.gz"
FILE_OUT="data/interim/dnms_38_gel_combined.vcf"

bcftools concat -a "${LIFTED}" "${RAW}" -Ou \
| bcftools sort -Ov -o "${FILE_OUT}"

bgzip --keep -f "${FILE_OUT}"
tabix "${FILE_OUT}.gz"
