#!/usr/bin/env bash
set -euo pipefail

# Count the number of GEL DNMs
GEL_LIFTED="data/interim/dnms_38_gel_lifted.vcf.gz"
GEL_RAW="data/interim/dnms_38_gel_raw.vcf.gz"
KAPLANIS_LIFTED="data/interim/dnms_38_kaplanis.vcf.gz"

WC_37=$(zcat $GEL_LIFTED | grep -v "^#" | wc -l)
WC_38=$(zcat $GEL_RAW | grep -v "^#" | wc -l)
WC_KAPLANIS=$(zcat $KAPLANIS_LIFTED | grep -v "^#" | wc -l)

echo "GEL GRCh37 variants after liftOver:" $WC_37
echo "GEL GRCh38 variants:" $WC_38
echo "GEL total DNMs:" $(($WC_37 + $WC_38))
echo "Kaplanis DNMs after liftover:" $WC_KAPLANIS
