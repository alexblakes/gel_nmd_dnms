#!/usr/bin/env bash

#BSUB -q short
#BSUB -P re_gecip_enhanced_interpretation
#BSUB -J "af_gel"
#BSUB -o data/logs/cluster/af_gel_%J.out
#BSUB -e data/logs/cluster/af_gel_%J.err
#BSUB -R "rusage[mem=64000]"
#BSUB -M 64000

set -euo pipefail
source activate bio

# Get allele frequency information from GEL

DIR="/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/additional_data/allele_frequencies/"
SCRATCH="/re_scratch/re_gecip/enhanced_interpretation/ablakes"
FILE_OUT="data/interim/allele_frequencies_gel.txt.gz"

# Get GEL whole-chort AF data
find ${DIR} -name *.tsv.gz -type f \
| sort -k1,1V \
| xargs zcat \
| cut -f 1-4,28,29,30 \
| grep -v "^#" \
| sort -k1,1V -k2,2n --stable --parallel 8 --buffer-size 50% -T $SCRATCH \
| bgzip > $FILE_OUT

tabix -s1 -b2 -e2 $FILE_OUT