#!/usr/bin/env bash
# Liftover GRCh37 VCF to GRCh38 VCF with Picard

# Load module
module load bio/picard/2.18.27-Java-1.8

# Picard command
java -jar $EBROOTPICARD/picard.jar \
LiftoverVcf \
I=data/interim/$1 \
O=data/interim/$2 \
CHAIN=data/raw/grch37_to_grch38.over.chain \
REJECT=data/logs/picard_liftover_rejected_$1 \
REFERENCE_SEQUENCE=/public_data_resources/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
WARN_ON_MISSING_CONTIG=true
