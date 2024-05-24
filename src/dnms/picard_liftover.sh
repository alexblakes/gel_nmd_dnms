#!/usr/bin/env bash
set -euo pipefail
# Liftover GRCh37 VCF to GRCh38 VCF with Picard

JAR_FILE="/tools/aws-workspace-apps/re_admin/source_code/picard/3.1.1/picard.jar"
CHAIN="data/raw/grch37_to_grch38.over.chain"
BASENAME=$(basename -s .gz $1)
FASTA="/public_data_resources/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# Load module
module load picard/3.1.1

# Picard command
java -jar $JAR_FILE \
LiftoverVcf \
I=$1 \
O=$2 \
CHAIN=$CHAIN \
REJECT="data/logs/picard_liftover_rejected_${BASENAME}" \
REFERENCE_SEQUENCE=$FASTA \
WARN_ON_MISSING_CONTIG=true
