#!/usr/bin/env bash
set -euo pipefail

# Annotate DNMs with VEP

DNMS="data/interim/dnms_38_combined_af.vcf.gz"
FILE_OUT="data/interim/dnms_38_combined_af_vep.vcf.gz"
VEP_CACHE="/public_data_resources/vep_resources/VEP_105"
FASTA="/public_data_resources/reference/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"

vep \
    --input_file $DNMS \
    --output_file STDOUT \
    --species homo_sapiens \
    --assembly GRCh38 \
    --dir ${VEP_CACHE} \
    --offline \
    --cache \
    --dir_cache ${VEP_CACHE} \
    --fork 4 \
    --buffer_size 10000 \
    --force_overwrite \
    --format vcf \
    --vcf \
    --fasta ${FASTA} \
    --no_stats \
    --minimal \
    --symbol \
    --canonical \
    --fields "Consequence,Feature,SYMBOL,CANONICAL" \
| filter_vep \
    --filter "CANONICAL is YES" \
    --filter "Consequence regex stop_gained|frameshift_variant" \
    --only_matched \
| bcftools +split-vep \
    --columns - \
    --duplicate \
| bcftools annotate \
    -x INFO/CSQ,INFO/CANONICAL \
    -o $FILE_OUT
