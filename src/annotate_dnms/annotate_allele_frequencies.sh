#!/usr/bin/env bash
set -euo pipefail

# Annotate DNMs with GEL and gnomAD allele frequencies

DNMS="data/interim/dnms_38_combined.vcf.gz"
GEL_AFS="data/interim/allele_frequencies_gel.txt.gz"
GNOMAD_AFS="data/interim/allele_frequencies_gnomad_v3.1.1_genomes.vcf.gz"
HEADER="data/manual/header_allele_frequencies_gel.txt"
FILE_OUT="data/interim/dnms_38_combined_af.vcf.gz"

# Write header
echo '##INFO=<ID=ac_gel,Number=1,Type=Integer,Description="Allele count in GEL genomes">' > $HEADER
echo '##INFO=<ID=an_gel,Number=1,Type=Integer,Description="Allele number in GEL genomes">' >> $HEADER
echo '##INFO=<ID=af_gel,Number=1,Type=Float,Description="Allele frequency in GEL genomes">' >> $HEADER

bcftools annotate \
    -a $GNOMAD_AFS \
    -c ac_gnomad:=AC,an_gnomad:=AN,af_gnomad:=AF \
    -Ou \
    $DNMS \
| bcftools annotate \
    -a $GEL_AFS \
    -c CHROM,POS,REF,ALT,ac_gel,an_gel,af_gel \
    -h $HEADER \
    -W=tbi \
    -o $FILE_OUT
