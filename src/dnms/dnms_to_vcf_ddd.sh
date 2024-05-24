#! /usr/bin/env bash
# Convert Kaplanis et al. DNMs to VCF

FASTA="/public_data_resources/reference/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
FILE_IN="data/interim/kaplanis_dnms_filled.tsv"
FILE_OUT="data/interim/dnms_37_kaplanis.vcf.gz"

awk -v OFS="\t" '{print $6"_"$1, $2, $3, $4, $5}' $FILE_IN \
| bcftools convert -c ID,CHROM,POS,REF,ALT -f "${FASTA}" --tsv2vcf - -Ou \
| bcftools sort --write-index=tbi -o $FILE_OUT
