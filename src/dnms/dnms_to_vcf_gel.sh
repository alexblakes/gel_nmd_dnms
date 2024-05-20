#! usr/bin/env bash

# Split stringent GEL DNMs by assembly

FILE_IN="data/raw/denovo_flagged_variants_2023-05-16_09-51-50.tsv"
FILE_OUT_37="data/interim/dnms_37_gel.vcf"
FILE_OUT_38="data/interim/dnms_38_gel_raw.vcf"
FASTA_37="/public_data_resources/reference/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
FASTA_38="/public_data_resources/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# Convert GRCh37 variants to VCF
tail -n +2 "${FILE_IN}" \
| cut -f 1,3-7,12 \
| awk -v OFS="\t" '($2=="GRCh37" && $NF==1) {print $0}' \
| bcftools convert -c ID,-,CHROM,POS,REF,ALT,- -f "${FASTA_37}" --tsv2vcf - -Ou \
| bcftools sort -Ov -o "${FILE_OUT_37}"

# Convert GRCh38 variants to VCF
tail -n +2 "${FILE_IN}" \
| cut -f 1,3-7,12 \
| awk -v OFS="\t" '($2=="GRCh38" && $NF==1) {print $0}' \
| sed 's/\t/\tchr/2' \
| bcftools convert -c ID,-,CHROM,POS,REF,ALT,- -f "${FASTA_38}" --tsv2vcf - -Ou \
| bcftools sort -Ov -o "${FILE_OUT_38}"

# Compress and index GRCh38 variants
bgzip --keep -f "${FILE_OUT_38}"
tabix "${FILE_OUT_38}.gz"