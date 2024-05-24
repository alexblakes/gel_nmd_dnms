#!/usr/bin/env bash
set -euo pipefail

# Filter by allele count, tidy VEP annotations, tidy cohort and ID.

FILE_IN="data/interim/dnms_38_combined_af_vep.vcf.gz"
FILE_OUT="data/interim/dnms_38_combined_af_vep_tidy.tsv"

bcftools view -i 'ac_gnomad <= 1 | ac_gnomad = "."' $FILE_IN \
|bcftools view -i 'af_gnomad <= 0.0001 | af_gnomad = "."' \
|bcftools view -i 'ac_gel <= 5 | ac_gel = "."' \
|bcftools annotate -x INFO/ac_gnomad,INFO/af_gnomad,INFO/an_gnomad,INFO/ac_gel,INFO/an_gel,INFO/af_gel \
|bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%Consequence\t%Feature\t%ID' \
| awk '
        $5~"stop_gained" {$5="nonsense"};
        $5~"frameshift_variant" {$5="frameshift"};
        1
' \
| awk -v OFS="\t" '
        $7 ~ /^GDX/ {$8 = "GDX"; split($7, a, "_"); $9 = a[2]};
        $7 ~ /^DDD/ {$8 = "DDD"; split($7, a, /_|\./); $9 = a[length(a)]};
        $7 ~ /^RUMC/ {$8 = "RUMC"; split($7, a, "_"); $9 = a[length(a)]};
        $7 !~ /^[GDX|DDD|RUMC]/ {$8 = "GEL"; split($7, a, "_"); $9 = a[1]};
        {$1 = $1}; 1
' \
| cut -f 1-6,8- \
> $FILE_OUT
