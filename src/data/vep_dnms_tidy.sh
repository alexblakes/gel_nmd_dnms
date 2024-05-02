#!/usr/bin/env bash

# Tidy the VEP-annotated DNMs.
#
# Tidy such that:
# 1. The header, starting with '#', is removed
# 2. Only variants with "stop_gained", "missense", "synonymous", or "frameshift" consequences are kept
# 3. The "location" column is split into "chr" and "pos" columns
# 4. Consequences are sanitised. For example "missense_variant,splice_region_variant" 
#    becomes "missense_variant".
# 5. Only variants in the canonical transcript list are kept.
# 6. For indels, only the start position is given.
# 7. The cohort (gel, genedx, rumc, ddd) is extracted.
#
# Save the result to a TSV.

cat data/interim/dnms_38_combined_vep.vcf | \
grep -v '^#' | \
grep -f data/uom_csf/transcript_ids.tsv | \
grep -E 'stop_gained|missense|synonymous|frameshift' | \
awk -F '[:\t]' '{ $1=$1; print $0 }' | \
awk '{ sub(".*missense.*", "missense_variant", $5); print $0 }' | \
awk '{ sub(".*synonymous.*", "synonymous_variant", $5); print $0 }' | \
awk '{ sub(".*stop_gained.*", "stop_gained", $5); print $0 }' | \
awk '{ sub(".*frameshift.*", "frameshift_variant", $5); print $0 }' | \
awk '{
    split($2, a, "-"); 
    ix7 = index($7, "_"); 
    split($7, b, "_"); 
    print($1, a[1], $3, $4, $5, $6, b[1], substr($7, ix7+1))
}' | \
awk '{sub("ddd", "gdx", $7); print($0)}' | \
awk '$8 ~ /^rumc/ {$7 = "rumc"}; {print($0)}' | \
awk 'BEGIN{OFS="\t"} $8 ~ /^DDD13k/ {$7 = "ddd"} {$1=$1} {print($0)}' > \
data/interim/dnms_38_combined_vep_tidy.tsv
