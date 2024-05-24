#!/usr/bin/env bash

#BSUB -q short
#BSUB -P re_gecip_enhanced_interpretation
#BSUB -J "af_gnomad"
#BSUB -o data/logs/cluster/af_gnomad_%J.out
#BSUB -e data/logs/cluster/af_gnomad_%J.err
#BSUB -R "rusage[mem=64000]"
#BSUB -M 64000

set -euo pipefail
source activate bio

# Get allele frequency data from gnomAD

DIR_GNOMAD_VCFS="/public_data_resources/gnomad/v3.1.1/vcf/"
FILE_GNOMAD_PATHS="data/scratch/gnomad_file_paths.txt"
FILE_OUT="data/interim/allele_frequencies_gnomad_v3.1.1_genomes.vcf.gz"
FILE_INTERIM_VCF_PATHS="data/scratch/gnomad_interim_vcf_paths.txt"
export FILE_DNMS="data/interim/dnms_38_combined.vcf.gz"

extract_afs() {
    # Extract AF data from each gnomAD VCF chunk
    bcftools view \
        -R $FILE_DNMS \
        -T $FILE_DNMS \
        -e 'FILTER ~ "AC0"' \
        -Ou \
        $1 \
    | bcftools annotate \
        -x ^INFO/AC,INFO/AN,INFO/AF \
        -W=tbi \
        -o $2
}

export -f extract_afs

# Create a text file with paths to the gnomAD VCF data
find ${DIR_GNOMAD_VCFS} -name *.vcf.bgz -type f \
| sort -k1,1V \
> $FILE_GNOMAD_PATHS

# Create a text file with paths to the interim processed VCF chunks
xargs -a $FILE_GNOMAD_PATHS -n1 basename \
| sed 's/^/data\/scratch\//' \
> $FILE_INTERIM_VCF_PATHS

# Extract AF data from each chunk
parallel \
    --link \
    --arg-file $FILE_GNOMAD_PATHS \
    --arg-file $FILE_INTERIM_VCF_PATHS \
    --keep-order \
    extract_afs

# Concatenate the interim VCFs
bcftools concat \
    --allow-overlaps \
    --file-list $FILE_INTERIM_VCF_PATHS \
    -Ou \
| bcftools sort \
    -m 4G \
    -W=tbi \
    -o $FILE_OUT

# Clean up
rm $FILE_GNOMAD_PATHS

while read -r LINE; do 
    rm $LINE "${LINE}.tbi"
done < $FILE_INTERIM_VCF_PATHS

rm $FILE_INTERIM_VCF_PATHS