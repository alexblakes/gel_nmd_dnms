#!/usr/bin/env bash

#BSUB -q inter
#BSUB -P re_gecip_enhanced_interpretation
#BSUB -J "vep_dnms"
#BSUB -o data/logs/cluster/vep_dnms_%J.out
#BSUB -e data/logs/cluster/vep_dnms_%J.err
#BSUB -n 4 # needs to be the same as the number of VEP forks
#BSUB -R "rusage[mem=10000]"
#BSUB -M 16000

# Load modules
module purge
module load tools/singularity/3.8.3

# Create empty VEP input file
touch data/interim/vep_inputs.txt

# Load environment variables
source src/data/vep.conf

# Run VEP
singularity exec \
    $MOUNT_WD \
    $MOUNT_GENOMES \
    $MOUNT_GEL_DATA_RESOURCES \
    $MOUNT_PUBLIC_DATA_RESOURCES \
    $MOUNT_SCRATCH \
    $MOUNT_RE_GECIP \
    $IMG \
vep \
    --input_file data/interim/dnms_38_combined.vcf \
    --output_file data/interim/dnms_38_combined_vep.vcf \
    --species homo_sapiens \
    --assembly GRCh38 \
    --dir /opt/vep/.vep \
    --offline \
    --cache \
    --dir_cache ${VEP_CACHE} \
    --fork 4 \
    --buffer_size 100000 \
    --force_overwrite \
    --format vcf \
    --tab \
    --fasta ${REFFASTA} \
    --no_stats \
    --show_ref_allele \
    --coding_only \
    --fields "Location,REF_ALLELE,Allele,Consequence,Feature,Uploaded_variation"