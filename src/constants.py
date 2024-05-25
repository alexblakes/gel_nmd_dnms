# Project-wide constants

# Visualisation
## Plotting
CM = 1 / 2.54  # cm to inches conversion
STYLE_DEFAULT = "src/visualisation/styles/default.mplstyle"
COLOR_VIBRANT = "src/visualisation/styles/color/vibrant.mplstyle"
COLOR_REGIONS = "src/visualisation/styles/color/regions.mplstyle"
COLOR_ENRICHMENT = "src/visualisation/styles/color/enrichment.mplstyle"

## Labels
REGIONS = "whole_transcript nmd_target start_proximal long_exon distal".split()
REGION_LABELS = [
    "Whole transcript",
    "NMD target",
    "Start proximal",
    "Long exon",
    "Distal",
]

# Files
## Data
### Raw
# GEL_DNMS = "data/raw/denovo_flagged_variants_2023-05-16_09-51-50.tsv"
# GEL_DE_NOVO_COHORT = "data/raw/denovo_cohort_information_2023-05-16_10-07-17.tsv"
# KAPLANIS_DNMS = "data/raw/kaplanis_dnms.tsv"
# KAPLANIS_DNMS_VCF = "data/interim/variants_37_kaplanis_dnms.vcf"
# LABKEY_EXIT_QUESTIONNAIRE = "data/raw/gmc_exit_questionnaire_2023-08-01_09-42-17.tsv"
# LABKEY_PARTICIPANTS = "data/raw/participant_2023-08-09_12-18-36.tsv"

### Interim
# DNMS_GRCH38_COMBINED = "data/interim/dnms_38_combined.vcf"
# GEL_DNMS_VCF_37 = "data/interim/variants_37_gel_stringent.vcf"
# GEL_DNMS_VCF_38 = "data/interim/variants_38_gel_stringent.vcf"
# GEL_37_LIFTED = "data/interim/liftover_38_gel.vcf"
# KAPLANIS_37_LIFTED = "data/interim/liftover_38_ddd.vcf"
# LABKEY_CLINICAL = "data/interim/labkey_participant_clinical.tsv"
# LABKEY_EQ_ACMG = "data/interim/labkey_exit_questionnaire_acmg.tsv"

### From UoM CSF
# GENE_IDS = "data/uom_csf/gene_ids.tsv"
# GENEMAP2_SIMPLE = "data/uom_csf/genemap2_simple.tsv"
# NMD_ANNOTATION = "data/uom_csf/nmd_annotations.tsv"
# REGIONAL_NONSENSE_CONSTRAINT = "data/uom_csf/regional_nonsense_constraint.tsv"
# REGIONAL_CONSTRAINT_STATS = "data/uom_csf/regional_constraint_stats.tsv"