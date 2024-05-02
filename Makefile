SHELL = bash

.ONESHELL :

.PHONY : all

all : data/interim/variants_37_kaplanis_dnms.vcf \
      data/interim/liftover_38_ddd.vcf \
      data/interim/gel_dnm_offspring_clean.tsv \
	  data/interim/variants_37_gel_stringent.vcf \
	  data/interim/liftover_38_gel.vcf \
	  data/interim/dnms_38_combined.vcf \
	  data/interim/dnms_38_combined_vep.vcf \
	  data/interim/dnms_38_combined_vep_tidy.tsv \
	  data/interim/labkey_exit_questionnaires_case_solved.tsv \
	  data/interim/labkey_exit_questionnaire_acmg.tsv \
	  data/interim/labkey_phenotypes_clean.tsv \
	  data/interim/labkey_tiering_clean.tsv \
	  data/interim/labkey_tiering_highest_tiers.tsv \
	  data/interim/labkey_participant_clinical.tsv \
	  data/interim/dnms_annotated.tsv \
	  data/interim/dnms_annotated_clinical.tsv \
	  data/statistics/case_solved_odds_ratios.tsv \

# Convert Kaplanis DNMs to VCF
data/interim/variants_37_kaplanis_dnms.vcf : src/data/ddd_dnms_to_vcf.py
	python3 -m src.data.ddd_dnms_to_vcf

# LiftOver Kaplanis DNMs
data/interim/liftover_38_ddd.vcf : data/interim/variants_37_kaplanis_dnms.vcf \
                                   src/data/picard_liftover_submit.sh
	bash src/data/picard_liftover_submit.sh variants_37_kaplanis_dnms.vcf liftover_38_ddd.vcf

# Clean GEL DNM cohort
data/interim/gel_dnm_offspring_clean.tsv : src/data/gel_dnms_offspring.py
	python3 -m src.data.gel_dnms_offspring

# Convert GEL DNMs (both GRCh37 and GRCh38) to VCF
data/interim/variants_37_gel_stringent.vcf : data/interim/gel_dnm_offspring_clean.tsv \
                                             src/data/gel_dnms_variants.py
	python3 -m src.data.gel_dnms_variants

# LiftOver GEL DNMs
data/interim/liftover_38_gel.vcf : data/interim/variants_37_gel_stringent.vcf
	bash src/data/picard_liftover_submit.sh variants_37_gel_stringent.vcf liftover_38_gel.vcf

# Combine DNMs after liftover
data/interim/dnms_38_combined.vcf : data/interim/variants_38_gel_stringent.vcf \
                                    data/interim/liftover_38_gel.vcf \
									data/interim/liftover_38_ddd.vcf \
									src/data/dnms_combine_after_liftover.py
	python3 -m src.data.dnms_combine_after_liftover

# Annotated DNMs with VEP
data/interim/dnms_38_combined_vep.vcf : data/interim/dnms_38_combined.vcf \
                                        src/data/vep.sh
	bsub < src/data/vep.sh
	bwait -w "done(vep_dnms)"

# Tidy DNMs after VEP
data/interim/dnms_38_combined_vep_tidy.tsv : src/data/vep_dnms_tidy.sh \
                                             data/interim/dnms_38_combined_vep.vcf \
											 src/data/dnms_tidy_logging.py
	bash src/data/vep_dnms_tidy.sh
	python3 -m src.data.dnms_tidy_logging

# Tidy LabKey exit questionnaire data
data/interim/labkey_exit_questionnaires_case_solved.tsv \
data/interim/labkey_exit_questionnaire_acmg.tsv : data/raw/gmc_exit_questionnaire_2023-08-01_09-42-17.tsv \
                                                  src/data/labkey_tidy_exit_questionnaires.py
	python3 -m src.data.labkey_tidy_exit_questionnaires

# Tidy labkey phenotype data
data/interim/labkey_phenotypes_clean.tsv : src/data/labkey_tidy_phenotypes.py \
                                           data/raw/rare_diseases_participant_phen_2023-08-01_09-40-58.tsv
	python3 -m src.data.labkey_tidy_phenotypes


# Tidy tiering data
data/interim/labkey_tiering_clean.tsv \
data/interim/labkey_tiering_highest_tiers.tsv : src/data/labkey_tidy_tiering.py \
                                                data/raw/tiering_data_2023-08-01_09-44-03.tsv
	python3 -m src.data.labkey_tidy_tiering

# Merge clinical annotations
data/interim/labkey_participant_clinical.tsv : data/interim/labkey_exit_questionnaires_case_solved.tsv \
                                               data/interim/labkey_tiering_highest_tiers.tsv \
											   data/interim/labkey_phenotypes_clean.tsv \
											   data/raw/participant_2023-08-09_12-18-36.tsv \
											   src/data/labkey_merge_clinical_annotations.py
	python3 -m src.data.labkey_merge_clinical_annotations

# Annotate DNMs with constraint and OMIM data
data/interim/dnms_annotated.tsv : data/interim/dnms_38_combined_vep_tidy.tsv \
                                  data/uom_csf/regional_nonsense_constraint.tsv \
								  data/uom_csf/nmd_annotations.tsv \
								  data/uom_csf/gene_ids.tsv \
								  data/uom_csf/genemap2_simple.tsv
	python3 -m src.data.dnms_annotate_constraint

# Merge DNMs with clinical annotations
data/interim/dnms_annotated_clinical.tsv : data/interim/gel_dnm_offspring_clean.tsv \
                                           data/interim/dnms_annotated.tsv \
										   data/interim/labkey_participant_clinical.tsv \
										   src/data/combine_labkey_and_dnm_data_for_clinical_review.py
	python3 -m src.data.combine_labkey_and_dnm_data_for_clinical_review

# Get case solved odds ratios
data/statistics/case_solved_odds_ratios.tsv : src/data/dnms_case_solved_odds.py \
                                              data/interim/dnms_annotated_clinical.tsv
	python3 -m src.data.dnms_case_solved_odds

# Next