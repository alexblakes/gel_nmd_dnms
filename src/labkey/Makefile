.ONESHELL :
.PHONY : all

SHELL = bash
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

exit_questionnaire = data/interim/labkey_exit_questionnaires_case_solved.tsv \
                     data/interim/labkey_exit_questionnaire_acmg.tsv \

tiering = data/interim/labkey_tiering_clean.tsv \
          data/interim/labkey_tiering_highest_tiers.tsv \

all : data/interim/gel_dnm_offspring_clean.tsv \
      $(exit_questionnaire) \
      data/interim/labkey_phenotypes_clean.tsv \
      $(tiering) \
      data/interim/labkey_participant_clinical.tsv \

# Clean GEL DNM cohort
data/interim/gel_dnm_offspring_clean.tsv : src/labkey/gel_dnms_offspring.py
	python3 -m src.labkey.gel_dnms_offspring

# Tidy exit questionnaire data
$(exit_questionnaire) &: data/raw/gmc_exit_questionnaire_2023-08-01_09-42-17.tsv \
                         src/labkey/labkey_tidy_exit_questionnaires.py
	python3 -m src.labkey.labkey_tidy_exit_questionnaires

# Tidy labkey phenotype data
data/interim/labkey_phenotypes_clean.tsv : src/labkey/labkey_tidy_phenotypes.py \
                                           data/raw/rare_diseases_participant_phen_2023-08-01_09-40-58.tsv
	python3 -m src.labkey.labkey_tidy_phenotypes

# Tidy tiering data
$(tiering) &: src/labkey/labkey_tidy_tiering.py \
             data/raw/tiering_data_2023-08-01_09-44-03.tsv
	python3 -m src.labkey.labkey_tidy_tiering

# Merge labkey annotations
data/interim/labkey_participant_clinical.tsv : data/interim/labkey_exit_questionnaires_case_solved.tsv \
                                               data/interim/labkey_tiering_highest_tiers.tsv \
						     data/interim/labkey_phenotypes_clean.tsv \
						     data/raw/participant_2023-08-09_12-18-36.tsv \
						     src/labkey/labkey_merge_clinical_annotations.py
	python3 -m src.labkey.labkey_merge_clinical_annotations
