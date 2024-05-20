.ONESHELL :
.PHONY : all

SHELL = bash
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

exit_questionnaire = data/interim/labkey_exit_questionnaires_case_solved.tsv \
                     data/interim/labkey_exit_questionnaire_acmg.tsv \

all : $(exit_questionnaire) \
       \

$(exit_questionnaire) : data/raw/gmc_exit_questionnaire_2023-08-01_09-42-17.tsv \
                        src/labkey/labkey_tidy_exit_questionnaires.py
	python3 -m src.labkey.labkey_tidy_exit_questionnaires

