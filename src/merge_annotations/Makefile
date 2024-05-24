.ONESHELL :
.PHONY : all

SHELL = bash

all : data/interim/dnms_annotated.tsv \
      data/interim/dnms_annotated_clinical.tsv \

# Annotate DNMs with constraint and OMIM data
data/interim/dnms_annotated.tsv : data/interim/dnms_38_combined_vep_tidy.tsv \
                                  data/uom_csf/regional_nonsense_constraint.tsv \
								  data/uom_csf/nmd_annotations.tsv \
								  data/uom_csf/gene_ids.tsv \
								  data/uom_csf/genemap2_simple.tsv \
								  src/merge_annotations/dnms_annotate_constraint.py
	source activate ukb
	python3 -m src.merge_annotations.dnms_annotate_constraint
	# Takes about 5 minutes to run

# Merge DNMs with clinical annotations
data/interim/dnms_annotated_clinical.tsv : data/interim/gel_dnm_offspring_clean.tsv \
                                           data/interim/dnms_annotated.tsv \
										   data/interim/labkey_participant_clinical.tsv \
										   src/merge_annotations/combine_labkey_and_dnm_data_for_clinical_review.py
	source activate ukb
	python3 -m src.merge_annotations.combine_labkey_and_dnm_data_for_clinical_review
