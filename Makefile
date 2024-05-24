SHELL = bash

.ONESHELL :

.PHONY : all recurrent_dnms

dnms :
	make -f src/dnms/Makefile all

annotate_dnms :
	make -f src/annotate_dnms/Makefile all

labkey :
	make -f src/labkey/Makefile all

merge_annotations :
	make -f src/merge_annotations/Makefile all

stats_recurrent_dnms :
	make -f src/stats_recurrent_dnms/Makefile all

stats_enrichment :
	make -f src/stats_enrichment/Makefile all

all : dnms \
      annotate_dnms \
	  labkey \
	  merge_annotations \
	  stats_recurrent_dnms \
	  stats_enrichment \
	  data/statistics/case_solved_odds_ratios.tsv \

# Get case solved odds ratios
data/statistics/case_solved_odds_ratios.tsv : src/data/dnms_case_solved_odds.py \
                                              data/interim/dnms_annotated_clinical.tsv
	python3 -m src.data.dnms_case_solved_odds

# Next