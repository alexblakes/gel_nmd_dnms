.ONESHELL :
.PHONY : all

SHELL = bash

all : data/statistics/case_solved_odds_ratios.tsv \

data/statistics/case_solved_odds_ratios.tsv : data/interim/dnms_annotated_clinical.tsv \
                                              src/stats_odds_ratios/dnms_case_solved_odds.py
	source activate ukb
	python3 -m src.stats_odds_ratios.dnms_case_solved_odds