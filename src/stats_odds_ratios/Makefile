.ONESHELL :
.PHONY : all

SHELL = bash

plots = data/plots/case_solved_odds_ratios.svg \
        data/plots/case_solved_odds_ratios.png

all : data/statistics/case_solved_odds_ratios.tsv \
      $(plots) \

# Calculate odds ratio statistics
data/statistics/case_solved_odds_ratios.tsv : data/interim/dnms_annotated_clinical.tsv \
                                              src/stats_odds_ratios/dnms_case_solved_odds.py
	source activate ukb
	python3 -m src.stats_odds_ratios.dnms_case_solved_odds

# Plot odds ratios
$(plots) : data/statistics/case_solved_odds_ratios.tsv \
           src/stats_odds_ratios/plot_odds_ratios.py
	source activate ukb
	python3 -m src.stats_odds_ratios.plot_odds_ratios