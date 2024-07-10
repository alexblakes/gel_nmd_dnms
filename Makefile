SHELL = bash

.ONESHELL :

.PHONY : dnms annotate_dnms labkey merge_annotations stats_recurrent_dnms stats_enrichment stats_odds_ratios all

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

stats_odds_ratios : 
	make -f src/stats_odds_ratios/Makefile all

all : dnms \
      annotate_dnms \
	  labkey \
	  merge_annotations \
	  stats_recurrent_dnms \
	  stats_enrichment \
	  stats_odds_ratios \
