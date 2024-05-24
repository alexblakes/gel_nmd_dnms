.ONESHELL :
.PHONY : all

SHELL = /bin/bash

all : data/final/recurrent_dnms.tsv \
      data/plots/recurrent_dnms.png \

# Create a table of recurrent DNMs in constrained regions
data/final/recurrent_dnms.tsv : data/interim/dnms_annotated_clinical.tsv \
                                src/stats_recurrent_dnms/table.py
	source activate ukb
      python3 -m src.stats_recurrent_dnms.table

# Plot recurrent DNMs in constrained regions
data/plots/recurrent_dnms.png : data/interim/dnms_annotated_clinical.tsv \
                                src/stats_recurrent_dnms/plot.py
	source activate ukb
      python3 -m src.stats_recurrent_dnms.plot
