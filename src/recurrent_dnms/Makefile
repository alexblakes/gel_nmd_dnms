.ONESHELL :
.PHONY : all

SHELL = /bin/bash

all : data/final/recurrent_dnms.tsv \
      

data/final/recurrent_dnms.tsv : data/interim/dnms_annotated_clinical.tsv \
                                src/recurrent_dnms/table.py
	python3 -m src.recurrent_dnms.table
