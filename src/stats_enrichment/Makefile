.ONESHELL :
.PHONY : all

SHELL = bash

all : data/statistics/dnms_enrichment_all_genes.tsv \

data/statistics/dnms_enrichment_all_genes.tsv : data/interim/dnms_annotated.tsv \
                                                data/uom_csf/regional_nonsense_constraint.tsv \
                                                data/uom_csf/regional_constraint_stats.tsv \
                                                data/uom_csf/gene_ids.tsv \
                                                data/uom_csf/genemap2_simple.tsv \
												src/stats_enrichment/dnms_enrichment.py
	source activate ukb
	python3 -m src.stats_enrichment.dnms_enrichment
