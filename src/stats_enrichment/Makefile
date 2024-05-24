.ONESHELL :
.PHONY : all

SHELL = bash

enrichment_stats = data/statistics/dnms_enrichment_all_genes.tsv \
                   data/statistics/dnms_enrichment_morbid_genes.tsv \
                   data/statistics/dnms_enrichment_non_morbid_genes.tsv \

plots = data/plots/dnm_enrichment_constrained_regions.png \
        data/plots/dnm_enrichment_constrained_regions.svg \

all : $(enrichment_stats) \
      $(plots) \

# Get DNM enrichment statistics
$(enrichment_stats) &: data/interim/dnms_annotated.tsv \
                       data/uom_csf/regional_nonsense_constraint.tsv \
                       data/uom_csf/regional_constraint_stats.tsv \
                       data/uom_csf/gene_ids.tsv \
                       data/uom_csf/genemap2_simple.tsv \
                       src/stats_enrichment/dnms_enrichment.py
	source activate ukb
	python3 -m src.stats_enrichment.dnms_enrichment

# Plot DNM enrichment statistics
$(plots) &: $(enrichment_stats) \
            src/stats_enrichment/plot_enrichment.py
	source activate ukb
	python3 -m src.stats_enrichment.plot_enrichment