.ONESHELL :
.PHONY : all

SHELL = bash

enrichment_stats = data/interim/dnms_enrichment_counts_all_genes.tsv \
                   data/interim/dnms_enrichment_counts_ad.tsv \
                   data/interim/dnms_enrichment_counts_ar.tsv \
                   data/interim/dnms_enrichment_counts_non_morbid.tsv

plots = data/plots/dnm_enrichment.png \
        data/plots/dnm_enrichment.svg \

all : data/interim/dnms_enrichment_counts.tsv \
      $(enrichment_stats) \
      $(plots) \

# Annotate DNMs with OMIM and constraint data
data/interim/dnms_enrichment_counts.tsv : data/interim/dnms_annotated.tsv \
                                          data/uom_csf/regional_nonsense_constraint.tsv \
                                          data/uom_csf/regional_constraint_stats.tsv \
                                          data/uom_csf/gene_ids.tsv \
                                          data/uom_csf/genemap2_simple.tsv \
                                          src/stats_enrichment/merge_annotations.py
	source activate ukb
	python3 -m src.stats_enrichment.merge_annotations


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