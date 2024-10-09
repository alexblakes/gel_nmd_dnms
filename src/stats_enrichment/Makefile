.ONESHELL :
.PHONY : all

SHELL = bash

enrichment_stats = data/interim/dnms_enrichment_counts_all_genes.tsv \
                   data/interim/dnms_enrichment_counts_ad.tsv \
                   data/interim/dnms_enrichment_counts_ar.tsv \
                   data/interim/dnms_enrichment_counts_non_morbid.tsv

enrichment_stats_tidy = data/statistics/dnms_enrichment_all_genes_tidy.tsv \
                        data/statistics/dnms_enrichment_ad_tidy.tsv \
                        data/statistics/dnms_enrichment_ar_tidy.tsv \
                        data/statistics/dnms_enrichment_non_morbid_tidy.tsv \

plots = data/plots/dnm_enrichment.png \
        data/plots/dnm_enrichment.svg \

all : data/interim/dnms_enrichment_counts.tsv \
      $(enrichment_stats) \
      $(enrichment_stats_tidy) \
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

# Tidy enrichment stats for plotting
$(enrichment_stats_tidy) &: $(enrichment_stats) \
                            src/stats_enrichment/statistics_for_plot.py
	source activate ukb
	python3 -m src.stats_enrichment.statistics_for_plot

# Plot DNM enrichment statistics
$(plots) &: $(enrichment_stats_tidy) \
            src/stats_enrichment/plot.py
	source activate ukb
	python3 -m src.stats_enrichment.plot