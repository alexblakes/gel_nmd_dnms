.ONESHELL :
.PHONY : all

SHELL = bash

allele_frequencies = data/interim/allele_frequencies_gnomad_v3.1.1_genomes.vcf.gz \
                     data/interim/allele_frequencies_gel.txt.gz \

all : $(allele_frequencies) \
      data/interim/dnms_38_combined_af.vcf.gz \
	  data/interim/dnms_38_combined_af_vep.vcf.gz \
	  data/interim/dnms_38_combined_af_vep_tidy.tsv \

# Extract allele frequency data from GEL and gnomAD
$(allele_frequencies) &: data/interim/dnms_38_combined.vcf.gz \
                         src/annotate_dnms/allele_frequencies_gnomad.sh \
                         src/annotate_dnms/allele_frequencies_gel.sh
	bsub < src/annotate_dnms/allele_frequencies_gnomad.sh
	bsub < src/annotate_dnms/allele_frequencies_gel.sh
	bwait -w "done(af_gel)"
	bwait -w "done(af_gnomad)"
	# Takes approximately 90 minutes to run

# Annotate DNMs with allele frequency data
data/interim/dnms_38_combined_af.vcf.gz : data/interim/dnms_38_combined.vcf.gz \
                                          $(allele_frequencies) \
                                          src/annotate_dnms/annotate_allele_frequencies.sh
	source activate bio
	bash src/annotate_dnms/annotate_allele_frequencies.sh
	source activate ukb
	# Takes approximately 5 minutes to run

# Annotate DNMs with VEP
data/interim/dnms_38_combined_af_vep.vcf.gz : data/interim/dnms_38_combined_af.vcf.gz \
                                              src/annotate_dnms/vep.sh
	source activate bio
	bash src/annotate_dnms/vep.sh
	source activate ukb
	# Takes about 5 minutes to run

# Filter DNMs by allele frequency. Tidy VEP consequences, cohort and participant ID.
data/interim/dnms_38_combined_af_vep_tidy.tsv : data/interim/dnms_38_combined_af_vep.vcf.gz \
                                                src/annotate_dnms/filter_and_tidy.sh \
												src/annotate_dnms/dnms_tidy_logging.py
	source activate bio
	bash src/annotate_dnms/filter_and_tidy.sh
	source activate ukb
	python3 -m src.annotate_dnms.dnms_tidy_logging