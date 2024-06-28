.ONESHELL :
.PHONY : all

SHELL = bash

gel_vcfs = data/interim/dnms_37_gel.vcf.gz \
           data/interim/dnms_38_gel_raw.vcf.gz \

liftover = data/interim/dnms_38_kaplanis.vcf.gz \
           data/interim/dnms_38_gel_lifted.vcf.gz \

all : data/interim/kaplanis_dnms_filled.tsv \
      data/interim/dnms_37_kaplanis.vcf.gz \
      $(gel_vcfs) \
      $(liftover) \
      data/interim/dnms_38_combined.vcf.gz \

data/interim/kaplanis_dnms_filled.tsv : data/raw/kaplanis_dnms.tsv \
                                        src/dnms/tidy_kaplanis_dnms.py
	source activate dnms
	python3 -m src.dnms.tidy_kaplanis_dnms

data/interim/dnms_37_kaplanis.vcf.gz : data/interim/kaplanis_dnms_filled.tsv \
                                    src/dnms/dnms_to_vcf_ddd.sh
	# Converting Kaplanis DNMs to VCF
	source activate bio
	bash src/dnms/dnms_to_vcf_ddd.sh
    
$(gel_vcfs) &: data/raw/denovo_flagged_variants_2023-05-16_09-51-50.tsv \
               src/dnms/dnms_to_vcf_gel.sh
	source activate bio
	bash src/dnms/dnms_to_vcf_gel.sh

$(liftover) &: data/interim/dnms_37_kaplanis.vcf.gz \
               data/interim/dnms_37_gel.vcf.gz \
               src/dnms/picard_liftover.sh
	source activate bio
	bash src/dnms/picard_liftover.sh data/interim/dnms_37_kaplanis.vcf.gz data/interim/dnms_38_kaplanis.vcf.gz
	bash src/dnms/picard_liftover.sh data/interim/dnms_37_gel.vcf.gz data/interim/dnms_38_gel_lifted.vcf.gz

data/interim/dnms_38_combined.vcf.gz : $(gel_vcfs) \
                                       $(liftover) \
                                       src/dnms/combine_dnms.sh
	source activate bio
	bash src/dnms/combine_dnms.sh
