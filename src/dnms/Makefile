.ONESHELL :
.PHONY : all

SHELL = bash
CONDA_ACTIVATE := source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

gel_vcfs = data/interim/dnms_37_gel.vcf \
           data/interim/dnms_38_gel_raw.vcf \
           data/interim/dnms_38_gel_raw.vcf.gz \
           data/interim/dnms_38_gel_raw.vcf.gz.tbi \

liftover = data/interim/dnms_38_kaplanis.vcf \
           data/interim/dnms_38_kaplanis.vcf.gz \
           data/interim/dnms_38_kaplanis.vcf.gz.tbi \
           data/interim/dnms_38_gel_lifted.vcf \
           data/interim/dnms_38_gel_lifted.vcf.gz \
           data/interim/dnms_38_gel_lifted.vcf.gz.tbi \

gel_combined = data/interim/dnms_38_gel_combined.vcf \
               data/interim/dnms_38_gel_combined.vcf.gz \
               data/interim/dnms_38_gel_combined.vcf.gz.tbi \

all : data/interim/kaplanis_dnms_filled.tsv \
      data/interim/dnms_37_kaplanis.vcf \
      $(gel_vcfs) \
      $(liftover) \
      $(gel_combined) \


data/interim/kaplanis_dnms_filled.tsv : data/raw/kaplanis_dnms.tsv \
                                        src/dnms/tidy_kaplanis_dnms.py
	source activate dnms
	python3 -m src.dnms.tidy_kaplanis_dnms
	source activate ukb

data/interim/dnms_37_kaplanis.vcf : data/interim/kaplanis_dnms_filled.tsv \
                                    src/dnms/dnms_to_vcf_ddd.sh
	# Converting Kaplanis DNMs to VCF
	$(CONDA_ACTIVATE) bio
	bash src/dnms/dnms_to_vcf_ddd.sh data/interim/kaplanis_dnms_filled.tsv data/interim/dnms_37_kaplanis.vcf
	$(CONDA_ACTIVATE) dnms
    
$(gel_vcfs) &: data/raw/denovo_flagged_variants_2023-05-16_09-51-50.tsv \
              src/dnms/dnms_to_vcf_gel.sh
	$(CONDA_ACTIVATE) bio
	bash src/dnms/dnms_to_vcf_gel.sh
	$(CONDA_ACTIVATE) dnms

$(liftover) &: data/interim/dnms_37_kaplanis.vcf \
              data/interim/dnms_37_gel.vcf \
              src/dnms/picard_liftover.sh \
              src/dnms/picard_liftover_submit.sh
	source activate bio
	bash src/dnms/picard_liftover.sh data/interim/dnms_37_kaplanis.vcf data/interim/dnms_38_kaplanis.vcf
	bash src/dnms/picard_liftover.sh data/interim/dnms_37_gel.vcf data/interim/dnms_38_gel_lifted.vcf
	source activate dnms
    
$(gel_combined) &: $(gel_vcfs) \
                  $(liftover) \
                  src/dnms/dnms_combine_gel.sh
	source activate bio
	bash src/dnms/dnms_combine_gel.sh
	source activate dnms