SHELL = bash

.ONESHELL :

.PHONY : all \
         count_grch37_trios \
		 count_dnms \

all : count_grch37_trios \
      count_dnms \

count_grch37_trios : 
	source activate ukb
	python3 -m src.statistics_for_paper.count_gel_trios_in_grch37

count_dnms :
	bash src/statistics_for_paper/count_dnms.sh
