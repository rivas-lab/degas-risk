#!/bin/bash

# get all filters
zcat /oak/stanford/groups/mrivas/private_data/ukbb/variant_filtering/variant_filter_table.tsv.gz | awk -F'\t' '($30 != 0){print $5}' > not_all_filters.txt

# MAF > 0.01%, missingness < 5%, all filters for array sites, LD independent
plink2 --pfile /oak/stanford/groups/mrivas/ukbb24983/array_combined/pgen/ukb24983_cal_hla_cnv \
       --maf 1e-4 --geno 0.05 --exclude not_all_filters.txt --indep-pairwise 50 5 0.5 --memory 30000 \
       --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe \
       --out ../reference/variant_qc 
