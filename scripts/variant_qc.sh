#!/bin/bash

# get "all" filters (PTV paper)
zcat /oak/stanford/groups/mrivas/private_data/ukbb/variant_filtering/variant_filter_table.tsv.gz | awk -F'\t' '($30 != 0){print $5}' > not_all_filters.txt


# 1a. array-specific LD pruning and missingness (Axiom)
plink2 --pfile /oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_cALL_v2 \
       --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe \
       --keep-fam /oak/stanford/groups/mrivas/ukbb24983/sqc/axiom_individuals.txt \
       --exclude not_all_filters.txt \
       --extract /oak/stanford/groups/mrivas/ukbb24983/sqc/axiom_specific_variants.txt \
       --maf 1e-4 --geno 0.05 \
       --indep-pairwise 50 5 0.5 \
       --out ../reference/temp_axiom
# keep-fam is equivalent to keep, and is used because there is no --keep-intersect

# 1b. array-specific LD pruning and missingness (BiLEVE)
plink2 --pfile /oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_cALL_v2 \
       --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe \
       --keep-fam /oak/stanford/groups/mrivas/ukbb24983/sqc/bileve_individuals.txt \
       --exclude not_all_filters.txt \
       --extract /oak/stanford/groups/mrivas/ukbb24983/sqc/bileve_specific_variants.txt \
       --maf 1e-4 --geno 0.05 \
       --indep-pairwise 50 5 0.5 \
       --out ../reference/temp_bileve

# 1c. array-specific LD pruning and missingness (CNV)
plink2 --pfile /oak/stanford/groups/mrivas/ukbb24983/cnv/pgen/cnv \
       --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe \
       --maf 1e-4 --geno 0.05 \
       --indep-pairwise 50 5 0.5 \
       --out ../reference/temp_cnv

# 1d. array-specific LD pruning and missingness (HLA)
plink2 --pfile /oak/stanford/groups/mrivas/ukbb24983/hla/pgen/ukb_hla_v3 \
       --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe \
       --maf 1e-4 --geno 0.05 \
       --indep-pairwise 50 5 0.5 \
       --out ../reference/temp_hla

# combine 
cat ../reference/temp_*.prune.in > ../reference/array_specific_filters.txt



# 2. missingness and AF for variants on both arrays
plink2 --pfile /oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_cALL_v2 \
       --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe \
       --exclude not_all_filters.txt \
       --extract /oak/stanford/groups/mrivas/ukbb24983/sqc/both_array_variants.txt \
       --maf 1e-4 --geno 0.05 --write-snplist \
       --out ../reference/shared_filter



# 3. Combine array/cnv/hla-specific filters, and missingness from shared variants
#     then add a more aggressitve MAF cutoff and re-do LD pruning 
plink2 --pfile /oak/stanford/groups/mrivas/ukbb24983/array_combined/pgen/ukb24983_cal_hla_cnv \
       --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe \
       --exclude not_all_filters.txt \
       --extract ../reference/array_specific_filters.txt ../reference/shared_filter.snplist \
       --indep-pairwise 50 5 0.5 \
       --maf 1e-4 \
       --memory 30000 \
       --out ../reference/variant_qc_v2

# clean up intermediates
rm ../reference/temp_*

