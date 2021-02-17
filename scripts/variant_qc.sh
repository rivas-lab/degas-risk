#!/bin/bash

# get "all" filters (PTV paper)
ukbb_dir=""
ukb2_dir=""
zcat ${ukbb_dir}/variant_filtering/variant_filter_table.tsv.gz | awk -F'\t' '($30 != 0){print $5}' > not_all_filters.txt


# 1a. array-specific LD pruning and missingness (Axiom)
plink2 --pfile ${ukb2_dir}/cal/pgen/ukb24983_cal_cALL_v2 \
       --keep ${ukb2_dir}/sqc/population_stratification/ukb24983_white_british.phe \
       --keep-fam ${ukb2_dir}/sqc/axiom_individuals.txt \
       --exclude not_all_filters.txt \
       --extract ${ukb2_dir}/sqc/axiom_specific_variants.txt \
       --maf 1e-4 --geno 0.05 \
       --indep-pairwise 50 5 0.5 \
       --out ../reference/temp_axiom
# keep-fam is equivalent to keep, and is used because there is no --keep-intersect

# 1b. array-specific LD pruning and missingness (BiLEVE)
plink2 --pfile ${ukb2_dir}/cal/pgen/ukb24983_cal_cALL_v2 \
       --keep ${ukb2_dir}/sqc/population_stratification/ukb24983_white_british.phe \
       --keep-fam ${ukb2_dir}/sqc/bileve_individuals.txt \
       --exclude not_all_filters.txt \
       --extract${ukb2_dir}/sqc/bileve_specific_variants.txt \
       --maf 1e-4 --geno 0.05 \
       --indep-pairwise 50 5 0.5 \
       --out ../reference/temp_bileve

# 1c. array-specific LD pruning and missingness (CNV)
plink2 --pfile ${ukb2_dir}/cnv/pgen/cnv \
       --keep ${ukb2_dir}/population_stratification/ukb24983_white_british.phe \
       --maf 1e-4 --geno 0.05 \
       --indep-pairwise 50 5 0.5 \
       --out ../reference/temp_cnv

# 1d. array-specific LD pruning and missingness (HLA)
plink2 --pfile ${ukb2_dir}/hla/pgen/ukb_hla_v3 \
       --keep ${ukb2_dir}/sqc/population_stratification/ukb24983_white_british.phe \
       --maf 1e-4 --geno 0.05 \
       --indep-pairwise 50 5 0.5 \
       --out ../reference/temp_hla

# combine 
cat ../reference/temp_*.prune.in > ../reference/array_specific_filters.txt



# 2. missingness and AF for variants on both arrays
plink2 --pfile ${ukb2_dir}/cal/pgen/ukb24983_cal_cALL_v2 \
       --keep ${ukb2_dir}/sqc/population_stratification/ukb24983_white_british.phe \
       --exclude not_all_filters.txt \
       --extract ${ukb2_dir}/sqc/both_array_variants.txt \
       --maf 1e-4 --geno 0.05 --write-snplist \
       --out ../reference/shared_filter



# 3. Combine array/cnv/hla-specific filters, and missingness from shared variants
#     then add a more aggressitve MAF cutoff and re-do LD pruning 
plink2 --pfile ${ukb2_dir}/array_combined/pgen/ukb24983_cal_hla_cnv \
       --keep ${ukb2_dir}/sqc/population_stratification/ukb24983_white_british.phe \
       --exclude not_all_filters.txt \
       --extract ../reference/array_specific_filters.txt ../reference/shared_filter.snplist \
       --indep-pairwise 50 5 0.5 \
       --maf 1e-4 \
       --memory 30000 \
       --out ../reference/variant_qc_v2

# clean up intermediates
rm ../reference/temp_*

