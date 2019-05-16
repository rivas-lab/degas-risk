#!/bin/bash
#SBATCH -J DEGWAS
#SBATCH -p normal,owners,mrivas,jpriest
#SBATCH --mem=26000
#SBATCH --cores=4
#SBATCH -t 1-00:00:00
#SBATCH -o logs/degas_train_gwas.%A_%a.out

# ensure usage
if [ $# != 1 ]; then
    if [ -z $SLURM_ARRAY_TASK_ID ]; then
        echo "usage: $0 GBE_ID"; exit 1
    else 
        phe=$(awk -v nr=$SLURM_ARRAY_TASK_ID '(NR==nr+1000){print $1}' ../reference/phenotypes.tsv )
    fi
else
    phe=$1
fi

# loop over variants on one/both arrays
out_prefix="/oak/stanford/groups/mrivas/projects/degas-risk/summary-stats/train/ukb24983_v2.degas-val"
for kind in "one_array" "both_array"; do 
    # run GWAS
    plink2 --bpfile /oak/stanford/groups/mrivas/private_data/ukbb/24983/cal/pgen/ukb24983_cal_cALL_v2 \
           --chr 1-22 \
           --covar /oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_GWAS_covar.phe \
           --covar-name age sex Array PC1-PC4 \
           --extract /oak/stanford/groups/mrivas/ukbb24983/sqc/${kind}_variants.txt \
           --glm firth-fallback hide-covar omit-ref \
           --keep /oak/stanford/groups/mrivas/projects/degas-risk/population-split/ukb24983_white_british_train.phe \
           --memory 25600 \
           --out "${out_prefix}.${kind}.${phe}" \
           --pheno /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.phe \
           --pheno-name $phe \
           --pheno-quantile-normalize \
           --threads 4
done

# combine summary stats
for suffix in "${phe}.glm.linear" "${phe}.glm.logistic.hybrid"; do
    if [ -f "${out_prefix}.both_array.${phe}.${suffix}" ]; then
        for kind in "both_array" "one_array"; do
            cat "${out_prefix}.${kind}.${phe}.${suffix}"
        done | sort -k1,1n -k2,2n -u > "${out_prefix}.${suffix}"
    fi
done

# combine logs
for kind in "both_array" "one_array"; do
    cat "${out_prefix}.${kind}.${phe}.log"
done > "${out_prefix}.${phe}.log"

# remove intermediates
for kind in "both_array" "one_array"; do
    rm ${out_prefix}.${kind}.${phe}.*
done

