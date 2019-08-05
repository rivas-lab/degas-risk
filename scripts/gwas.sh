#!/bin/bash
#SBATCH -J DEGWAS
#SBATCH -p owners,mrivas,normal
#SBATCH --mem=25600
#SBATCH --cores=4
#SBATCH -t 1-00:00:00
#SBATCH -o old_logs/degas_train_gwas_chrsplit.%A_%a.out

# ensure usage
if [ $# != 1 ]; then
    if [ -z $SLURM_ARRAY_TASK_ID ]; then
        echo "usage: $0 GBE_ID"; exit 1
    else 
        phe=$(awk '{for (i=0; i<22; i++){print}}' rerun_train_gwas_20190728.tsv | awk -v nr=$SLURM_ARRAY_TASK_ID '(NR==nr){print $1}')
        chr=$(expr $(expr $SLURM_ARRAY_TASK_ID % 22) + 1)
    fi
else
    phe=$1
    chr=22
fi

# load an appropriate plink version
if grep -q "CPU_GEN:HSW\|CPU_GEN:BDW\|CPU_GEN:SKX" <(a=$(hostname); sinfo -N -n ${a::-4} --format "%50f"); then
   # AVX2 is suitable for use on this node if CPU is recent enough
   ml load plink2/20190402
else
   ml load plink2/20190402-non-AVX2
fi

# directories and output file prefix, for convenience
ukbb_dir="/oak/stanford/groups/mrivas/private_data/ukbb/24983"
proj_dir="/oak/stanford/groups/mrivas/projects/degas-risk"
out_prefix="${proj_dir}/summary-stats/train/chrsplit/ukb24983_v2.degas-val.chr${chr}"

# loop over variants on one/both arrays
for kind in "one_array" "both_array"; do 
    # what to do with the variants on one array
    if [ $kind == "one_array" ]; then 
        mode="--extract";
    else 
        mode="--exclude";
    fi
    # run GWAS
    plink2 --bpfile "${ukbb_dir}/array_combined/pgen/ukb24983_cal_hla_cnv" \
           --chr ${chr} \
           --covar "${ukbb_dir}/sqc/ukb24983_GWAS_covar.phe" \
           --covar-name age sex Array PC1-PC4 \
           --covar-variance-standardize \
           $mode "${ukbb_dir}/sqc/one_array_variants.txt" \
           --glm firth-fallback hide-covar omit-ref \
           --keep "${proj_dir}/population-split/ukb24983_white_british_train.phe" \
           --memory 25000 \
           --out "${out_prefix}.${kind}.${phe}" \
           --pheno "${ukbb_dir}/phenotypedata/master_phe/master.phe" \
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

