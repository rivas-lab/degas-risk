#!/bin/bash
#SBATCH -J DEGWAS
#SBATCH -p owners,mrivas,normal
#SBATCH --mem=51200
#SBATCH --cores=8
#SBATCH -t 2-00:00:00
#SBATCH -o logs/degas_train_gwas_v2_last.%A_%a.out

# ensure usage
if [ $# != 1 ]; then
    if [ -z $SLURM_ARRAY_TASK_ID ]; then
        echo "usage: $0 GBE_ID"; exit 1
    else 
        # phe=$(awk '{for (i=0; i<22; i++){print}}' rerun_train_gwas_20190728.tsv | awk -v nr=$SLURM_ARRAY_TASK_ID '(NR==nr){print $1}')
        # chr=$(expr $(expr $SLURM_ARRAY_TASK_ID % 22) + 1)
	phe=$(awk -v nr=$SLURM_ARRAY_TASK_ID 'NR==nr' ../reference/final_rerun_v2.txt )
    fi
else
    phe=$1
    chr="1-22"
fi

# load an appropriate plink version
if grep -q "CPU_GEN:HSW\|CPU_GEN:BDW\|CPU_GEN:SKX" <(a=$(hostname); sinfo -N -n ${a::-4} --format "%50f"); then
   # AVX2 is suitable for use on this node if CPU is recent enough
   ml load plink2/20190402
else
   ml load plink2/20190402-non-AVX2
fi

# directories and output file prefix, for convenience -- redacted
ukbb_dir=""
proj_dir=""
out_prefix="" 
phe_dirs="" 
phe_file=$(find ${phe_dirs} -name "${phe}.phe")


# loop over variants on one/both arrays
for kind in "one_array" "both_array"; do 
    # which variants to exclude (as oppose to extract)
    if [ $kind == "one_array" ]; then 
        exclude="${ukbb_dir}/sqc/both_array_variants.txt"
        array=""
    else 
        exclude="${ukbb_dir}/sqc/one_array_variants.txt"
        array="Array"
    fi
    # run GWAS
    plink2 --bpfile "${ukbb_dir}/array_combined/pgen/ukb24983_cal_hla_cnv" \
           --chr 1-22 \
           --covar "${ukbb_dir}/sqc/ukb24983_GWAS_covar.phe" \
           --covar-name age sex $array PC1-PC4 \
           --covar-variance-standardize \
           --exclude $exclude \
           --extract "../reference/variant_qc_v2.prune.in" \
           --keep "${proj_dir}/population-split/ukb24983_white_british_train.phe" \
           --remove "${ukbb_dir}/sqc/w24983_20200204.phe" \
           --memory 48000 \
           --pheno ${phe_file} --pheno-quantile-normalize \
           --glm firth-fallback hide-covar omit-ref \
           --threads 8 \
           --out "${out_prefix}.${kind}.${phe}" 
done

# combine summary stats
phe2="PHENO1"
for suffix in "glm.linear" "glm.logistic.hybrid"; do
    if [ -f "${out_prefix}.both_array.${phe}.${phe2}.${suffix}" ]; then
        for kind in "both_array" "one_array"; do
            cat "${out_prefix}.${kind}.${phe}.${phe2}.${suffix}"
        done | sort -k1,1n -k2,2n -u > "${out_prefix}.${phe}.${suffix}"
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

