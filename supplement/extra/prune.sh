#!/bin/bash
#SBATCH -J ld
#SBATCH --mem=16000
#SBATCH -t 6:00:00
#SBATCH -p mrivas,normal,owners
#SBATCH -o logs/batch.neale.ld.chr%a.%A.out
#SBATCH --array=1-22

# software
ml load biology; ml load plink2
which plink2

# shortnames
ukb="/oak/stanford/groups/mrivas/ukbb24983"
repo="/oak/stanford/groups/mrivas/users/magu/repos/rivas-lab/degas-risk"


# check if we're in a batch (array) job
if -z $SLURM_ARRAY_TASK_ID; then
	# if so, loop over all chromosomes
	c1=1
	c2=22
else
	# otherwise just do the one specified by the array index
	c1=$SLURM_ARRAY_TASK_ID
	c2=$SLURM_ARRAY_TASK_ID
fi


# loop 
for c in $( seq $c1 $c2 ); do 
	# ld prune
	plink2 --pgen ${ukb}/imp/pgen/ukb24983_imp_chr${c}_v3.pgen \
	       --psam ${ukb}/imp/pgen/ukb24983_imp_chr${c}_v3.psam \
	       --pvar ${ukb}/imp/pgen/ukb24983_imp_chr${c}_v3.pvar.zst \
	       --keep ${ukb}/sqc/population_stratification/ukb24983_white_british.phe \
	       --extract <( cut -f1 ${repo}/scripts/neale_lab_variants_in_ukb.txt ) \
	       --indep-pairwise 50 5 0.5 --memory 16000 --out neale_ld_indep.chr${c} 
done 
