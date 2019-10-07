#!/bin/bash
#SBATCH -J dPRS
#SBATCH -p owners,mrivas
#SBATCH -t 5:00:00
#SBATCH --mem=12000
#SBATCH -o logs/phe_test.20191006.%A_%a.out
#SBATCH --array=0-999

# which set of traits do we want?
ix=$(expr ${SLURM_ARRAY_TASK_ID:=1} - 0)

# ok now run them
for phe_id in $(awk -v nr=$ix '(NR > 1 && NR % 1000 == nr){print $1}' ../reference/final_phe_codes.txt); do
    python pipeline.py $phe_id ;
done 
