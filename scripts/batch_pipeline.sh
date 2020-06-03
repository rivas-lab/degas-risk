#!/bin/bash
#SBATCH -J dPRS
#SBATCH -p owners,mrivas
#SBATCH -t 2:00:00
#SBATCH --mem=16000
#SBATCH -o logs/phe_pipe.20200506.%A_%a.out
#SBATCH --array=1-667

# which set of traits do we want?
ix=$(expr ${SLURM_ARRAY_TASK_ID:=1} + 0)

# ok now run them
for phe_id in $(awk -v nr=$ix '(NR % 1000 == nr){print $1}' ../reference/final_phe_codes_v2.txt); do
    echo $phe_id
    python pipeline.py $phe_id ;
done 
