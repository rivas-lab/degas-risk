#!/bin/bash
#SBATCH -J dPRS
#SBATCH -p owners,mrivas,jpriest
#SBATCH -t 2:00:00
#SBATCH --mem=12000
#SBATCH -o logs/phe_test.20191011.%A_%a.out
#SBATCH --array=252,265,318,331,337,398,420,425,807,808,976,979 #,1825,1868
# #SBATCH --array=825,868

# which set of traits do we want?
ix=$(expr ${SLURM_ARRAY_TASK_ID:=1} + 0)

# ok now run them
for phe_id in $(awk -v nr=$ix '(NR == nr){print $1}' ../reference/final_phe_codes.txt); do
    echo $phe_id
    python pipeline.py $phe_id ;
done 
