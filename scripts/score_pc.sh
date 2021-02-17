#!/bin/bash
#SBATCH -p normal,owners,mrivas
#SBATCH --mem=16000
#SBATCH -t 10:00:00
#SBATCH -J scorePCs
#SBATCH -o logs/score_pc.%A_%a.out
set -eo pipefail

# dataset info (path to numpy file, number of components)
npz=$1
npc=500
# note: in this formation, U is variants and V is phenotypes
name=$(basename $npz )

# plink spec
ml load plink
bfile="" # redacted: ukb24983_cal_hla_cnv
outD="" # redacted scorefiles dir
mkdir -p ${outD}/temp_score

# mode 1: score the specified PC
if [ $# -gt 1 ] || [ ! -z "$SLURM_ARRAY_TASK_ID" ] ; then 

if [ $# -gt 1 ]; then
pcx=$2
else
pcx=$SLURM_ARRAY_TASK_ID
fi

pytemp="
import numpy as np;
z = np.load('"$npz"'); 
print('\n'.join(['\t'.join(tup) for tup in zip(*[map(str,l) for l in [z['label_var'],z['label_var_minor_allele'],z['D']["$pcx"-1]*z['U'][:,"$pcx"-1]]])]))"

plink --bfile $bfile \
      --score <( python <( echo $pytemp ) ) sum center double-dosage \
      --out ${outD}/temp_score/${name::-4}_PC${pcx} --memory 16000

# mode 2: score all PCs up to number specified above
else

for pcx in $( seq 1 $npc ); do
continue # only combine
pytemp="
import numpy as np;
z = np.load('"$npz"'); 
print('\n'.join(['\t'.join(tup) for tup in zip(*[map(str,l) for l in [z['label_var'],z['label_var_minor_allele'],z['D']["$pcx"-1]*z['U'][:,"$pcx"-1]]])]))"

plink --bfile $bfile \
      --score <( python <( echo $pytemp ) ) sum center double-dosage \
      --out ${outD}/temp_score/${npz::-4}_PC${pcx} --memory 16000 
done

# combine result files
pytemp="import pandas as pd; print('ok'); 
df = reduce(lambda x,y:pd.merge(x,y,left_index=True,right_index=True), [pd.read_table('"${outD}"/temp_score/"${name::-4}"_PC'+str(i+1)+'.profile', usecols=['IID', 'SCORESUM'], index_col='IID', sep='\s+') for i in range("${npc}")]); 
df.columns = ['SCORE_PC{}'.format(i+1) for i in range("${npc}")]; 
df.to_csv('"${outD}"/"${name::-4}"_full.profile', sep='\t')"

echo $pytemp
# do it then remove the files we made
python <(echo $pytemp)
# rm "temp_score/${npz::-4}_PC*"

fi
