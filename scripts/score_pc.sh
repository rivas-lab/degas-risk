#!/usr/bin/bash
#SBATCH -p normal,owners,mrivas
#SBATCH --mem=16000
#SBATCH -t 2:00:00
#SBATCH -J scorePCs
#SBATCH -o logs/score_pc.%A_%a.out
set -euo pipefail

# dataset info (path to numpy file, number of components)
npz=$1
npc=500
# note: in this formation, U is variants and V is phenotypes

# plink spec
ml load plink
bgen=$OAK/ukbb24983/cal/pgen/ukb24983_cal_cALL_v2
outD=$OAK/projects/degas-risk/scorefiles/
mkdir -p ${outD}/temp_score

# mode 1: score the specified PC
if [ $# -gt 1 ] || [ ! -z $SLURM_ARRAY_TASK_ID ] ; then 

if [ $# -gt 1 ]; then
pcx=$2
else
pcx=$SLURM_ARRAY_TASK_ID
fi

pytemp="
import numpy as np;
z = np.load('"$npz"'); 
print('\n'.join(['\t'.join(tup) for tup in zip(*[map(str,l) for l in [z['label_var'],z['label_var_minor_allele'],z['D']["$pcx"-1]*z['U'][:,"$pcx"-1]]])]))"

name=$(basename $npz )
plink --bfile "$OAK/private_data/ukbb/24983/cal/pgen/ukb24983_cal_cALL_v2" \
      --score <( python <( echo $pytemp ) ) sum center double-dosage \
      --out ${outD}/temp_score/${name::-4}_PC${pcx}

# mode 2: score all PCs up to number specified above
else

for pcx in $( seq 1 $npc ); do
pytemp="
import numpy as np;
z = np.load('"$npz"'); 
print('\n'.join(['\t'.join(tup) for tup in zip(*[map(str,l) for l in [z['label_var'],z['label_var_minor_allele'],z['D']["$pcx"-1]*z['U'][:,"$pcx"-1]]])]))"

name=$(basename $npz )
plink --bfile "$OAK/private_data/ukbb/24983/cal/pgen/ukb24983_cal_cALL_v2" \
      --score <( python <( echo $pytemp ) ) sum center double-dosage \
      --out ${outD}/temp_score/${npz::-4}_PC${pcx} 
done

# combine result files
pytemp="import pandas as pd; 
df = reduce(lambda x,y:pd.merge(x,y,left_index=True,right_index=True), [pd.read_table('temp_score/${npz::-4}_PC'+str(i+1)+'.profile', usecols=['IID', 'SCORESUM'], index_col='IID', sep='\s+') for i in range("${npc}")]); 
df.columns = ['SCORE_PC{}'.format(i+1) for i in range("${npc}")]; 
df.to_csv('"${npz::-4}"_full.profile', sep='\t')"

# do it then remove the files we made
python <(echo $pytemp)
rm "temp_score/${npz::-4}_PC*"

fi
