#!/bin/python
import sys,os
import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD
from scipy.sparse import csr_matrix
from copy import deepcopy

# ensure usage
if len(sys.argv) < 3:
	print("usage: python tsvd.py /path/to/dataset.pkl.gz n_pcs")
	sys.exit(2)
score = True
center = True

# parse
dataset = sys.argv[1]
n = int(sys.argv[2])
dataset_info = os.path.basename(dataset).split('.')[0].split('_')
dataset_name = '_'.join([datum for datum in dataset_info] + [str(n)+'PCs'])

# load, do analysis
data = pd.read_pickle(dataset)

# subset to phenotypes in reference list with at least 2 hits at p < p_star
exclude=pd.read_table("../reference/blacklist.txt", header=None).iloc[:,0].tolist()
phe_ref=pd.read_csv("../reference/final_sumstats_v3.tsv", 
                    header=None, sep=" ").iloc[:,0].tolist()

print(data.shape)
data = data[[i for i in phe_ref if i not in exclude and i in data.columns]]
print(data.shape)
data = data.dropna(axis=0, how='all').dropna(axis=1, how='all')
print(data.shape)

# center data
if center:
    dataset_name = dataset_name.replace('nonCenter','center')
    data = data.subtract(data.mean()).divide(data.std())
# impute zeros -- for some reason df.fillna() leaves some nulls
data = pd.DataFrame(np.nan_to_num(data.values), index=data.index, columns=data.columns)
print(data.shape)

print(np.sum(np.sum(np.asarray_chkfinite(data))))


# do TSVD with these parameters
matt = TruncatedSVD(n_components=n, n_iter=20, random_state=24983)
US = matt.fit_transform(csr_matrix(data.values)) 

# necessary for allele scoring
ukb_root="" # redacted
bim_file=ukb_root+'array_combined/pgen/ukb24983_cal_hla_cnv.pvar'
with open(bim_file, 'r') as f:
    id2alt = {line.split()[2]:line.rstrip().split()[-1] for line in f}

# save the results
np.savez(os.path.join(os.path.dirname(dataset), 'tsvd', dataset_name),
         U = US/matt.singular_values_,
         V = matt.components_.T,
         D = matt.singular_values_,
         variance_explained = matt.explained_variance_,
         variance_explained_ratio = matt.explained_variance_ratio_,
         label_phe_code = np.array(data.columns),
         label_var = np.array(data.index),
         label_var_minor_allele = np.array(data.index.map(id2alt.get))
)

# optional -- do allele scoring
if score:
    proj_dir='' # redacted
    prs_dir=proj_dir+'scorefiles/v2/'
    weights=prs_dir+dataset_name+'.prs.weights.txt'
    pd.DataFrame(np.hstack((np.array(data.index).reshape(-1,1),
                            np.array(data.index.map(id2alt.get)).reshape(-1,1),
                            US))
                ).to_csv(weights, sep='\t', header=False, index=False)
    os.system('ml load plink2')
    os.system(' '.join(['plink2 --bfile ', bim_file[:-5], '--memory 30000', 
                        '--score', weights, 'center', 'cols=nmissallele,denom,scoresums',
                        '--score-col-nums', '3-{}'.format(2+n),
                        '--out', prs_dir+dataset_name]))
