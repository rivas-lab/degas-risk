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
cca = True 

# parse
dataset = sys.argv[1]
n = int(sys.argv[2])
dataset_info = os.path.basename(dataset).split('.')[0].split('_')
dataset_name = '_'.join([datum for datum in dataset_info] + [str(n)+'PCs'])
dataset_name += '_cca' if cca else ''

# these are useful
phe_corr='/oak/stanford/groups/mrivas/projects/degas-risk/covars/all_white_british_phe_corr.pkl.gz'
bim_file='/oak/stanford/groups/mrivas/ukbb24983/array_combined/pgen/ukb24983_cal_hla_cnv.pvar'

# load, do analysis
data = pd.read_pickle(dataset).fillna(value=0)
# process input -- this is a helper for CCA
def mat_sqrt_inv(x):
    # x = a.diag(d).aT
    d,a = np.linalg.eig(x)
    d = np.diag(map(lambda i:i**(-1/2) if i > 0 else 0, d))
    return a.dot(d).dot(a)
if cca:
    phes = deepcopy(sorted(data.columns))
    data = data[phes]
    yty  = pd.read_pickle(phe_corr).sort_index().fillna(value=0)
    yty  = yty[phes]
    data = data.dot(mat_sqrt_inv(yty + (0.99 * np.eye(yty.shape[0]))))
    data.columns = phes

# parameters
matt = TruncatedSVD(n_components=n, n_iter=20, random_state=24983)
# CCA isn't sparse because of the matrix multiplication above
US = matt.fit_transform(data.values if cca else csr_matrix(data.values)) 

# necessary for allele scoring
with open(bim_file, 'r') as f:
    id2alt = {line.split()[2]:line.rstrip().split()[-1] for line in f}

# save the results
np.savez(os.path.join(os.path.dirname(dataset), dataset_name),
         U = US/matt.singular_values_,
         V = matt.components_.T,
         D = matt.singular_values_,
         variance_explained = matt.explained_variance_,
         variance_explained_ratio = matt.explained_variance_ratio_,
         label_phe_code = np.array(data.columns),
         label_var = np.array(data.index),
         label_var_minor_allele = np.array(data.index.map(id2alt.get))
)
