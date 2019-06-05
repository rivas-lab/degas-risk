#!/bin/python
import sys

# ensure usage
if len(sys.argv) < 3:
	print("usage: python tsvd.py /path/to/dataset.pkl.gz n_pcs")
	sys.exit(2)
cca = True

# parse
import os
dataset = sys.argv[1]
n = int(sys.argv[2])
dataset_info = os.path.basename(dataset).split('.')[0].split('_')
dataset_name = '_'.join([datum for datum in dataset_info] + [str(n)+'PCs'])
dataset_name += '_cca' if cca else ''
print(dataset_name)

# these are useful
import numpy as np
import pandas as pd
from scipy.linalg import inv,sqrtm
from sklearn.decomposition import TruncatedSVD
from scipy.sparse import csr_matrix
phe_corr='/oak/stanford/groups/mrivas/projects/degas-risk/covars/all_white_british_phe_corr.pkl.gz'

# load, do analysis
data = pd.read_pickle(dataset)
print("loaded")
if cca:
    data = data[sorted(data.columns)]
    yty  = pd.read_pickle(phe_corr).sort_index()
    data = data.fillna(value=0).dot(inv(sqrtm(yty + (0.99 * np.identity(yty.shape[0])))))
print("woohoo")
matt = TruncatedSVD(n_components=n, n_iter=20, random_state=24983)
US = matt.fit_transform(csr_matrix(data.fillna(value=0).values))

# save the results
np.savez(os.path.join(os.path.dirname(dataset), dataset_name),
         U = US/matt.singular_values_,
         V = matt.components_.T,
         D = matt.singular_values_,
         variance_explained = matt.explained_variance_,
         variance_explained_ratio = matt.explained_variance_ratio_,
         label_phe_code = np.array(data.columns),
         label_var = np.array(data.index),
)
