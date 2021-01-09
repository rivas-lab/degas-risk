#!/bin/python3
import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD
from scipy.sparse import csr_matrix
from copy import deepcopy


# shorten paths
projects='/oak/stanford/groups/mrivas/projects/degas-risk/'

# load input genetic data
df=pd.read_pickle('neale_beta_p1e-6.ld-indep.manual_combine.sumstats.pkl.gz')
print(df.shape)

# filter to specified phenotypes
keep=pd.read_table(projects+'neale/phenos_include.tsv', header=None, index_col=0).iloc[:,0]
name=keep.to_dict()
df=df.loc[:,keep.index.astype(str)]

# remove phenotypes with < 2 genetic associations
df=df.loc[:,(~df.isnull()).sum(axis=0) > 1]

# remove variants with no phenotype associations
df=df.loc[(~df.isnull()).sum(axis=1) > 0,:]

# confirm matrix shape, center, fill null values with zeros
print(df.shape)
df=df.subtract(df.mean()).divide(df.std())
df=pd.DataFrame(np.nan_to_num(df.values), index=df.index, columns=df.columns)
print(np.sum(np.sum(np.asarray_chkfinite(df))))


# do tsvd
matt = TruncatedSVD(n_components=500, n_iter=20, random_state=24983)
US = matt.fit_transform(csr_matrix(df.values)) 

# helper to map rsid back to alleles (for prs scoring)
v=pd.read_table('neale_lab_variants_in_ukb.txt', index_col=1, header=None).iloc[:,0].to_dict()

# save the results
np.savez(projects+'datasets/neale/nl-degas_beta_p1e-6_500PCs_20201221.interactive',
         U = US/matt.singular_values_,
         V = matt.components_.T,
         D = matt.singular_values_,
         variance_explained = matt.explained_variance_,
         variance_explained_ratio = matt.explained_variance_ratio_,
         label_phe_code = df.columns.array,
         label_phe_name = np.array([name.get(i) for i in df.columns]),
	 label_var_rsid = df.index.array,
         label_var = np.array([v.get(i) for i in df.index]),
	 label_var_minor_allele = np.array([v.get(i).split(':')[-1] for i in df.index])
)
