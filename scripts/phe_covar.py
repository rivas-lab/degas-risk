#!/bin/python
import numpy as np
import pandas as pd

phe_in='/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.phe'
phenos=pd.read_table('../reference/phenotypes.tsv').iloc[:,0].tolist()

phe = pd.read_table(phe_in, usecols=['IID'] + phenos, index_col='IID', na_values=-9)

# subset to individuals of interest
pop='/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe'
pop_name='all_white_british'
with open(pop, 'r') as f:
    inds = list(filter(lambda x: x > 0, [int(line.split()[0]) for line in f]))
phe = phe.loc[inds,:]

# re-binarize phenotypes
for pheno in phe.columns:
    if len(phe[pheno].value_counts()) == 2:
        phe[pheno] -= 1

# compute correlation
corr = phe.corr(method='pearson')

# write to file
out_prefix='/oak/stanford/groups/mrivas/projects/degas-risk/covars/' + pop_name + '_phe_corr'
corr.to_pickle(out_prefix + '.pkl.gz')
corr.to_csv(out_prefix + '.txt', sep='\t') 
