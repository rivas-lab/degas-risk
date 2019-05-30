#!/bin/python
import glob
import numpy as np
import pandas as pd

# pvalue cutoff for var::phe association
p = 0.01
# z-statistic instead of beta?
z = True
# z-score values by phenotype?
c = False

# find files to import
with open('../reference/summary_stats_gbe.tsv', 'r') as f:
	file_ref = [(line.split()[0],line.rstrip().split()[1:]) for line in f]

# only keep variants not in LD, MAF > 0.01%, and QC'd
with open('../reference/variant_qc.prune.in', 'r') as f:
	var_set = set([line.rstrip() for line in f])

# get alt alleles -- need to standardize summary stats ALT allele
var_alt = {var:None for var in var_set}
with open('/oak/stanford/groups/mrivas/ukbb24983/array_combined/pgen/ukb24983_cal_hla_cnv.pvar', 'r') as f:
	for line in f:
		chrom,pos,var,ref,alt = line.rstrip().split()
		var_alt[var] = alt	

# helper function to load files
def get_summary_stats(files, p=0.01, z=True):
	# load data, perform p-value and LD pruning
	binary,qt = 'logistic' in files[0], 'linear' in files[0]
	maxerr = 0.8 if qt else 0.2 if binary else 0 # error catch
	x = pd.concat([pd.read_table(f, index_col=2) for f in files]).dropna()
	x = x[(x['P'] < p) & (x['TEST'] == "ADD") & (x.index.isin(var_set)) & (x['SE'] < maxerr)]
	x['SIGN'] = 2*x['ALT'].eq(x['A1']) - 1 # 1 if ALT is A1 else -1
	# take LOR if binary, manage z-scoring
	if qt:
		return x['T_STAT' if z else 'BETA'].divide(x['SIGN'])
	elif binary:
		return x['Z_STAT' if z else 'OR'].apply(np.log if z else lambda x:x).divide(x['SIGN'])
	else:
		exit(1)

# load data one phenotype at a time
data = pd.SparseDataFrame(index=var_set)
for phe_id,files in file_ref:
	data[phe_id] = get_summary_stats(files, p=p, z=z)
# center betas/z-stats by rows (phenotype-level normalization)
if c:
	data = data.subtract(data.mean()).divide(data.std())

# name and save
dataset_name = '_'.join(('all',
                         'z' if z else 'beta',
                         'center' if c else 'nonCenter',
                         'p'+str(p).replace('.',''),
                         str(pd.Timestamp.today()).split()[0].replace('-','')))

path='/oak/stanford/groups/mrivas/projects/degas-risk/datasets/all_pop/'
data.to_pickle(path + dataset_name + '.full_df.pkl.gz', compression='gzip')
