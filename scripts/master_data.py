#!/bin/python3
import sys
import numpy as np
import pandas as pd
import joblib
# import sklearn.externals.joblib as joblib

# ensure usage, pick statistic to keep, ensure usage
usage="usage: python master_data.py [z OR b OR p OR se]"
if len(sys.argv) < 2:
    print(usage); exit(1)

z,b,p,se = tuple(map(lambda x:sys.argv[1].lower() == x, ['z','b','p','se']))
if z + b + p + se != 1:
    print(usage); exit(1)
else:
    print("Extracting {} values...".format("z" if z else "beta" if b else "p")) 

# files to import be here
with open('../reference/final_sumstats_v3.tsv', 'r') as f:
	file_ref = [(line.split()[0],line.rstrip().split()[1:]) for line in f]

# only keep variants not in LD, MAF > 0.01%, and QC'd
with open('../reference/variant_qc_v2.prune.in', 'r') as f:
	var_set = set([line.rstrip() for line in f])

# get alt alleles -- need to standardize summary stats ALT allele
ref_bim = '/'.join(('/oak/stanford/groups/mrivas/ukbb24983',
                    'array_combined/pgen/ukb24983_cal_hla_cnv.pvar'))
var_alt = {var:None for var in var_set}
with open(ref_bim, 'r') as f:
	for line in f:
		chrom,pos,var,ref,alt = line.rstrip().split()
		var_alt[var] = alt	

# helper function to load files
def get_summary_stats(files, z=True, b=False, p=False, se=False):
	# identify file type by path
	binary,qt = 'logistic' in files[0], 'linear' in files[0]
	# load data, perform p-value and LD pruning
	x = pd.concat([pd.read_table(f, index_col=2) for f in files]).fillna(0)
	x = x[(x['TEST'] == "ADD") & (x.index.isin(var_set))]
	x['SIGN'] = 2*x['ALT'].eq(x['A1']) - 1 # 1 if ALT is A1 else -1
	# take LOR if binary, manage z-scoring
	if qt:
		stat = 'T_STAT' if z else 'BETA' if b else 'P' if p else 'SE'
		return x[stat].astype(float).divide(x['SIGN'] if not (p or se) else 1)
	elif binary:
		stat = 'Z_STAT' if z else 'OR' if b else 'P' if p else 'SE'
		func = lambda x: np.log(float(x)) if z else float(x)
		return x[stat].apply(func).divide(x['SIGN'] if not (p or se) else 1)
	else:
		exit(1)

# load data one phenotype at a time
data = pd.DataFrame(index=var_set)
for phe_id,files in file_ref:
	data[phe_id] = get_summary_stats(files, z=z, b=b, p=p, se=se)

# name and save
dataset_name = '_'.join(('all',
                         'z' if z else 'beta' if b else 'p' if p else 'se',
                         'nonCenter',
                         str(pd.Timestamp.today()).split()[0].replace('-','')))

path='/oak/stanford/groups/mrivas/projects/degas-risk/datasets/train/v2/'
joblib.dump(data, path + dataset_name + '.full_df.pkl.gz', compress=5)
