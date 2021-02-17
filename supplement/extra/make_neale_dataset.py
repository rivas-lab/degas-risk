#!/bin/python3
import glob
import numpy as np
import pandas as pd

collab_dir=""
proj_dir=""

sumstats=glob.glob(collab_dir+'associations/*.gz')
variants=pd.read_table(proj_dir+'neale/ld/neale_ld_indep.dual_id.prune.in', header=None).iloc[:,1]
filemaps=pd.read_table(collab_dir+'scripts/rename/ukb1189_1859_phenosummary_final.HC.lookup.tsv')
filemaps.iloc[:,1]=filemaps.iloc[:,1].astype(str)
print(filemaps.iloc[:,1].head())

data={}
for f in sumstats:
	print(f)
	pheno_t=f.split('/')[-1].split('.')[0]
	if pheno_t not in filemaps.iloc[:,1].values or filemaps[filemaps.iloc[:,1].values==pheno_t].iloc[0,2] < 100:
		continue
	x=pd.read_table(f)
	subset=x.query('(pval < 1e-6) & (se < 0.08) & (rsid in @variants)')
	subset.index=subset['rsid']	
	data[pheno_t]=subset['beta'].rename(pheno_t)	

pd.DataFrame(data).to_pickle('neale_beta_p1e-6.ld-indep.unfiltered-phenos.sumstats.pkl.gz')
