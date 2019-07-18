#!/bin/python
import numpy as np
import pandas as pd
import sys,os
from statsmodels.formula.api import ols,Logit

# ensure usage, process input
if len(sys.argv) < 3:
    print("usage: python validate.py [dataset.npz] [dataset.sscore] <phenos>")
    sys.exit(1)
else:
    # load zipped output from tsvd
    npz=np.load(sys.argv[1])
    # excludes NMISS_ALLELE_CT and DENOM
    pcs=pd.read_table(sys.argv[2], index_col='#IID')
    # if no phenotypes are provided, use everything in the reference
    if len(sys.argv) > 3:
        phe=sys.argv[3:]
    else:
        phe=list(npz['label_phe_code']) 

# load data: train (wbr) and test (nbw) populations and phenotype values
wbr=pd.read_table('/oak/stanford/groups/mrivas/ukbb24983/sqc/population_'+
                  'stratification/ukb24983_white_british.phe',
                   header=None).iloc[:,0].tolist()
nbw=pd.read_table('/oak/stanford/groups/mrivas/ukbb24983/sqc/population_'+
                  'stratification/ukb24983_non_british_white.phe',
                   header=None).iloc[:,0].tolist()
phenos=pd.read_table('/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/'+
                     'master_phe/master.phe', 
                     usecols=['IID','age','sex','PC1','PC2','PC3','PC4']+phe,
                     index_col='IID', na_values=-9)

# combine data and redefine cohorts 
data=pd.merge(phenos, pcs, left_index=True, right_index=True)
wbr=list(filter(lambda i: i in data.index, wbr))
nbw=list(filter(lambda i: i in data.index, nbw))

# set up output file, initialize scores, metrics, and null model formula
out=os.path.join(os.path.dirname(sys.argv[2]), 'results', 
           os.path.splitext(os.path.basename(sys.argv[2]))[0]+'.pearsonr.tsv')
data['SCORE']=0
cols=[x+'_'+y for x in ['WBR','NBW'] for y in ['RAW','RESID','JOINT','COVAR']]
covs='age+sex+PC1+PC2+PC3+PC4'

# compute validation for all phenotypes
with open(out, 'w') as o:
    o.write('\t'.join(['PHE']+cols)+'\n')
    for phe_id in phe:
        o.write(phe_id+'\t')
        model=ols if 'INI' in phe_id or 'QT' in phe_id else Logit
        weights=npz['V'][np.where(npz['label_phe_code'] == phe_id),:].flatten()
        data['SCORE']=data.iloc[:,-501:-1].dot(weights)
        for pop in [wbr,nbw]:
            r1=data.loc[pop,['SCORE',phe_id]].corr().iloc[0,1]
            nm=model(formula=phe_id+'~'+covs, data=data.loc[pop]).fit()
            r2=ols(nm.pearson_resid, data.loc[pop,'SCORE']).fit().rsquared_adj
            r3=model(formula=phe_id+'~'+covs+'+SCORE', 
                        data=data.loc[pop]).fit().rsquared_adj
            r4=nm.rsquared_adj
            o.write('\t'.join(map(str,[r1,r2,r3,r4])))
        o.write('\n')
