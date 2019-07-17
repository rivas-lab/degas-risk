#!/bin/python
import numpy as np
import pandas as pd
import sys,os

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
                     na_values=-9, usecols=['IID']+phe, index_col='IID')

# combine data and redefine cohorts 
data=pd.merge(phenos, pcs, left_index=True, right_index=True)
wbr=list(filter(lambda i: i in data.index, wbr))
nbw=list(filter(lambda i: i in data.index, nbw))

# compute validation for all phenotypes
out=os.path.join(os.path.dirname(sys.argv[2]), 'results', 
                os.path.splitext(os.path.basename(sys.argv[2]))[0]+'.pearsonr.tsv')
data['SCORE']=0
with open(out, 'w') as o:
    for phe_id in phe:
        weights=npz['V'][np.where(npz['label_phe_code'] == phe_id),:].flatten()
        data['SCORE']=data.iloc[:,-501:-1].dot(weights)
        r1=data.loc[wbr,['SCORE',phe_id]].corr(method='pearson').iloc[0,1]
        r2=data.loc[nbw,['SCORE',phe_id]].corr(method='pearson').iloc[0,1]
        o.write('\t'.join([phe_id,str(r1),str(r2)])+'\n')
