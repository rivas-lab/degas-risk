#!/bin/python
import numpy as np
import pandas as pd
import sys,os
from statsmodels.formula.api import logit,ols
from statsmodels.api import OLS
from scipy.stats import pearsonr,spearmanr

# ensure usage, process input
if len(sys.argv) < 3:
    print("usage: python validate.py [file.npz] [file.sscore] <npc> <phenos>")
    sys.exit(1)
else:
    # process
    npz_f, pc_f = sys.argv[1:3]
    # load zipped output from tsvd
    npz=np.load(npz_f)
    # excludes NMISS_ALLELE_CT and DENOM
    pcs=pd.read_table(pc_f, index_col='#IID')
    # if no phenotypes are provided, use everything in the reference
    if len(sys.argv) > 3:
        if sys.argv[3].isdigit():
            npc=int(sys.argv[3])
            if len(sys.argv) > 4:
                phe=[i for i in sys.argv[4:] if i in npz['label_phe_code']]
            else:
                phe=list(npz['label_phe_code'])
        else:
            npc=500
            phe=[i for i in sys.argv[3:] if i in npz['label_phe_code']]
    else:
        phe=list(npz['label_phe_code']) 

# load data: train and test populations and phenotype values
train=pd.read_table('/oak/stanford/groups/mrivas/projects/degas-risk/'+
                    'population-split/ukb24983_white_british_train.phe',
                    usecols=[0]).values.flatten().tolist()
# evaluation is in a test set (valid.phe) of 10% of WBR
test=pd.read_table('/oak/stanford/groups/mrivas/projects/degas-risk/'+
                    'population-split/ukb24983_white_british_valid.phe',
                    usecols=[0]).values.flatten().tolist()
phenos=pd.read_table('/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/'+
                     'master_phe/master.phe', 
                     usecols=['IID','age','sex','PC1','PC2','PC3','PC4']+phe,
                     index_col='IID', na_values=-9)
# remove duplicate IID; adjust binary traits from 2/1 to 1/0
phenos=phenos.loc[~phenos.index.duplicated(keep='first')]
phenos[[p for p in phe if 'INI' not in p and 'QT' not in p]] -= 1

# combine data and redefine cohorts 
data=pd.merge(phenos, pcs, left_index=True, right_index=True)
train=list(filter(lambda i: i in data.index, train))
test=list(filter(lambda i: i in data.index, test))

# set up output file, initialize scores, metrics, and null model formula
m='spearman'
correlate=lambda x,y: spearmanr(x,y)[0] if m=='spearman' else pearsonr(x,y)[0]
out=os.path.join(os.path.dirname(pc_f), 'results', 
           os.path.splitext(os.path.basename(pc_f.replace('500PC',
                                                          str(npc)+'PC')))[0]+
           '.'+m+'r.tsv')
data['SCORE']=0
cols=[x+'_'+y for x in ['TRAIN','TEST'] for y in ['RAW','RESID','JOINT','COVAR']]

# compute validation for all phenotypes
with open(out, 'w') as o:
    o.write('\t'.join(['PHE']+cols)+'\n')
    for phe_id in phe:
        model=ols if 'INI' in phe_id or 'QT' in phe_id else logit
        w=npz['V'][np.where(npz['label_phe_code']==phe_id),:].flatten()[:npc]
        data['SCORE']=data.iloc[:,-501:-(501-npc)].dot(w)
        for ix,pop in enumerate([train,test]):
            pop=data.loc[pop,[phe_id,'SCORE']].dropna().index.tolist()
            r1=data.loc[pop,['SCORE',phe_id]].corr(method=m).iloc[0,1]
            try:
                nm=model(formula=phe_id+'~ age+sex+PC1+PC2+PC3+PC4', 
                            data=data.loc[pop,:]).fit()
                jm=model(formula=phe_id+'~ age+sex+PC1+PC2+PC3+PC4+SCORE', 
                            data=data.loc[pop,:]).fit()
                r2=correlate(nm.resid_pearson,data.loc[pop,'SCORE'])
                r3=correlate(jm.fittedvalues,data.loc[pop,'SCORE'])
                r4=correlate(nm.fittedvalues,data.loc[pop,'SCORE'])
            except:
                r2,r3,r4='NA','NA','NA'
            if ix == 0:
                o.write('\t'.join(map(str,[phe_id,r1,r2,r3,r4,""])))
            if ix == 1:
                o.write('\t'.join(map(str,[r1,r2,r3,r4]))+'\n')
