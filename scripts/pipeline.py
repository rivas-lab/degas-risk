#!/bin/python
import os
import numpy as np
import pandas as pd
from statsmodels.api import OLS,Logit
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import argparse
from scipy.stats import zscore,pearsonr,spearmanr
from scipy.spatial.distance import euclidean
from sklearn.metrics import roc_curve,roc_auc_score
from sklearn.preprocessing import normalize
from sklearn.cluster import KMeans

# 1. load data
dataset=''
npc=0
phe_codes=[]
covariates=['age','sex']+['PC'+str(i+1) for i in range(4)]+['1']
score_pcs=['SCORE{}_SUM'.format(pc+1) for pc in range(npc)]
profl_pcs=['PROF_PC{}'.format(pc+1) for pc in range(npc)]
z=np.load(dataset)
scores=pd.read_table('/oak/stanford/groups/mrivas/projects/degas-risk/scorefiles/'+
                      os.path.splitext(os.path.basename(dataset))[0]+'.sscore',
                     index_col='#IID')
phenos=pd.read_table('/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.phe',
                     usecols=['IID']+phe_codes+covariates[:-1],
                     index_col='IID',
                     na_values=-9)
phenos['1']=1
prs=reduce(lambda x,y: pd.merge(x,y,left_index=True,right_index=True),
           map(lambda phe: pd.read_table(os.path.join('/oak/stanford/groups/mrivas/projects/degas-risk/PRS/',
                                             os.path.splitext(os.path.basename(dataset))[0][:-7], 
                                             phe+'_PRS.profile'),
                                         usecols=['IID','SCORESUM'], 
                                         index_col='IID').rename({'SCORESUM':phe+'_PRS'}, axis=1),
               phe_codes))
with open('/oak/stanford/groups/mrivas/users/magu/repos/rivas-lab/ukbb-tools/05_gbe/icdinfo.txt','r') as f:
    code_to_name = {line.split()[0]:line.split()[2].replace('_',' ').split('(')[0].capitalize() for line in f}

# 2. define population groupings
train=set(pd.read_table('/oak/stanford/groups/mrivas/projects/degas-risk/population-split/'+
                         'ukb24983_white_british_train.phe').iloc[:,0].astype(float).tolist())
valid=set(pd.read_table('/oak/stanford/groups/mrivas/projects/degas-risk/population-split/'+
                         'ukb24983_white_british_valid.phe').iloc[:,0].astype(float).tolist())
test=set(pd.read_table('/oak/stanford/groups/mrivas/projects/degas-risk/population-split/'+
                         'ukb24983_white_british_test.phe').iloc[:,0].astype(float).tolist())
nbw=set(pd.read_table('/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/'+
                         'ukb24983_non_british_white.phe').iloc[:,0].astype(float).tolist())

# subset data, package for use
df=phenos.loc[train | valid | test | nbw, phe_codes+covariates].dropna()
df=df.merge(scores.merge(prs, left_index=True, right_index=True), left_index=True, right_index=True)

# 3. analysis
for phe_code in phe_codes:
    # setup: names, model weights/scores, trait data type
    phe_name=code_to_name[phe_code]
    weights=z['V'][np.where(z['label_phe_code'] == phe_code),:].flatten()[:npc]
    df['dPRS']=df[score_pcs].dot(weights)
    df.sort_values(by='dPRS', inplace=True, ascending=True)
    is_bin=(len(df[phe_code].value_counts()) == 2)
    if is_bin:
        # adjust from plink 1/2 coding 
        df[phe_code]-=1
        regress=Logit
    else:
        regress=OLS
    # define relevant models and performance summary statistics
    models={m:{} for m in ['dPRS','COVAR','RESID','JOINT']}
    params={'dPRS':['dPRS'], 'COVAR':covariates, 'JOINT':covariates+['dPRS']}
    stats =pd.DataFrame(index=['train','valid','test','nbw'], 
                        columns=['n','beta2','auc','pearsonr','spearmanr'])
    # iterate over populations
    for pop_id,pop in zip(['train','valid','test','nbw'],[train,valid,test,nbw]):
        # normalize dPRS and (if qt) trait values
        df2=df.loc[pop,:].dropna()
        df2['dPRS']=zscore(df2['dPRS'])
        if not is_bin:
            df2[phe_code]=zscore(df2[phe_code])
        pop=[ind for ind in pop if ind in df2.index]
        # iterate over models, compute statistics from relevant groupings 
        # joint model for train/valid fit in train; for test/nbw fit in valid
        for m in ['dPRS','COVAR','RESID','JOINT']:
            if m=='RESID':
                models[m][pop_id]=OLS(models['COVAR'][pop_id].resid_pearson, 
                                      df2.loc[pop,'dPRS']).fit(disp=0)
            elif m=='JOINT' and pop_id in ['test','nbw']:
                models[m][pop_id]=models[m]['valid']
            elif m=='JOINT' and pop_id == 'valid':
                models[m][pop_id]=models[m]['train']
            else:
                models[m][pop_id]=regress(df2.loc[pop,phe_code], 
                                          df2.loc[pop,params[m]]).fit(disp=0)
        # compute statistics from relevant grouping
        # beta2 is OR/beta of 2%ile against population; auc uses joint model (bin) 
        # pearsonr is on residuals (qt); spearmanr is dPRS only
        if is_bin:
            stats[pop_id,'beta2']=df2[phe_code].iloc[-int(0.02*len(pop)):].mean() / df2[phe_code].mean()
            stats[pop_id,'auc']=roc_auc_score(df2.loc[pop,phe_code],
                                              models['JOINT'][pop_id].predict(df2.loc[pop,:]))
            stats[pop_id,'pearsonr']='na'
            stats[pop_id,'n']=df2.loc[pop,phe_code].count(1)
        else:
            stats[pop_id,'beta2']=df2[phe_code].iloc[-int(0.02*len(pop)):].mean() - df2[phe_code].mean()
            stats[pop_id,'auc']='na'
            stats[pop_id,'pearsonr']=pearsonr(models['COVAR'][pop_id].resid,
                                              df2.loc[pop,'dPRS'])[0]
            stats[pop_id,'n']=df2.loc[pop,phe_code].shape[0]
        stats[pop_id,'spearmanr']=spearmanr(df2.loc[pop,phe_code],df2.loc[pop,'dPRS'])[0]
        # paint profiles and make figures for these two populations
        if pop_id == 'test' or pop_id == 'nbw':
            # profiles
            for pc in range(npc):
                df2[profl_pcs[pc]]=(df2[score_pcs[pc]] * weights[pc] * df['dPRS']).clip_lower(0)
            df2[profl_pcs] = normalize(df2[profl_pcs], norm='l1')
            # median risk profile (for clustering) and outlier identification
            centroid=df2[profl_pcs].median()
            df3=df2.iloc[-int(0.1*len(pop)):]
            df3['mahal']=df3[profl_pcs].apply(lambda x: euclidean(x, centroid), axis=1)
            m_star=df3['mahal'].mean() + 2*df3['mahal'].std()
            outliers=df3.query('mahal > @m_star').index
            # perform clustering
            k=1
            cluster = KMeans(n_clusters=k, n_init=25).fit(df3.loc[outliers,profl_pcs])
            pre_frac = 0.8 - 10.0/len(outliers)
            errors = [cluster.inertia_]
            for new_k in [2,3,4,5]:
                new_cluster = KMeans(n_clusters=new_k, n_init=25).fit(df3.loc[outliers,profl_pcs])
                errors.append(new_cluster.inertia_)
                if new_cluster.inertia_ / cluster.inertia_ < pre_frac**(new_k - k): 
                    cluster,k = new_cluster,new_k

