#!/bin/python
import sys
import os
import numpy as np
import pandas as pd
from sklearn.externals import joblib

_README_="""
A script to subset pickled dataframes (output from master_dataset.py) for use in DEGAS (tsvd.py).
(you'll probably want at least 64Gb memory for this)

Author: Matthew Aguirre (SUNET: magu)
"""

def load_p_and_se(dataset):
    rename=lambda s,k:s.replace('_z_','_'+k+'_').replace('_beta_','_'+k+'_')
    pf,sef=rename(dataset,'p'),rename(dataset,'se')
    return joblib.load(pf), joblib.load(sef)


def make_dataset(dataset, out, p_star=0.01, center=True):
    # keep variants not in LD, MAF > 0.01%, and QC'd; get their alt alleles
    with open('../reference/variant_qc.prune.in', 'r') as f:
        var_set = set([line.rstrip() for line in f])
    var_alt = {var:None for var in var_set}
    with open('/oak/stanford/groups/mrivas/ukbb24983/array_combined/pgen/ukb24983_cal_hla_cnv.pvar', 'r') as f:
        for line in f:
            chrom,pos,var,ref,alt = line.rstrip().split()
            var_alt[var] = alt    
    # load data
    data=joblib.load(dataset)
    p,se=load_p_and_se(dataset)
    # filter by p-value and SE (0.2 if binary else 0.08); center if specified
    qts=[x for x in se.columns if 'INI' in x or 'QT' in x]
    df=data[(p < p_star) & (se < 0.2) & ~((se.isin(se[qts])) & (se >= 0.08))]
    if center:
        df=df.subtract(df.mean()).divide(data.std())
    # sparsify and write to file
    df.to_sparse().to_pickle(out, compression='gzip')
    return


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=_README_
    )
    def p_value(p):
        if float(p) < 0 or float(p) > 1:
            raise argparse.ArgumentTypeError("%s is an invalid p-value" % p)
        return float(p)    
    parser.add_argument('--p', dest="p", required=True, nargs=1, type=p_value,
                            default=0.01,
                            help='Minimum p-value threshold for inclusion.')
    parser.add_argument('--center', dest="center", action='store_true',
                            help='Standardize columns (phenotypes).') 
    bz=parser.add_mutually_exclusive_group()
    bz.add_argument('--beta', dest="beta", action='store_false',
                        help='Use regression coefficients for dataset.') 
    bz.add_argument('--z', dest="z", action='store_true',
                        help='Use regression z-statistics for dataset.')
    args=parser.parse_args()
    c,p=args.center,float(args.p[0])
    # input/output naming
    path=os.path.join('/oak/stanford/groups/mrivas/projects/degas-risk/',
                      'datasets/all_pop',
                      '_'.join(('all','z' if args.z else 'beta',
                                'nonCenter_20190621.full_df.pkl.gz')))
    outP=os.path.join(os.path.dirname(path), 
                      '_'.join(('all','z' if args.z else 'beta',
                                'center' if c else 'nonCenter',
                                'p'+str(p).replace('.',''), 
                                '20190621.full_df.pkl.gz')))
    # everything else goes here
    make_dataset(dataset=path, out=outP, p_star=p, center=c)

