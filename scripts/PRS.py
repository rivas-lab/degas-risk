#!/bin/python
import sys
import os
import errno
import numpy as np
import pandas as pd

_README_="""
A script to compute vanilla (threshold and clump) from a dataset that would be
input to DEGAS (tsvd.py). No clumping is necessary since the input file only
contains variants which are in LE with one another (pruned).

Author: Matthew Aguirre (SUNET: magu)
"""

# function used in pipeline
def PRS(dataset, pheno):
    # work here
    splitext = lambda f: os.path.splitext(f)[0] # convenient
    out_dir = splitext(splitext(splitext(dataset))).replace('datasets','PRS')
    try:
        os.makedirs(out_dir)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(out_dir):
            pass
        else:
            raise

    # load data
    bfile="/oak/stanford/groups/mrivas/ukbb24983/array_combined/pgen/ukb24983_cal_hla_cnv"
    bim = pd.read_table(bfile+'.pvar', index_col='ID')
    coefs = pd.read_pickle(dataset)
    coefs = coefs.join(bim['ALT'])

    # extract coefficients
    weight_file = os.path.join(out_dir, pheno+'_weights.txt')
    coefs[['ALT',pheno]].dropna().to_csv(weight_file, sep='\t', header=None)

    # score alleles with plink
    m_phe="/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.phe"
    os.system("ml load biology; ml load plink")
    os.system(" ".join(["plink", "--bfile", bfile,
                                 "--score", weight_file, "sum center double-dosage",
                                 "--pheno", m_phe, "--pheno-name", pheno,
                                 "--memory 8000",
                                 "--out", os.path.join(out_dir, pheno+"_PRS")]))

# 
if __name__ == '__main__':
    # ensure usage
    if len(sys.argv) < 3:
        print("usage: python " + sys.argv[0] + " dataset.pkl.gz PHE_ID1 ...")
        sys.exit(1)
    # do the thing
    for phe in sys.argv[2:]:
        PRS(sys.argv[1], phe)
