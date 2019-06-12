#!/bin/python
import sys
import os
import numpy as np
import pandas as pd

# ensure usage
if len(sys.argv) < 3:
    print("usage: python " + sys.argv[0] + " dataset.pkl.gz PHE_ID")
    sys.exit(1)

# input
dataset = sys.argv[1]
pheno = sys.argv[2]

# work here
out_dir = os.path.join(os.path.dirname(dataset), 'PRS')
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# load data
bfile="/oak/stanford/groups/mrivas/ukbb24983/array_combined/pgen/ukb24983_cal_hla_cnv"
bim = pd.read_table(bfile+'.pvar', index_col='ID')
coefs = pd.read_pickle(dataset)
coefs = coefs.join(bim['ALT'])

# extract coefficients
weight_file = os.path.join(out_dir, pheno+'_weights.txt')
coefs[['ALT',pheno]].dropna().to_dense().to_csv(weight_file, sep='\t', header=None)

# score alleles with plink
os.system("ml load biology; ml load plink")
os.system(" ".join(["plink", "--bfile", bfile,
                             "--score", weight_file, "sum center double-dosage",
                             "--out", os.path.join(out_dir, pheno+"_PRS")]))
