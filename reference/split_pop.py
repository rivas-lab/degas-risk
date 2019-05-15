#!/bin/python
import random

# init
random.seed(24983)

# get samples
in_pop="/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe"
with open(in_pop, 'r') as f:
    fid_iid = [line.rstrip().split() for line in f]

# make splits
random.shuffle(fid_iid)
n1 = int(0.6*len(fid_iid))
n2 = int(0.2*len(fid_iid))

# write to files
outD="/oak/stanford/groups/mrivas/projects/degas-risk/population-split/"

with open(outD + "ukb24983_white_british_train.phe", "w") as o:
    o.write("\n".join(map(lambda l:"\t".join(l), fid_iid[:n1])))

with open(outD + "ukb24983_white_british_test.phe", "w") as o:
    o.write("\n".join(map(lambda l:"\t".join(l), fid_iid[n1:n1+n2])))

with open(outD + "ukb24983_white_british_valid.phe", "w") as o:
    o.write("\n".join(map(lambda l:"\t".join(l), fid_iid[n1+n2:])))
