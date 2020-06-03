#!/bin/python3
import numpy as np 
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import gzip
import sys
import random

# parameters
npc=500

# dataset and variant:gene info
z=np.load('/oak/stanford/groups/mrivas/projects/degas-risk/datasets/train/v2/tsvd/all_beta_center_p1e-06_20200506_500PCs.npz', allow_pickle=True)
var2gene={}
with gzip.open('/oak/stanford/groups/mrivas/private_data/ukbb/variant_filtering/variant_filter_table.new.tsv.gz', 'r') as f:
    for line in f:
        arr=line.decode().rstrip().split('\t')
        if len(arr) > 30:
            var2gene[arr[4]]=arr[30]

# mapping
# var2gene={v:random.choice(g.split(',')) if g and g != '""' else v for v,g in var2gene.items()}
var2gene={v:random.choice(g.split(',')) if g else '""' for v,g in var2gene.items()}
varcont=pd.DataFrame(z['U'][:,:npc], columns=['PC{}'.format(i+1) for i in range(npc)])
# varcont.index=pd.MultiIndex.from_tuples([(var2gene[v] if v in var2gene else v,v) for v in z['label_var']])
varcont.index=pd.MultiIndex.from_tuples([(var2gene[v] if v in var2gene else '""',v) for v in z['label_var']])


# loop over pcs
patterns = ('-', '+', 'x', '\\', '*', '.')
for pc in range(npc):
    # compute contribution scores
    gene2cont={}
    for i in range(varcont.shape[0]):
        if varcont.index[i][0] not in gene2cont:
            gene2cont[varcont.index[i][0]] = 0
        gene2cont[varcont.index[i][0]] += varcont.iloc[i,pc]**2
    # draw the figure
    plt.figure(figsize=(5,5), dpi=300)
    base=0
    for ix,(gene,c2) in enumerate(sorted(gene2cont.items(),key=lambda x:int(x[0]=='""')-x[1])):
        if c2 > 1e-3 and ix < 30:
            plt.bar(list(range(3)), [c2,0,0], bottom=[base,0,0], label=gene,
                    hatch=patterns[(ix + 4) % len(patterns)])
        else:
            plt.bar(list(range(3)), [1-base,0,0], bottom=base,  color='grey')
            plt.xticks(list(range(3)), ['PC{}'.format(pc+1),'',''])
            plt.ylabel('Gene contribution score')
            plt.legend(fontsize=6, loc='upper right')#, bbox_to_anchor=(1, 1))
            plt.title('Gene contribution to PC'+str(pc+1))
            break
        base+=c2
    plt.savefig('pdfs/gene_v2/gene_pc'+str(pc+1)+'.pdf')
    plt.close()
    print(pc)
