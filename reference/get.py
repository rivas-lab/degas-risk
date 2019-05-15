#!/bin/python

# get phenotype info
ids,names,info = set(),set(),{}
with open("../../ukbb-tools/05_gbe/phenotype_info.tsv", "r") as i:
    for n,line in enumerate(i):
        # skip header
        if n > 0:
            gbe_id,name,_,_,_,_,_,n_gbe,_,_,_,_,_,path = line.rstrip().split()
            names.add(name)
            ids.add(gbe_id)
            info[gbe_id] = (name,n_gbe,path)

# write subset to file
with open("phenotypes.tsv", "w") as o:
    # write header
    o.write("\t".join(["#GBE_ID","GBE_N","GBE_NAME","PATH\n"]))
    # write data
    for gbe_id,(name,n,path) in sorted(info.items(), key=lambda x:x[0]):
        # filter on number of observations/cases
        if int(n) < 100:
            continue
        # remove pilot phenotypes
        if "(pilot)" in name:
            continue
        # remove medications
        if "MED" in gbe_id:
            continue
        # remove family history and rohit's phenotypes
        if "RH" in gbe_id or "FH" in gbe_id:
            continue
        # remove acceleration phenotypes
        if "acceleration" in name:
            continue 
        # remove unadjusted biomarkers
        if name + "_covariate_and_statin_adjusted" in names:
            continue
        # remove biomarkers not adjusted for statins
        if "_covariate_adjutsed" in name:
            continue
        # remove phenotypes which are angles
        if "angle" in name:
            continue
        # remove acquisition phase for MRI images
        if gbe_id == "INI25780":
            continue
        o.write("\t".join([gbe_id,n,name,path+'\n']))

