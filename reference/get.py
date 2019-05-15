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

# hard code biomarker fields and ids
biomarker_ids = ["INI20030620", "INI20030600", "INI10030610", "INI20030630",
                 "INI20030640", "INI20030650", "INI20030680", "INI20030690",
                 "INI20030710", "INI20030700", "INI20030720", "INI20030660",
                 "INI20030730", "INI10030740", "INI20030750", "INI20030760",
                 "INI20030770", "INI20030780", "INI20030790", "INI20030810",
                 "INI20030830", "INI20030850", "INI10030840", "INI20030860",
                 "INI20030870", "INI10030880", "INI20030670", "INI20030890",
                 "INI10030510", "INI40030700", "INI10030500", "INI10030520",
                 "INI10030530", "INI10030910", "BIN10030510", "BIN10030500",
                 "BIN10030820", "BIN10030800"]

biomarker_fields = ["30620", "30600", "30610", "30630", "30640", "30650", 
                    "30680", "30690", "30710", "30700", "30720", "30660", 
                    "30730", "30740", "30750", "30760", "30770", "30780", 
                    "30790", "30810", "30830", "30850", "30840", "30860", 
                    "30870", "30880", "30670", "30890", "30510", "30700", 
                    "30500", "30520", "30530", "30910", "30510", "30500", 
                    "30820", "30800"]

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
        # remove acceleration phenotypes
        if "acceleration" in name:
            continue 
        # filter biomarkers by the above tables
        if gbe_id[-5:] in biomarker_fields and gbe_id not in biomarker_ids:
            continue
        # remove phenotypes which are angles
        if "angle" in name:
            continue
        # remove phenotypes which mark data acquisition methods
        if "method" in name:
            continue
        # remove environmental phenotypes (pollution/traffic/etc.)
        if gbe_id in ["INI"+str(24000+i) for i in range(3,25)]:
            continue
        # remove acquisition phase, safety, and completion info for brain MRIs
        if gbe_id in ["INI25780","BIN12139","BIN12188"]:
            continue
        # finally
        o.write("\t".join([gbe_id,n,name,path+'\n']))

