## Scripts for DeGAs and (d)PRS

This directory contains the scripts used to run primary analyses for this project. These include genome-wide association scans (GWAS) using UK Biobank data, data processing, DeGAs, and polygenic risk score building and evaluation. Each of these steps (and its corresponding scripts) is detailed below.

### Genome-wide associations

UK Biobank individuals of white British ancestry are split at random (70\% train, 10\% validation, 20\% test) for all analyses, using `split_pop.sh`. Variant quality control and GWAS are both conducted in the train set, respectively using `variant_qc.sh` and `gwas.sh`. Phenotypes were manually curated -- please see the [reference](../reference) folder for more information.

### Data Assembly

As summary stats from \~1000 GWAS are rather cumbersome to work with, two scripts are used to pre-process output from GWAS using PLINK, from the previous step. Complete matrixes of all summary stats are assembled using `master_data.py`, set to query either betas, z-scores, standard errors, or p-values. These full data matrixes are then subsetted using `make_dataset.py` to generate input for matrix decomposition with DeGAs.

### DeGAs and (d)PRS

Decomposition of genetic associations (DeGAs) is run using `tsvd.py`. Note that DeGAs hyperparameters, other than the number of components to compute, are specified by the input matrix. This script also optionally computes component polygenic scores (cPRS), which can also be computed separately using `score_pc.sh`.

Polygenic risk scores (PRS) using DeGAs input data is performed using `PRS.py`. Helper functions in this script are also used to compute the DeGAs polygenic risk score, along with other analyses using DeGAs (e.g. risk profile clustering), for each trait. These analyses are run using `pipeline.py`, which has an accompanying bash script (`batch_pipeline.sh`) to assist parallelization, and jupyter notebook (`test_phe_pipeline.ipynb`) for testing. Performance measures for PRS and dPRS were computed using `validate.py`. 
