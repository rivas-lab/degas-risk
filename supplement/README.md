## Supplementary Figures and Results

Like the [figures](../figures) folder, this folder contains jupyter notebooks which generated the supplementary figures for the manuscript. As the same caveat with respect to individual-level UK Biobank data holds here, it is similarly not possible to use these notebooks as-is until other analyses are re-run locally. They are provided as-is for replication and future analysis. 

### Notes

With few exceptions, there is a one-to-one correspondence between notebooks and supplementary figures, tables, and data. Exceptions are listed below:

- Data for Table S1 are derived from the notebook for Figure S1
- The "Figure S1b" notebook results are not present in the manuscript version of the supplement
- Figure S3 was derived using the phenotype pipeline from the [scripts](../scripts) directory
- Visualizations for Data S2 and S3 (phenotype and gene contribution scores) were respectively made using the notebook for Figure S1 and `gene_plots.py`. Numeric values for these contribution scores are available for download on the [Global Biobank Engine](biobankengine.stanford.edu/downloads)

### Model Performance

Results from profiling the predictive performance of polygenic risk scores (PRS), including the DeGAs polygenic risk score (dPRS), are in the `results` folder. Each file contains results from either PRS or dPRS in one of the four analysis populations: train, dev, and test sets of white British individuals, and an additional test set of non-British white individuals (`nbw`) from UK Biobank. 

Columns in these files are as follow:
- `pheno`: Global Engine Phenotype code (see [reference](../reference) folder for more information).
- `n`: Number of cases, or individuals with quantitative measurement
- `beta2`: Regression coefficient for being in the top 2% of PRS or dPRS for a given trait (adjusted for covariates)
- `auc`: Area under the receiver operating curve using PRS or dPRS (unadjusted) as the classifying score, for binary traits
- `pearsonr`: Pearson's r between PRS and dPRS and quantitative trait values
- `spearmanr`: Spearman's rho (rank correlation) between PRS and dPRS and trait values 

Additional information on these measures, including covariate adjustment, can be found in the Methods section of our manuscript.

### Supplementary Results

We also computed DeGAs on an additional publicly available [dataset](http://www.nealelab.is/blog/2017/7/19/rapid-gwas-of-thousands-of-phenotypes-for-337000-samples-in-the-uk-biobank) of summary statistics from UK Biobank GWAS. Please refer to the documentation in the `extra` folder for more details.
