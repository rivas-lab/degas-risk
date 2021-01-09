## Supplementary Figures

Like the [figures](../figures) folder, this folder contains jupyter notebooks which generated the supplementary figures for the manuscript. As the same caveat with respect to individual-level UK Biobank data holds here, it is similarly not possible to use these notebooks as-is until other analyses are re-run locally. They are provided as-is for replication and future analysis. 

### Notes

With few exceptions, there is a one-to-one correspondence between notebooks and supplementary figures, tables, and data. Exceptions are listed below:

- Data for Table S1 are derived from the notebook for Figure S1
- The "Figure S1b" notebook results are not present in the manuscript version of the supplement, and are here provided as bonus results
- Figure S3 was derived using the phenotype pipeline from the [scripts](../scripts) directory
- Visualizations for Data S2 and S3 (phenotype and gene contribution scores) were respectively made using `phe_plots.py` and `gene_plots.py`. Numeric values for these contribution scores are available for download on the [Global Biobank Engine](biobankengine.stanford.edu/downloads)
