## DeGAs on Neale Lab Rapid GWAS

Here are scripts used to compute DeGAs on version 1 of summary statistics from the [Neale Lab Rapid GWAS] of UK Biobank data. These data were downloaded, and subject to the following analysis steps:
- Filtering for LD-independence (`prune.sh`)
- Compiling into an accessible format for python (`make_neale_dataset.py`)
- Computing DeGAs (`neale_tsvd.py`).

We have not extensively considered the results of this application of DeGAs, though some initial analysis is in Figure S8 (see [notebook](../supplement/FigureS8-Neale-DeGAs-Instance.ipynb) in the [parent directory](../supplement)). Likewise, we have not benchmarked the performance of PRS or dPRS using these data -- though please feel free to do so! This part of the repository is mainly indended to serve as an additional example application of DeGAs on public data, with starter code for future analysis. Happy hunting!
