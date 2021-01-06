## Reference data

This folder contains files used for reference in the project — specifically, variant filtration and phenotype selection. Variant filtration is performed by a script in the [scripts folder](../scripts). Traits were selected by filtering on sample size from genome-wide association studies used for the [Global Biobank Engine](biobankengine.stanford.edu), and were then subject to manual curation. Both processes are described in detail in the manuscript.

### Files

- Variants included: `variant_qc_v2.prune.in`
- Traits included: `final_sumstats_v3.txt`
- Traits excluded: `blacklist.txt`
- Traits in final DeGAs instance: `final_phe_codes_v2.txt`

The above files for traits often reference traits by their Global Biobank Engine IDs, which are not always human readable. For example, `HC382` corresponds to [Asthma](https://biobankengine.stanford.edu/RIVAS_HG19/coding/HC382), which was defined using both hospital records and verbal questionnaire data as part of a [digital phenotyping project](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7212271/). A map between GBE IDs and names can be downloaded [here](https://biobankengine.stanford.edu/downloads/).

NB: Due to phenotype inclusion criteria, not all analysis traits present in `final_sumstats_v3.txt` were present in the final DeGAs instance (`final_phe_codes_v2.txt`), due to lack of genetic associations at `p<1e-6` and to manual curation. 
