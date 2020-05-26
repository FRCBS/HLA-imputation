# Finnish HLA HIBAG imputation models

Repository for R code, summary statistic data and trained HLA models

 
Ritari J, Hyvärinen K, Clancy J, FinnGen, Partanen J, Koskela S. _Increasing accuracy of HLA imputation by a population-specific reference panel in a Finngen biobank cohort._ NAR Genomics and Bioinformatics, Volume 2, Issue 2, June 2020, lqaa030, https://doi.org/10.1093/nargab/lqaa030 


## code (./src)
`results.R` R code that generates the plots for figures 2-5

`functions.R` helper functions for e.g. plotting

`HLA_imputation.R` example of running imputation with Finnish reference panel


## reference panels (./models)
Trained HLA imputation models for hg19 and hg38 human genome builds

## summary statistic data (./data)
Contains imputation error rates, autoimmune association summaries and HLA allele frequency information
