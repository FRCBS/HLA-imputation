# Finnish HLA HIBAG imputation models

Repository for code, summary statistic data and trained HLA models


Publication:

Ritari J, Hyvärinen K, Clancy J, FinnGen, Partanen J, Koskela S. _Increasing accuracy of HLA imputation by a population-specific reference panel in a Finngen biobank cohort._ NAR Genomics and Bioinformatics, Volume 2, Issue 2, June 2020, lqaa030, https://doi.org/10.1093/nargab/lqaa030 


## code (./src)
`results.R` R code that generates the plots for figures 2-5

`functions.R` helper functions for e.g. plotting

`HLA_imputation.R` example of running imputation with Finnish reference panel

`FG_HLA_impute.R` R code for running HLA imputation for FinnGen R9, and result VCF file production


## reference panels (./models)
Trained HLA imputation models for hg19 and hg38 human genome builds.
The models for DRB3-5 genes are in hg38 build only. The 'ng' in DRB3-5 is not an allele as such, but  indicates that the gene is missing.

## summary statistic data (./data)
Contains imputation error rates, autoimmune association summaries and HLA allele frequency information

## FinnGen HLA imputation methods description

MHC region SNPs were first converted to plink format (.bed, .bim, .fam) from VCF genotypes using [tabix](https://www.htslib.org/doc/tabix.html) and [plink](https://www.cog-genomics.org/plink/), and subequently processed with [R](https://www.r-project.org/). HLA alleles for HLA-A, -B, -C, -DPB1, -DQA1, -DQB1, -DRB1, and -DRB3-5 genes were imputed with the [HIBAG](https://bioconductor.org/packages/release/bioc/html/HIBAG.html) R library using a Finnish reference panel. Imputation posterior probabilities (pp) for each imputed allele were extracted from the HIBAG output by summing the pp values from all imputed HLA genotypes of which a given allele was a part of. The pp values were imported to plink2 using `--import-dosage` command and converted to BGEN format with `--export bgen-1.2 bits=16 ref-first`.

