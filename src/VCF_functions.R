
## libs
library(data.table)
library(tidyverse)

## opts
options(stringsAsFactors=F, scipen=100, digits=4)


## -----------------------------------------
## functions
## -----------------------------------------

# function to convert allele dosage table into VCF
# allele dosage is coded into the GP field
imputation2vcf <- function(i.data, hla.gene, outfile) {
  
  # i.data    = HIBAG imputation post prob table
  # hla.gene  = HLA gene
  # outfile   = output VCF path/filename
  
  # HLA gene start positions in hg38
  hla.genepos <- data.frame(gene=c('A', 'B', 'C', 'DPB1' ,'DQA1', 'DQB1', 'DRB1', 'DRB3', 'DRB4', 'DRB5'), 
                            start=c(29.94126, 31.26949, 31.26875, 33.07599, 
                                    32.62818, 32.65947, 32.57877, 32.5, 32.5, 32.5)*1e6) %>% 
    filter(gene==hla.gene)
  
  
  # imputed alleles
  out.alleles      <- rownames(i.data) %>% str_split_fixed(., '/', 2) 
  out.alleles.list <- out.alleles %>% as.vector() %>% unique()
  out.alleles.list <- out.alleles.list[nchar(out.alleles.list)>4]
  
  # collect imputation pp values per allele
  out.alleles.pp <- map_dfr(1:length(out.alleles.list), function(x) {
    i.data[c(which(out.alleles[, 1] %in% out.alleles.list[x]), which(out.alleles[, 2] %in% out.alleles.list[x])), ] %>% 
      colSums
  }) %>% data.frame
  
  # plink dosage format
  out <- data.frame(SNP=out.alleles.list,
                    A1='<present>',
                    A2='<absent>',
                    out.alleles.pp)
  out$SNP <- paste0(hla.gene, '*', out$SNP)
  
  # fam file
  tmp.map <- data.frame(colnames(out)[-c(1:3)], colnames(out)[-c(1:3)], 0, 0, 0, 0)
  
  # sample IDs to dosage file
  colnames(out)[-c(1:3)] <- paste(colnames(out)[-c(1:3)], colnames(out)[-c(1:3)]) 
  
  # write to tmp
  fwrite(out, 'tmp/tmp.dosage', sep='\t')
  fwrite(tmp.map, 'tmp/tmp.dosage.map', sep='\t', col.names=F)
  
  # coversion to VCF
  system(paste0("plink2 --import-dosage tmp/tmp.dosage single-chr=6 ref-last ",
                "--fam tmp/tmp.dosage.map --export vcf vcf-dosage=GP --out tmp/tmp"))
  
  # save VCF header
  system("head -n 6 ./tmp/tmp.vcf > tmp/tmp.header")
  
  # read the tmp VCF
  tmpvcf <- fread('./tmp/tmp.vcf', data.table=F)
  # fix sample names
  colnames(tmpvcf)[10:ncol(tmpvcf)] <- 
    str_split_fixed(colnames(tmpvcf)[10:ncol(tmpvcf)], '_', 2)[, 1]
  # insert position
  tmpvcf$POS <- hla.genepos$start
  
  # write output
  fwrite(tmpvcf, './tmp/tmp.vcf', sep='\t')
  system(paste("cat tmp/tmp.header tmp/tmp.vcf >", outfile))
  
}

# alleles into one column
formatImpTab <- function(x) {
  return(
    rbind(
      data.frame(sample.id=x$sample.id, allele=x$allele1, prob=x$prob),
      data.frame(sample.id=x$sample.id, allele=x$allele2, prob=x$prob)
    )
  )
}


