
## HLA imputation for R9

# libs
library(HIBAG)
library(data.table)
library(tidyverse)
library(parallel)

# functions
source('./src/VCF_functions.R')

# computing nodes
cl <- makeCluster(4)


## Extract MHC SNPs needed in imputation

system(paste0("plink --bfile /finngen/library-red/finngen_R9/genotype_plink_1.0/data/finngen_R9 ", 
              "--allow-no-sex --extract /home/ivm/Documents/models/hla.snp.list.2 --make-bed --out data/genotypes/R9_MHC_HLA"))


## Read MHC genotype data 

genotypes <- hlaBED2Geno(bed.fn='/home/ivm/Documents/R9/data/genotypes/R9_MHC_HLA.bed', 
                         bim.fn='/home/ivm/Documents/R9/data/genotypes/R9_MHC_HLA.bim', 
                         fam.fn='/home/ivm/Documents/R9/data/genotypes/R9_MHC_HLA.fam', 
                         assembly='hg38')

## Finnish hg38 models

fin38model.A <- hlaModelFromObj(get(load('/home/ivm/Documents/models/Fin_hg38_model_A.RData'))[['A']])
imputed.a    <- predict(fin38model.A, genotypes, type='response+prob', match.type='Position', cl=cl)
write.table(imputed.a$value, 'results/R9_A_imputed.tsv', quote=F, sep='\t', row.names=F)
saveRDS(imputed.a, file='results/R9_A_imputed.rds')

fin38model.B <- hlaModelFromObj(get(load('/home/ivm/Documents/models/Fin_hg38_model_B.RData'))[['B']])
imputed.b    <- predict(fin38model.B, genotypes, type='response+prob', match.type='Position', cl=cl)
write.table(imputed.b$value, 'results/R9_B_imputed.tsv', quote=F, sep='\t', row.names=F)
saveRDS(imputed.b, file='results/R9_B_imputed.rds')

fin38model.C <- hlaModelFromObj(get(load('/home/ivm/Documents/models/Fin_hg38_model_C.RData'))[['C']])
imputed.c    <- predict(fin38model.C, genotypes, type='response+prob', match.type='Position', cl=cl)
write.table(imputed.c$value, 'results/R9_C_imputed.tsv', quote=F, sep='\t', row.names=F)
saveRDS(imputed.c, file='results/R9_C_imputed.rds')

fin38model.DPB1 <- hlaModelFromObj(get(load('/home/ivm/Documents/models/Fin_hg38_model_DPB1.RData'))[['DPB1']])
imputed.dpb1    <- predict(fin38model.DPB1, genotypes, type='response+prob', match.type='Position', cl=cl)
write.table(imputed.dpb1$value, 'results/R9_DPB1_imputed.tsv', quote=F, sep='\t', row.names=F)
saveRDS(imputed.dpb1, file='results/R9_DPB1_imputed.rds')

fin38model.DRB1 <- hlaModelFromObj(get(load('/home/ivm/Documents/models/Fin_hg38_model_DRB1.RData'))[['DRB1']])
imputed.drb1    <- predict(fin38model.DRB1, genotypes, type='response+prob', match.type='Position', cl=cl)
write.table(imputed.drb1$value, 'results/R9_DRB1_imputed.tsv', quote=F, sep='\t', row.names=F)
saveRDS(imputed.drb1, file='results/R9_DRB1_imputed.rds')

fin38model.DQA1 <- hlaModelFromObj(get(load('/home/ivm/Documents/models/Fin_hg38_model_DQA1.RData'))[['DQA1']])
imputed.dqa1    <- predict(fin38model.DQA1, genotypes, type='response+prob', match.type='Position', cl=cl)
write.table(imputed.dqa1$value, 'results/R9_DQA1_imputed.tsv', quote=F, sep='\t', row.names=F)
saveRDS(imputed.dqa1, file='results/R9_DQA1_imputed.rds')

fin38model.DQB1 <- hlaModelFromObj(get(load('/home/ivm/Documents/models/Fin_hg38_model_DQB1.RData'))[['DQB1']])
imputed.dqb1    <- predict(fin38model.DQB1, genotypes, type='response+prob', match.type='Position', cl=cl)
write.table(imputed.dqb1$value, 'results/R9_DQB1_imputed.tsv', quote=F, sep='\t', row.names=F)
saveRDS(imputed.dqb1, file='results/R9_DQB1_imputed.rds')

fin38model.DRB3 <- hlaModelFromObj(get(load('/home/ivm/Documents/models/DRB3_model.RData'))[['DRB3']])
imputed.drb3    <- predict(fin38model.DRB3, genotypes, type='response+prob', match.type='Position', cl=cl)
imputed.drb3$locus <- 'DRB3'
write.table(imputed.drb3$value, 'results/R9_DRB3_imputed.tsv', quote=F, sep='\t', row.names=F)
saveRDS(imputed.drb3, file='results/R9_DRB3_imputed.rds')

fin38model.DRB4 <- hlaModelFromObj(get(load('/home/ivm/Documents/models/DRB4_model.RData'))[['DRB4']])
imputed.drb4    <- predict(fin38model.DRB4, genotypes, type='response+prob', match.type='Position', cl=cl)
imputed.drb4$locus <- 'DRB4'
write.table(imputed.drb4$value, 'results/R9_DRB4_imputed.tsv', quote=F, sep='\t', row.names=F)
saveRDS(imputed.drb4, file='results/R9_DRB4_imputed.rds')

fin38model.DRB5 <- hlaModelFromObj(get(load('/home/ivm/Documents/models/DRB5_model.RData'))[['DRB5']])
imputed.drb5    <- predict(fin38model.DRB5, genotypes, type='response+prob', match.type='Position', cl=cl)
imputed.drb5$locus <- 'DRB5'
write.table(imputed.drb5$value, 'results/R9_DRB5_imputed.tsv', quote=F, sep='\t', row.names=F)
saveRDS(imputed.drb5, file='results/R9_DRB5_imputed.rds')


## convert imputed HLAs to VCF

# imputed.a    <- readRDS('results/R9_A_imputed.rds')
# imputed.b    <- readRDS('results/R9_B_imputed.rds')
# imputed.c    <- readRDS('results/R9_C_imputed.rds')
# imputed.dpb1 <- readRDS('results/R9_DPB1_imputed.rds')
# imputed.dqa1 <- readRDS('results/R9_DQA1_imputed.rds')
# imputed.dqb1 <- readRDS('results/R9_DQB1_imputed.rds')
# imputed.drb1 <- readRDS('results/R9_DRB1_imputed.rds')
# imputed.drb3 <- readRDS('results/R9_DRB3_imputed.rds')
# imputed.drb4 <- readRDS('results/R9_DRB4_imputed.rds')
# imputed.drb5 <- readRDS('results/R9_DRB5_imputed.rds')

imputation2vcf(imputed.a$postprob, 'A', 'results/R9_A_imputed.vcf')
imputation2vcf(imputed.b$postprob, 'B', 'results/R9_B_imputed.vcf')
imputation2vcf(imputed.c$postprob, 'C', 'results/R9_C_imputed.vcf')
imputation2vcf(imputed.dpb1$postprob, 'DPB1', 'results/R9_DPB1_imputed.vcf')
imputation2vcf(imputed.dqa1$postprob, 'DQA1', 'results/R9_DQA1_imputed.vcf')
imputation2vcf(imputed.dqb1$postprob, 'DQB1', 'results/R9_DQB1_imputed.vcf')
imputation2vcf(imputed.drb1$postprob, 'DRB1', 'results/R9_DRB1_imputed.vcf')
imputation2vcf(imputed.drb3$postprob, 'DRB3', 'results/R9_DRB3_imputed.vcf')
imputation2vcf(imputed.drb4$postprob, 'DRB4', 'results/R9_DRB4_imputed.vcf')
imputation2vcf(imputed.drb5$postprob, 'DRB5', 'results/R9_DRB5_imputed.vcf')

# collect into a single VCF file, compress
system("./src/R9_vcf_finalize.sh")



