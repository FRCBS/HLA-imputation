
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

# set finngen data release
FG.release <- 'R10'

# path to SNPs to be extracted
hla.snp.list <- '/finngen/shared/hla.snp.list.2/20220808_115746/files/ritarij/hla.snp.list.2' 



## Extract MHC SNPs needed in imputation

system(paste0("plink --bfile /finngen/library-red/finngen_", FG.release, "/genotype_plink_1.0/data/finngen_", FG.release, 
              " --allow-no-sex --extract ", hla.snp.list, " --make-bed --out data/genotypes/", FG.release, "_MHC_HLA"))



## Read MHC genotype data (may need the full path here)

genotypes <- hlaBED2Geno(bed.fn=paste0('data/genotypes/', FG.release, '_MHC_HLA.bed'), 
                         bim.fn=paste0('data/genotypes/', FG.release, '_MHC_HLA.bim'), 
                         fam.fn=paste0('data/genotypes/', FG.release, '_MHC_HLA.fam'), 
                         assembly='hg38')



## impute using Finnish hg38 models

fin38model.A <- hlaModelFromObj(get(load('models/hg38/Fin_hg38_model_A.RData'))[['A']])
imputed.a    <- predict(fin38model.A, genotypes, type='response+prob', match.type='Position', cl=cl)
write.table(imputed.a$value, paste0('results/', FG.release, '_A_imputed.tsv'), quote=F, sep='\t', row.names=F)
saveRDS(imputed.a, file=paste0('results/', FG.release, '_A_imputed.rds'))

fin38model.B <- hlaModelFromObj(get(load('models/hg38/Fin_hg38_model_B.RData'))[['B']])
imputed.b    <- predict(fin38model.B, genotypes, type='response+prob', match.type='Position', cl=cl)
write.table(imputed.b$value, paste0('results/', FG.release, '_B_imputed.tsv'), quote=F, sep='\t', row.names=F)
saveRDS(imputed.b, file=paste0('results/', FG.release, '_B_imputed.rds'))

fin38model.C <- hlaModelFromObj(get(load('models/hg38/Fin_hg38_model_C.RData'))[['C']])
imputed.c    <- predict(fin38model.C, genotypes, type='response+prob', match.type='Position', cl=cl)
write.table(imputed.c$value, paste0('results/', FG.release, '_C_imputed.tsv'), quote=F, sep='\t', row.names=F)
saveRDS(imputed.c, file=paste0('results/', FG.release, '_C_imputed.rds'))

fin38model.DPB1 <- hlaModelFromObj(get(load('models/hg38/Fin_hg38_model_DPB1.RData'))[['DPB1']])
imputed.dpb1    <- predict(fin38model.DPB1, genotypes, type='response+prob', match.type='Position', cl=cl)
write.table(imputed.dpb1$value, paste0('results/', FG.release, '_DPB1_imputed.tsv'), quote=F, sep='\t', row.names=F)
saveRDS(imputed.dpb1, file=paste0('results/', FG.release, '_DPB1_imputed.rds'))

fin38model.DRB1 <- hlaModelFromObj(get(load('models/hg38/Fin_hg38_model_DRB1.RData'))[['DRB1']])
imputed.drb1    <- predict(fin38model.DRB1, genotypes, type='response+prob', match.type='Position', cl=cl)
write.table(imputed.drb1$value, paste0('results/', FG.release, '_DRB1_imputed.tsv'), quote=F, sep='\t', row.names=F)
saveRDS(imputed.drb1, file=paste0('results/', FG.release, '_DRB1_imputed.rds'))

fin38model.DQA1 <- hlaModelFromObj(get(load('models/hg38/Fin_hg38_model_DQA1.RData'))[['DQA1']])
imputed.dqa1    <- predict(fin38model.DQA1, genotypes, type='response+prob', match.type='Position', cl=cl)
write.table(imputed.dqa1$value, paste0('results/', FG.release, '_DQA1_imputed.tsv'), quote=F, sep='\t', row.names=F)
saveRDS(imputed.dqa1, file=paste0('results/', FG.release, '_DQA1_imputed.rds'))

fin38model.DQB1 <- hlaModelFromObj(get(load('models/hg38/Fin_hg38_model_DQB1.RData'))[['DQB1']])
imputed.dqb1    <- predict(fin38model.DQB1, genotypes, type='response+prob', match.type='Position', cl=cl)
write.table(imputed.dqb1$value, paste0('results/', FG.release, '_DQB1_imputed.tsv'), quote=F, sep='\t', row.names=F)
saveRDS(imputed.dqb1, file=paste0('results/', FG.release, '_DQB1_imputed.rds'))

fin38model.DRB3 <- hlaModelFromObj(get(load('models/hg38/DRB3_model.RData'))[['DRB3']])
imputed.drb3    <- predict(fin38model.DRB3, genotypes, type='response+prob', match.type='Position', cl=cl)
imputed.drb3$locus <- 'DRB3'
write.table(imputed.drb3$value, paste0('results/', FG.release, '_DRB3_imputed.tsv'), quote=F, sep='\t', row.names=F)
saveRDS(imputed.drb3, file=paste0('results/', FG.release, '_DRB3_imputed.rds'))

fin38model.DRB4 <- hlaModelFromObj(get(load('models/hg38/DRB4_model.RData'))[['DRB4']])
imputed.drb4    <- predict(fin38model.DRB4, genotypes, type='response+prob', match.type='Position', cl=cl)
imputed.drb4$locus <- 'DRB4'
write.table(imputed.drb4$value, paste0('results/', FG.release, '_DRB4_imputed.tsv'), quote=F, sep='\t', row.names=F)
saveRDS(imputed.drb4, file=paste0('results/', FG.release, '_DRB4_imputed.rds'))

fin38model.DRB5 <- hlaModelFromObj(get(load('models/hg38/DRB5_model.RData'))[['DRB5']])
imputed.drb5    <- predict(fin38model.DRB5, genotypes, type='response+prob', match.type='Position', cl=cl)
imputed.drb5$locus <- 'DRB5'
write.table(imputed.drb5$value, paste0('results/', FG.release, '_DRB5_imputed.tsv'), quote=F, sep='\t', row.names=F)
saveRDS(imputed.drb5, file=paste0('results/', FG.release, '_DRB5_imputed.rds'))



## convert imputed HLAs to BGEN

# load imputation result data

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

imputation2bgen(imputed.a$postprob,    'A',    paste0('results/', FG.release,'_HLA_A'))
imputation2bgen(imputed.b$postprob,    'B',    paste0('results/', FG.release,'_HLA_B'))
imputation2bgen(imputed.c$postprob,    'C',    paste0('results/', FG.release,'_HLA_C'))
imputation2bgen(imputed.dpb1$postprob, 'DPB1', paste0('results/', FG.release,'_HLA_DPB1'))
imputation2bgen(imputed.dqa1$postprob, 'DQA1', paste0('results/', FG.release,'_HLA_DQA1'))
imputation2bgen(imputed.dqb1$postprob, 'DQB1', paste0('results/', FG.release,'_HLA_DQB1'))
imputation2bgen(imputed.drb1$postprob, 'DRB1', paste0('results/', FG.release,'_HLA_DRB1'))
imputation2bgen(imputed.drb3$postprob, 'DRB3', paste0('results/', FG.release,'_HLA_DRB3'))
imputation2bgen(imputed.drb4$postprob, 'DRB4', paste0('results/', FG.release,'_HLA_DRB4'))
imputation2bgen(imputed.drb5$postprob, 'DRB5', paste0('results/', FG.release,'_HLA_DRB5'))


