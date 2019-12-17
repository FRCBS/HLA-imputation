

library(HIBAG)


## --------------------------------------------------
## Impute with trained HIBAG models
## --------------------------------------------------

## Load Plink formatted genotype data

hlageno <- hlaBED2Geno(bed.fn="./data/MHC.bed", 
                       fam.fn="./data/MHC.fam", 
                       bim.fn="./data/.bim", assembly='hg19')

## Load models

# Finnish hg19/GRCh37 models
fin19.model.a    <- hlaModelFromObj(get(load('./models/hg19/Fin_hg19_model_A.RData'))[['A']])
fin19.model.b    <- hlaModelFromObj(get(load('./models/hg19/Fin_hg19_model_B.RData'))[['B']])
fin19.model.c    <- hlaModelFromObj(get(load('./models/hg19/Fin_hg19_model_C.RData'))[['C']])
fin19.model.drb1 <- hlaModelFromObj(get(load('./models/hg19/Fin_hg19_model_DRB1.RData'))[['DRB1']])
fin19.model.dqa1 <- hlaModelFromObj(get(load('./models/hg19/Fin_hg19_model_DQA1.RData'))[['DQA1']])
fin19.model.dqb1 <- hlaModelFromObj(get(load('./models/hg19/Fin_hg19_model_DQB1.RData'))[['DQB1']])
fin19.model.dpb1 <- hlaModelFromObj(get(load('./models/hg19/Fin_hg19_model_DPB1.RData'))[['DPB1']])

# default European
model.list <- get(load('./models/hg19/European-HLA4-hg19.RData'))
model.a    <- hlaModelFromObj(model.list[['A']])
model.b    <- hlaModelFromObj(model.list[['B']])
model.c    <- hlaModelFromObj(model.list[['C']])
model.drb1 <- hlaModelFromObj(model.list[['DRB1']])
model.dqa1 <- hlaModelFromObj(model.list[['DQA1']])
model.dqb1 <- hlaModelFromObj(model.list[['DQB1']])
model.dpb1 <- hlaModelFromObj(model.list[['DPB1']])


## Imputation

# Fin 
fin19.hlageno.a    <- predict(fin19.model.a,    hlaGenoSubset(hlageno), type="response+prob", match.type="Position")
fin19.hlageno.b    <- predict(fin19.model.b,    hlaGenoSubset(hlageno), type="response+prob", match.type="Position")
fin19.hlageno.c    <- predict(fin19.model.c,    hlaGenoSubset(hlageno), type="response+prob", match.type="Position")
fin19.hlageno.drb1 <- predict(fin19.model.drb1, hlaGenoSubset(hlageno), type="response+prob", match.type="Position")
fin19.hlageno.dqa1 <- predict(fin19.model.dqa1, hlaGenoSubset(hlageno), type="response+prob", match.type="Position")
fin19.hlageno.dqb1 <- predict(fin19.model.dqb1, hlaGenoSubset(hlageno), type="response+prob", match.type="Position")
fin19.hlageno.dpb1 <- predict(fin19.model.dpb1, hlaGenoSubset(hlageno), type="response+prob", match.type="Position")

# default European
eur.hlageno.a    <- predict(model.a,    hlaGenoSubset(hlageno), type="response+prob", match.type="Position")
eur.hlageno.b    <- predict(model.b,    hlaGenoSubset(hlageno), type="response+prob", match.type="Position")
eur.hlageno.c    <- predict(model.c,    hlaGenoSubset(hlageno), type="response+prob", match.type="Position")
eur.hlageno.drb1 <- predict(model.drb1, hlaGenoSubset(hlageno), type="response+prob", match.type="Position")
eur.hlageno.dqa1 <- predict(model.dqa1, hlaGenoSubset(hlageno), type="response+prob", match.type="Position")
eur.hlageno.dqb1 <- predict(model.dqb1, hlaGenoSubset(hlageno), type="response+prob", match.type="Position")
eur.hlageno.dpb1 <- predict(model.dpb1, hlaGenoSubset(hlageno), type="response+prob", match.type="Position")


