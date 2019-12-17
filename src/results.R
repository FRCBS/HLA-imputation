
## Libs
library(tidyverse)
library(data.table)
library(HIBAG)
library(ggrepel)
library(ggpubr)

theme_set(theme_grey())

## functions
source('./src/Functions.R')


## --------------------------------------------------
## Imputation reference error analysis
## Figure 2
## --------------------------------------------------

# errors and post probs
reg37.hlageno.r.err <- fread('./data/reg37.hlageno.r.err.tsv', data.table=F) 
reg37.hlageno.d.err <- fread('./data/reg37.hlageno.d.err.tsv', data.table=F) 
reg38.hlageno.r.err <- fread('./data/reg38.hlageno.r.err.tsv', data.table=F)
reg38.hlageno.d.err <- fread('./data/reg38.hlageno.d.err.tsv', data.table=F) 
def.hlageno.r.err   <- fread('./data/def.hlageno.r.err.tsv', data.table=F) 
def.hlageno.d.err   <- fread('./data/def.hlageno.d.err.tsv', data.table=F) 

# CV for selection of imputation ref by best imputation probability score 
imputed.cv <- lapply(1:100, function(iter) {
  
  sample.ind <- sample(1:nrow(reg37.hlageno.d.err), nrow(reg37.hlageno.d.err), T)

  out <- rbind(sapply(list(1:2, 3:4, 5:6, 7:8, 9:10, 11:12, 13:14), correctErrs02D, wvec=c(1,1,1), inds=sample.ind),
                sapply(list(1:2, 3:4, 5:6, 7:8, 9:10, 11:12, 13:14), correctErrs02R, wvec=c(1,1,1), inds=sample.ind))
  
  colnames(out) <- c('A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPB1')
  rownames(out) <- c('All_D', 'Fin38_D', 'Fin37_D', 'Eur37_D', 'All_R', 'Fin38_R', 'Fin37_R', 'Eur37_R')
  out.list <- list(out)
  names(out.list) <- c('Bootstrap')
  out.list
  
})

# CV errors
imputed.cv.df <- lapply(1:100, function(i) imputed.cv[[i]] %>% melt %>% filter(., L1=='Bootstrap') %>% 
                          data.frame(., Fold=i)) %>% do.call(rbind, .) 
imputed.cv.df$Var1 <- gsub('_R|_D', '', imputed.cv.df$Var1)
imputed.cv.df <- imputed.cv.df %>% group_by(Var1, Var2, L1, Fold) %>% summarize(value=mean(value)) 

imputed.cv.df$Var1[(imputed.cv.df$L1=='Test_imputed' & imputed.cv.df$Var1=='All')] <- 'W'
imputed.cv.df$Var1[imputed.cv.df$Var1=='All'] <- 'bestPP'
imputed.cv.df <- imputed.cv.df[!(imputed.cv.df$L1=='Test_imputed' & (imputed.cv.df$Var1!='W' & imputed.cv.df$Var1!='Eq')), ]
imputed.cv.df$Var1 <- imputed.cv.df$Var1 %>% factor(., levels=c('bestPP', 'Fin37', 'Fin38', 'Eur37'), ordered=F)
colnames(imputed.cv.df)[1] <- 'Method'

# plot figure 2
pdf('./results/Fig2_Imputation_errors_by_method.pdf', width=6, height=3.5)
ggplot(imputed.cv.df %>% filter(Method!='W') %>% data.frame %>% 
         mutate(., value=value/212), aes(Var2, value*100, fill=Method)) +
  geom_boxplot(outlier.size=0.5, outlier.shape=1, size=0.15, alpha=0.6, color='black', 
               position=position_dodge2(preserve='single')) +
  xlab('HLA gene') +
  ylab('% errors') +
  labs(color='Reference') +
  scale_fill_viridis_d(option='D', name='reference') +
  theme(panel.background=element_blank(), axis.line=element_line(colour="black", size=.2), legend.key=element_blank())
dev.off()



## --------------------------------------------------
## Imputation reference error analysis
## Figure 3
## --------------------------------------------------


# table of error and probs per gene and reference
error.prob.all <- rbind(
  pmap(list(seq(1, 14, by=2), c('A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPB1')), function(b, i) {
    data.frame(rbind(reg37.hlageno.d.err, reg37.hlageno.r.err)[, c(b, b+1)], Gene=i) }) %>% do.call(rbind, .) %>% 
    data.frame(., Ref='Fin37'),
  pmap(list(seq(1, 14, by=2), c('A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPB1')), function(b, i) {
    data.frame(rbind(reg38.hlageno.d.err, reg38.hlageno.r.err)[, c(b, b+1)], Gene=i) }) %>% do.call(rbind, .) %>% 
    data.frame(., Ref='Fin38'),
  pmap(list(seq(1, 14, by=2), c('A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPB1')), function(b, i) {
    data.frame(rbind(def.hlageno.d.err, def.hlageno.r.err)[, c(b, b+1)], Gene=i) }) %>% do.call(rbind, .)  %>% 
    data.frame(., Ref='Eur37'))

# calculate erros over probability cutoff values
# HIBAG PP > i --> imputation errors
error.prob.all.cutoffs <- data.frame( # Fin37 ref
  CUTOFF=seq(0, 1, 0.01),
  ERROR_PROB=sapply(seq(0, 1, 0.01), function(i) {
    tmp <- error.prob.all %>% filter(Ref=='Fin37', PROB>=i)
    sum(tmp$ERRS) / nrow(tmp) }),
  RETAIN=sapply(seq(0, 1, 0.01), function(i) {
    tmp <- error.prob.all %>% filter(Ref=='Fin37', PROB>=i)
    nrow(tmp) / error.prob.all %>% filter(Ref=='Fin37') %>% nrow })
)
error.prob.all.cutoffs.2 <- data.frame( # Eur37 ref
  CUTOFF=seq(0, 1, 0.01),
  ERROR_PROB=sapply(seq(0, 1, 0.01), function(i) {
    tmp <- error.prob.all %>% filter(Ref=='Eur37', PROB>=i)
    sum(tmp$ERRS) / nrow(tmp) }),
  RETAIN=sapply(seq(0, 1, 0.01), function(i) {
    tmp <- error.prob.all %>% filter(Ref=='Eur37', PROB>=i)
    nrow(tmp) / error.prob.all %>% filter(Ref=='Eur37') %>% nrow })
)

# add linear models
error.prob.all.cutoffs <- data.frame(
  error.prob.all.cutoffs,
  ERROR_PROB_PRED=glm(ERROR_PROB ~ I(1/(1+CUTOFF^2)) + I(exp(CUTOFF)), 
                      data=error.prob.all.cutoffs) %>% predict(., error.prob.all.cutoffs),
  RETAIN_PRED=smooth.spline(error.prob.all.cutoffs$CUTOFF, error.prob.all.cutoffs$RETAIN, spar=0.7) %>% .$y)

error.prob.all.cutoffs.2 <- data.frame(
  error.prob.all.cutoffs.2,
  ERROR_PROB_PRED=glm(ERROR_PROB ~ I(1/(1+CUTOFF^2)) + I(exp(CUTOFF)), 
                      data=error.prob.all.cutoffs.2) %>% predict(., error.prob.all.cutoffs.2),
  RETAIN_PRED=smooth.spline(error.prob.all.cutoffs.2$CUTOFF, error.prob.all.cutoffs.2$RETAIN, spar=0.7) %>% .$y )

colnames(error.prob.all.cutoffs.2) <- paste0('EUR_', colnames(error.prob.all.cutoffs.2))

error.prob.all.cutoffs  <- cbind(error.prob.all.cutoffs, error.prob.all.cutoffs.2)


# difference between Fin and Eur error proportions over cutoff thresholds
sum(error.prob.all.cutoffs$ERROR_PROB - error.prob.all.cutoffs$EUR_ERROR_PROB, na.rm=T)

# permutate refrence labels and compute difference
error.prob.all.cutoffs.perm <- sapply(1:10000, function(i) {
  tmp <- error.prob.all.cutoffs
  ind <- sample(1:nrow(tmp)/2, nrow(tmp)/2, F)
  tmp$ERROR_PROB[ind] <- tmp$EUR_ERROR_PROB[ind]
  sum(tmp$ERROR_PROB*100 - tmp$EUR_ERROR_PROB*100, na.rm=T)
})

# permutate labels and compute total number of errors between Fin and Eur refs
error.prob.perm <- sapply(1:10000, function(i) {
  tmp <- error.prob.all
  tmp$ERRS <- error.prob.all[base::sample(1:nrow(error.prob.all), nrow(error.prob.all), F), 'ERRS']
  sum(tmp %>% filter(Ref=='Fin37') %>% .$ERRS) - sum(tmp %>% filter(Ref=='Eur37') %>% .$ERRS)
})


# Plots

# plot error rate over prob thresholds for Fin and Eur refs
p.error.prob.all.cutoffs <- ggplot(error.prob.all.cutoffs %>% 
                                     select(., CUTOFF, Finnish=ERROR_PROB_PRED, European=EUR_ERROR_PROB_PRED) %>% 
                                     melt(., c('CUTOFF')),
                                   aes(CUTOFF, value*100, linetype=variable, color=variable)) +
  geom_line(size=0.5) +
  geom_point(x=c(error.prob.all.cutoffs$CUTOFF, error.prob.all.cutoffs$CUTOFF), 
             y=c(error.prob.all.cutoffs$ERROR_PROB*100, error.prob.all.cutoffs$EUR_ERROR_PROB*100),
             shape=21, color='grey50', size=0.8, alpha=0.8) +
  xlab('posterior probability cutoff') + ylab('% errors') +
  scale_color_manual(values=c('grey55', 'black')) +
  theme(panel.background=element_blank(), axis.line=element_line(colour="black", size=.2)) +
  theme(legend.position='top', legend.title=element_blank(), legend.text=element_text(size=11), legend.key=element_blank())

# plot proportion of samples left over prob thresholds
p.error.prob.all.retained <- ggplot(error.prob.all.cutoffs %>% 
                                      select(., CUTOFF, Finnish=RETAIN_PRED, European=EUR_RETAIN_PRED) %>% 
                                      melt(., c('CUTOFF')),
                                    aes(CUTOFF, value*100, linetype=variable, color=variable)) +
  geom_line(size=0.5) +
  xlab('posterior probability cutoff') + ylab('% samples left') +
  scale_color_manual(values=c('grey55', 'black')) +
  theme(panel.background=element_blank(), axis.line=element_line(colour="black", size=.2)) +
  theme(legend.position='top', legend.title=element_blank(), legend.text=element_text(size=11), legend.key=element_blank())


# plot total error difference permutation test result histogram
p.error.prob.perm <- ggplot(data.frame(PERMS=error.prob.perm), aes(PERMS)) +
  geom_histogram(alpha=0.8, bins=25, fill='grey60') + 
  xlab('difference in number of errors') +
  ggtitle('permutation test Fin vs. Eur\n(p = 0.015)') +
  geom_vline(xintercept=-35, color='red', alpha=0.7) +
  theme(panel.background=element_blank(), axis.line=element_line(colour="black", size=.2)) +
  theme(legend.position='none', legend.title=element_blank(), plot.title=element_text(size=11))


# plot permutation value histogram and actual value as red vertical line
p.error.prob.all.cutoffs.perm <- ggplot(data.frame(PERMS=error.prob.all.cutoffs.perm), aes(PERMS)) +
  geom_histogram(alpha=0.8, bins=55, fill='grey60') + 
  xlab('difference in sum of error %') +
  ggtitle('permutation test Fin vs. Eur\n(p = 1e-4)') +
  geom_vline(xintercept=-95.37677, color='red', alpha=0.7) +
  theme(panel.background=element_blank(), axis.line=element_line(colour="black", size=.2)) +
  theme(legend.position='none', legend.title=element_blank(), plot.title=element_text(size=11))


# plot PP value vs. error number for Fin and Eur refs
p.err.ref <- ggplot(error.prob.all %>% filter(Ref!='Fin38'), 
                    aes(y=PROB, x=ERRS %>% factor)) +
  geom_quasirandom(alpha=0.7, varwidth=T, shape=19, size=.45, color='grey50') + 
  geom_boxplot(alpha=0.5, varwidth=F, outlier.shape=NA, size=0.2, color='black') +
  xlab('number of errors per locus') + ylab('posterior probability') +
  stat_summary(fun.data=function(x) return(c(label=length(x), y=1.05)), 
               geom="text", position=position_dodge(width=0.75)) +
  facet_wrap(~Ref, labeller=labeller(Ref=c(Fin37='Finnish\n', Eur37='European\n'))) +
  theme(panel.background=element_blank(), axis.line=element_line(colour="black", size=.2)) +
  theme(strip.text.x=element_text(size=11), panel.grid.major.x=element_blank(), strip.background=element_blank(), 
        legend.position='none')


# compute and plot AUC/ROC for FIN and Eur references 
# 0: 0 errors (controls), 1: >=1 erros (cases)
p.err.roc <- ggroc(list(`Finnish\n(AUC = 0.93)`=roc(error.prob.all %>% filter(Ref=='Fin37') %>% 
                                                      transmute(., ERRS=replace(.$ERRS, .$ERRS==2, 1)) %>% .[,1], 
                                                    error.prob.all %>% filter(Ref=='Fin37') %>% .$PROB), 
                        `European\n(AUC = 0.89)`=roc(error.prob.all %>% filter(Ref=='Eur37') %>% 
                                                       transmute(., ERRS=replace(.$ERRS, .$ERRS==2, 1)) %>% .[,1], 
                                                     error.prob.all %>% filter(Ref=='Eur37') %>% .$PROB)), alpha=0.8, size=0.5, 
                   aes=c('linetype', 'color')) +
  geom_abline(slope=1, intercept=1, color='red', alpha=0.8, size=0.5, linetype='dotted') +
  scale_color_manual(values=c('grey55', 'black')) +
  theme(panel.background=element_blank(), axis.line=element_line(colour="black", size=.2)) +
  theme(legend.position='top', legend.title=element_blank(), legend.text=element_text(size=11), legend.key=element_blank())


# collect plots and arrange into second panel
p.error.prob.all.cutoffs.panel <- ggarrange(ggplot()+theme_void(),
                                            p.error.prob.all.retained,
                                            ggplot()+theme_void(),
                                            p.error.prob.all.cutoffs,
                                            ggplot()+theme_void(),
                                            p.error.prob.all.cutoffs.perm,
                                            ggplot()+theme_void(),
                                            nrow=1, ncol=7, widths=c(0.1, 1, 0.1, 1, 0.1, 1, 0.1),
                                            labels=c('', 'd', '', 'e', '', 'f', ''), font.label=list(size=16))

# collect plots and print into pdf
pdf('./results/Fig3_PP_ROC_Err_by-ref.pdf', height=7, width=11)
ggarrange(
  ggarrange(p.err.ref, 
            ggplot()+theme_void(), 
            p.error.prob.perm,
            ggplot()+theme_void(),
            p.err.roc,
            nrow=1, ncol=5, widths=c(1.1, 0.05, 0.8, 0.05, 0.8), 
            labels=c('a', '', 'b', '', 'c'), font.label=list(size=16)),
  ggplot()+theme_void(), 
  p.error.prob.all.cutoffs.panel, 
  nrow=3, ncol=1, heights=c(1, 0.1, 1))
dev.off()



## --------------------------------------------------
## FinnGen HLA imputation results
## Figures 4 & 5
## --------------------------------------------------

# imputed allele frequencies by batch in R2 FinnGen
impfrq.batch <- c('A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPB1') %>% map(function(x) {
  print(x)
  tmp <- fread(paste0('./data/FinnGen_R2/R2_', x, '_batchfreqs.tsv'), data.table=F, fill=T)
  tmp <- tmp[nchar(tmp[, 'allele2'])>3, ] # accept only 4-digit imputed
  tmp$allele2 <- paste0(x, '*', tmp$allele2)
  tmp <- left_join(tmp, ref[, 1:2], by=c('allele2'='Var1'))
  data.frame(tmp, Gene=x)
}) %>% do.call(rbind, .)
impfrq.batch$allele2 <- gsub('68167', '68:167', impfrq.batch$allele2)


# calculate difference of batch frequency of each allele to ref frequency
impfrq.batch <- impfrq.batch %>% group_by(allele2) %>% mutate(diff=(batch.freq-value))
ggplot(impfrq.batch, aes(log10(n), log(diff))) + geom_point()

# calculate n and mean freq deviation in batches
impfrq.batch.means <- data.frame(group_by(impfrq.batch, batch) %>% summarise(sum(n)) %>% data.frame,
                                 MEANDIFF=group_by(impfrq.batch, batch) %>% summarise(mean(diff, na.rm=T)) %>% data.frame %>% .[, 2])


## Plot 

# plot figure 4
cairo_pdf('./results/Fig4_R2_AlleleFreqs.pdf', width=8.1, height=10)
ggplot(impfrq.batch %>% filter(., batch!='AxiomGT1_b16.calls', batch!='AxiomGT1_b17.calls'), 
       aes(reorder(factor(allele2), batch.freq*(-1), FUN=median), batch.freq, color=Gene, fill=Gene)) +
  geom_point(shape=19, alpha=0.4, position=position_jitter(w=0.1, h=0), size=2.5) + 
  geom_text(data=impfrq.batch, aes(label='\u2014', y=value), col='grey20', alpha=0.8) +
  facet_wrap(~ Gene, scales='free', ncol=1, strip.position='right') +
  xlab('Gene*Allele') +
  ylab('Allele frequency') +
  theme(legend.position='none', axis.text.x=element_text(angle=25, hjust=1, size=6.3), strip.background=element_rect(fill="grey90")) +
  theme(panel.background=element_blank(), axis.line=element_line(colour="black", size=.2))
dev.off()


## HLA association data

# reference allele associations for tested autoimmune diseases
rep.alleles <- fread('./data/imputation_alleles.txt', data.table=F, na.strings='') %>% tidyr::fill(., Disease)

# read association data
hla.assoc.dat <- list.files('./data/imputation_summaries', 'assoc$', full.names=T) %>% map(function(x) {
  fread(x, data.table=F)
})
names(hla.assoc.dat) <- gsub('_HLA.assoc', '', list.files('./data/imputation_summaries', 'assoc$'), fixed=T)

hla.assoc.dat.2 <- map2(tapply(rep.alleles$Allele, rep.alleles$Disease, function(x) x), c(1, 6, 3, 4, 5, 2), function(x, y) {
  out <- hla.assoc.dat[[y]]
  out <- out[out$Allele %in% x, ]
  out$Allele <- str_pad(out$Allele, ifelse(nchar(out$Allele)==10, 0, 13))
  out
})
hla.assoc.dat.2 <- hla.assoc.dat.2 %>% melt %>% spread(., variable, value) %>% 
  mutate(., MinusError=Estimate-Std..Error, PlusError=Estimate+Std..Error, OR=exp(Estimate)) %>% 
  mutate(., MinusErrorOR=exp(MinusError), PlusErrorOR=exp(PlusError)) %>% arrange(., L1)
colnames(hla.assoc.dat.2)[2] <- 'Disease'
hla.assoc.dat.2 <- mutate(hla.assoc.dat.2, logP=log10(P) %>% `*`(-1))

# plot figure 5
cairo_pdf('./results/Fig5_R2_HLA_assoc.pdf', width=8, height=4)
ggarrange(
  plotHLAassoc(hla.assoc.dat.2 %>% filter(Disease=='T1D'), '    T1D'),
  plotHLAassoc(hla.assoc.dat.2 %>% filter(Disease=='Coeliac disease'), 'Coeliac disease'),
  plotHLAassoc(hla.assoc.dat.2 %>% filter(Disease=='Psoriatic arthritis'), ' Psoriatic arthritis'),
  plotHLAassoc(hla.assoc.dat.2 %>% filter(Disease=='Rheumatoid arthritis'), 'Rheumatoid arthritis'),
  plotHLAassoc(hla.assoc.dat.2 %>% filter(Disease=='Psoriasis'), '  Psoriasis'),
  plotHLAassoc(hla.assoc.dat.2 %>% filter(Disease=='MS'), '    MS'),
  nrow=3, ncol=2, label.x=0
)
dev.off()

