

## --------------------------------------------------
## Functions
## --------------------------------------------------

# uses donor imputation error data to select the best post prob imputation ref
correctErrs02D <- function(cc, wvec=c(1,1,1), inds) {
  tmp <- cbind(reg38.hlageno.d.err[inds, cc], reg37.hlageno.d.err[inds, cc], def.hlageno.d.err[inds, cc])
  tmp[, 2] <- tmp[, 2]*wvec[1]
  tmp[, 4] <- tmp[, 4]*wvec[2]
  tmp[, 6] <- tmp[, 6]*wvec[3]
  tmp.err <- sapply(1:nrow(tmp), function(x) {
    c(tmp[x, 1], tmp[x, 3], tmp[x, 5])[which.max(c(tmp[x, 2], tmp[x, 4], tmp[x, 6]))]  }) %>% t
  c(sum(tmp.err), sum(tmp[, 1]) , sum(tmp[, 3]), sum(tmp[, 5]))
}

# uses recipient imputation error data to select the best post prob imputation ref
correctErrs02R <- function(cc, wvec=c(1,1,1), inds) {
  tmp <- cbind(reg38.hlageno.r.err[inds, cc], reg37.hlageno.r.err[inds, cc], def.hlageno.r.err[inds, cc])
  tmp[, 2] <- tmp[, 2]*wvec[1]
  tmp[, 4] <- tmp[, 4]*wvec[2]
  tmp[, 6] <- tmp[, 6]*wvec[3]
  tmp.err <- sapply(1:nrow(tmp), function(x) {
    c(tmp[x, 1], tmp[x, 3], tmp[x, 5])[which.max(c(tmp[x, 2], tmp[x, 4], tmp[x, 6]))]  }) %>% t
  c(sum(tmp.err), sum(tmp[, 1]) , sum(tmp[, 3]), sum(tmp[, 5]))
}


# plot HLA association data
plotHLAassoc <- function(d, main.text) {
  
  d <- arrange(d, OR)
  d$Allele <- factor(d$Allele, levels=d$Allele)
  
  p1.1 <- ggplot(d, aes(x=Allele, y=OR)) +
    geom_pointrange(aes(ymin=MinusErrorOR, ymax=PlusErrorOR), alpha=0.7, fatten=3) +
    ylab('OR') + xlab('') +
    ggtitle('') +
    geom_hline(yintercept=0, linetype='dotted', color='grey30', size=0.25) +
    theme(axis.text.y=element_text(angle=0, hjust=1, size=7.3, color='red')) + 
    theme(panel.background=element_blank(), axis.line=element_line(colour="black", size=.2)) +
    coord_flip() 
  
  p1.2 <- ggplot(d, aes(x=Allele, y=logP)) +
    geom_bar(stat='identity', alpha=0.5, width=0.95) +
    ylab( expression('-log'[10]*'('*italic(p)*')') ) +
    xlab('') +
    ggtitle('') +
    coord_flip() +
    ylim(0, 130) +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
    theme(panel.background=element_blank(), axis.line=element_line(colour="black", size=.2))
  
  ggarrange(p1.1, p1.2, nrow=1, ncol=2, widths=c(1, 0.7),
            labels=main.text, font.label=list(face='italic', size=12), hjust=-0.2)

}



