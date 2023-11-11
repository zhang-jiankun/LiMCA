#!/usr/bin/env Rscript

# See https://stephanecastel.wordpress.com/2017/02/15/how-to-generate-ase-data-with-phaser/

args <- commandArgs()

ase_table <- args[6]
output <- args[7]

# load haplotypic counts
ase_counts = read.delim(ase_table)

# select only genes with sufficient coverage
cov = subset(ase_counts, ase_counts$totalCount>=10)

# perform binomial test for deviation from 0.5
cov$binom_p = apply(cov[, c("aCount", "bCount")],1,function(x){binom.test(x[1], x[1]+x[2], p=0.5)$p.value})

# perform multiple testing correction with FDR
cov$binom_q = p.adjust(cov$binom_p, method="fdr")

pdf(output, width=5, height=5)
# Plot haplotype A versus haplotype, highlight sites with significant imbalance (FDR<5%)
plot(
  cov$bCount, 
  cov$aCount, 
  log="xy", 
  col=(cov$binom_q<0.05)+1, 
  ylab="Haplotype B count",
  xlab="Haplotype A count",
  pch=19,
  cex=0.5
)
abline(0, 1, col="grey")
legend("topleft", c("No significant imbalance", "Significant imbalance"), pch=c(15,15), col=c(1,2))

# QQ plot of pvalues

ggd.qqplot = function(pvector, main=NULL, ...){
  o = -log10(sort(pvector, decreasing=F))
  e = -log10(1:length(o)/length(o))
  plot(e, o, pch=19, cex=1, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0, max(e)), ylim=c(0,max(o)))
       lines(e,e, col="red")
}
ggd.qqplot(cov[cov$binom_p>0,]$binom_p)
dev.off()