library(ggplot2)
library(reshape)
library(rtracklayer)
library(effsize)

wilcoxon.subset <- function(global, subset, col) {
  exclude <- subsetByOverlaps(global, setdiff(global, subset))
  df.global <- mcols(exclude)[, col]
  df.subset <- mcols(subset)[, col]
  
  # Drop non finite values (if log was applied)
  g <- df.global[is.finite(df.global)]
  s <- df.subset[is.finite(df.subset)]
  wilcox.test(s, g)
}

compute.effect.size.vs.global <- function(global, subset, col) {
  exclude <- subsetByOverlaps(global, setdiff(global, subset))
  df.global <- mcols(exclude)[, col]
  df.subset <- mcols(subset)[, col]
  
  # Drop non finite values (if log was applied)
  g <- df.global[is.finite(df.global)]
  s <- df.subset[is.finite(df.subset)]
  
  cohen.d(s, g)
}

args = commandArgs(trailingOnly=TRUE)

in.bins.file <- args[1]
in.iap.file <- args[2]
outfile <- args[3]

df <- read.table(in.bins.file, sep='\t', header=T)
values <- df[, c('chr', 'start', 'end', 'Martire2019.ATAC.H33WT','Martire2019.ATAC.H33KO', 'H3K9me3.WT')]
values$ATAC.logfc <- log2(values$Martire2019.ATAC.H33KO/values$Martire2019.ATAC.H33WT)

iap.ranges <- import(in.iap.file)
bins.grang <- makeGRangesFromDataFrame(values, keep.extra.columns = T)
iap.bins <- subsetByOverlaps(bins.grang, iap.ranges, minoverlap=2500)
iap.df <- data.frame(mcols(iap.bins))

df.melted <- melt(values, id.vars=c('chr','start','end'))
iap.melted <- melt(iap.df)

p <- ggplot(values, aes(x=ATAC.logfc, y=H3K9me3.WT)) +
  geom_point(color='#aaaaaa', alpha=0.2, size=0.8) +
  geom_point(data=iap.df, aes(x=ATAC.logfc, y=H3K9me3.WT), color='#409c70', alpha=0.7, size=0.8) +
  ylab('H3K9me3 RPGC') +
  xlab('ATAC-seq log2 H3.3KO/WT') +
  ggtitle('5kb bins') +
  theme_classic() +
  theme(legend.position='none') +
  geom_vline(xintercept=0, linetype='dotted') +
  ylim(0,30) +
  xlim(-5,5)

ggsave(outfile, plot=p, dpi=300)

# Wilcoxon ranked sum test significance for global values vs values at IAP bins
wilcoxon.subset(bins.grang, iap.bins, 'ATAC.logfc')
wilcoxon.subset(bins.grang, iap.bins, 'H3K9me3.WT')

# Effect size calculation
compute.effect.size.vs.global(bins.grang, iap.bins, 'ATAC.logfc')
compute.effect.size.vs.global(bins.grang, iap.bins, 'H3K9me3.WT')

