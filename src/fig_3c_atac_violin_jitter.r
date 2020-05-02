library(ggplot2)
library(reshape)
library(rtracklayer)
library(effsize)

args = commandArgs(trailingOnly=TRUE)

in.bins.file <- args[1]
in.iap.file <- args[2]
outfile <- args[3]

df <- read.table(in.bins.file, sep='\t', header=T)
values <- df[, c('chr', 'start', 'end', 'Martire2019.ATAC.H33WT','Martire2019.ATAC.H33KO')]

iap.ranges <- import(in.iap.file)
bins.grang <- makeGRangesFromDataFrame(values, keep.extra.columns = T)
iap.bins <- subsetByOverlaps(bins.grang, iap.ranges, minoverlap=2500)
iap.df <- data.frame(mcols(iap.bins))

df.melted <- melt(values, id.vars=c('chr','start','end'))
iap.melted <- melt(iap.df)

p <- ggplot(df.melted, aes(x=variable, y=value)) +
  geom_violin(fill='#444444') +
  geom_jitter(data=iap.melted, aes(x=variable, y=value, color=variable), alpha=0.7) +
  ylab('ATAC-seq RPGC') +
  xlab('') +
  ggtitle('ATAC-seq 5kb bins H3.3 KO vd WT at IAP') +
  theme_classic() +
  theme(legend.position='none') +
  ylim(0,10)

ggsave(outfile, plot=p, dpi=300)

# Wilcoxon ranked sum test significance for global ATAC H33KO vs WT
wilcox.test(values$Martire2019.ATAC.H33KO, values$Martire2019.ATAC.H33WT)

# Wilcoxon ranked sum test significance for iap ATAC H33KO vs WT
wilcox.test(iap.df$Martire2019.ATAC.H33KO, iap.df$Martire2019.ATAC.H33WT)

# Effect sizes
cohen.d(iap.df$Martire2019.ATAC.H33KO, iap.df$Martire2019.ATAC.H33WT)
cohen.d(values$Martire2019.ATAC.H33KO, values$Martire2019.ATAC.H33WT)

