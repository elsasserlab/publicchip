library(ggplot2)
library(reshape)
library(rtracklayer)
library(effsize)
library(ggpubr)


remove.outliers <- function(df, column) {
  # Remove outliers by Tukey's approach (the one implemented in boxplot)
  outliers <- boxplot(df[,column], plot=F)$out
  result <- df[! df[,column] %in% outliers, ]
  result
}

remove.outliers.melted <- function(dfmelt) {
  df <- split(dfmelt, dfmelt[,c('variable')])
  noout <- lapply(df, remove.outliers, column='value')
  do.call(rbind, noout)
}

args = commandArgs(trailingOnly=TRUE)

in.bins.file <- args[1]
in.iap.file <- args[2]
outfile <- args[3]

df <- read.table(in.bins.file, sep='\t', header=T)
values <- df[, c('chr', 'start', 'end', "Navarro2019.ATAC.H33KO", "Navarro2019.ATAC.H33KO.H33resc", "Navarro2019.ATAC.H33KO.smarcad1KD")]

colnames(values) <- c('chr', 'start', 'end', "-", "+H33", "smarcad1KD")

iap.ranges <- import(in.iap.file)
bins.grang <- makeGRangesFromDataFrame(values, keep.extra.columns = T)
iap.bins <- subsetByOverlaps(bins.grang, iap.ranges, minoverlap=2500)
iap.df <- data.frame(mcols(iap.bins))
colnames(iap.df) <- c("-", "+H33", "smarcad1KD")

df.melted <- melt(values, id.vars=c('chr','start','end'))
iap.melted <- melt(iap.df)

df.noout <- remove.outliers.melted(df.melted)
iap.melted.noout <- remove.outliers.melted(iap.melted)

p <- ggplot(df.noout, aes(x=variable, y=value)) +
  geom_violin(fill='#bbbbbb') +
  geom_jitter(data=iap.melted.noout, aes(x=variable, y=value, color=variable), alpha=0.7, size=0.8) +
  ylab('RPGC') +
  xlab('H33KO') +
  ggtitle('ATAC-seq 5kb bins H3.3 KO') +
  theme_classic() +
  stat_compare_means() +
  theme(legend.position='none')

ggsave(outfile, plot=p, dpi=300)

