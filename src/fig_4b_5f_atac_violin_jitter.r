library(ggplot2)
library(reshape2)
library(rtracklayer)
library(effsize)
library(purrr)
library(BSgenome.Mmusculus.UCSC.mm9)
library(ggpubr)
library(RColorBrewer)

wilcoxon.category.enriched <- function(global, subset, col) {
  exclude <- subsetByOverlaps(global, setdiff(global, subset))
  df.global <- mcols(exclude)[, col]
  df.subset <- mcols(subset)[, col]
  
  # Drop non finite values (if log was applied)
  g <- df.global[is.finite(df.global)]
  s <- df.subset[is.finite(df.subset)]
  wilcox.test(s, g)
}

# This must be a path with 5kb bins for the datasets mentioned below. ATAC seq mostly
# in.bins.file <- ''
in.iap.file <- '../data/Navarro_2020_IAPEz_consensus.bed'
df <- read.table(in.bins.file, sep='\t', header=T)

iap.ranges <- import(in.iap.file)
bins.grang <- makeGRangesFromDataFrame(df, keep.extra.columns = T)

iap.bins <- subsetByOverlaps(bins.grang, iap.ranges, minoverlap=2000)
iap.df <- data.frame(iap.bins)

select.cols <- c("Martire2019.ESC.H33WT.ATAC","Martire2019.ESC.H33KO.ATAC","Navarro2020.ATAC.H33KO.H32resc",
                 "Navarro2020.ATAC.H33KO.H33mut","Navarro2020.ATAC.H33KO.H33resc","Navarro2020.ATAC.H33KO",
                 "Navarro2020.ATAC.H33KO.smarcad1KD","Navarro2020.ATAC.H33WT","Navarro2020.ATAC.H33WT.smarcad1KD")

# Fig 3c
select.cols <- c("Martire2019.ESC.H33WT.ATAC", "Martire2019.ESC.H33KO.ATAC")

df.selected <- df[, c('chr', 'start','end', select.cols)]
iap.selected <- iap.df[, c('seqnames','start','end', select.cols)]

colnames(df.selected) <- c('chr', 'start','end', 'H3.3WT','H3.3KO')
colnames(iap.selected) <- c('seqnames','start','end', 'H3.3WT','H3.3KO')

df.melted <- melt(df.selected, id.vars=c('chr','start','end'))
iap.melted <- melt(iap.selected, id.vars=c('seqnames','start','end'))

df.melted$variable <- factor(df.melted$variable, levels <- c('H3.3WT','H3.3KO'))
iap.melted$variable <- factor(iap.melted$variable, levels <- c('H3.3WT','H3.3KO'))


fsize <- 7.5

colors <- c("#19c5c9", "#f8766d")
names(colors) <-  c('H3.3WT','H3.3KO')
colScale <- scale_colour_manual(name = "variable",values = colors)

p <- ggplot(df.melted, aes(x=variable, y=value)) +
  stat_compare_means(data=df.melted, method='wilcox.test', label.x=0.7, label.y=10, size=fsize, color='#666666') +
  geom_violin(fill='#888888') +
  geom_jitter(data=iap.melted, aes(x=variable, y=value, color=variable), alpha=0.6, size=0.8) +

  ylab('ATAC-seq RPGC') +
  xlab('') +
  ggtitle('ATAC-seq 5kb bins') +
  stat_compare_means(data=iap.melted, method='wilcox.test', label.x=0.7, label.y=9.4, size=fsize, color="#f8766d") +
  theme_classic() +
  theme(legend.position='none') +
  annotate("text", label="Cohen's d:", x=0.8, y=9, size=fsize) +
  # annotate("text", label="2.51 (large) [2.44,2.58] CI", x=1.5, y=8, color="#19c5ca") +
  # annotate("text", label="9.3e-06 (negligible) [-0.0037, 0.0037] CI", x=2, y=7.5, color="#888888") +
  annotate("text", label="2.51 [2.44,2.58] CI", x=1, y=7.9, color="#f8766d", size=fsize) +
  annotate("text", label="9.3e-06 [-0.0037, 0.0037] CI", x=1.3, y=8.45, color="#666666", size=fsize) +
  colScale +
  ylim(0,10)

ggsave(plot=p, 'fig3c_stats_annotated.pdf', height=9, width=8, dpi=300)


# Fig 4b
select.cols <- c("Navarro2020.ATAC.H33KO.H32resc",
                 "Navarro2020.ATAC.H33KO.H33mut","Navarro2020.ATAC.H33KO.H33resc","Navarro2020.ATAC.H33KO")

df.selected <- df[, c('chr', 'start','end', select.cols)]
iap.selected <- iap.df[, c('seqnames','start','end', select.cols)]

colnames(df.selected) <- c('chr', 'start','end', '+H3.2','+H3.3mut', '+H3.3', '-')
colnames(iap.selected) <- c('seqnames','start','end', '+H3.2','+H3.3mut', '+H3.3', '-')

df.melted <- melt(df.selected, id.vars=c('chr','start','end'))
iap.melted <- melt(iap.selected, id.vars=c('seqnames','start','end'))

df.melted$variable <- factor(df.melted$variable, levels <- c(  '-', '+H3.3', '+H3.2','+H3.3mut'))
iap.melted$variable <- factor(iap.melted$variable, levels <- c(  '-', '+H3.3', '+H3.2','+H3.3mut'))


fsize <- 7.5

p <- ggplot(df.melted, aes(x=variable, y=value)) +
  stat_compare_means(data=df.melted, method='wilcox.test', size=fsize, label='p.signif', color='#888888', ref.group='-') +
  geom_violin(fill='#aaaaaa') +
  geom_jitter(data=iap.melted, aes(x=variable, y=value, color=variable), alpha=0.6, size=0.8) +
  
  ylab('ATAC-seq RPGC') +
  xlab('') +
  ggtitle('ATAC-seq 5kb bins') +
  stat_compare_means(data=iap.melted, method='wilcox.test', label='p.signif', label.y=9.4, size=fsize, color="#eb786e", ref.group='-') +
  theme_classic() +
  theme(legend.position='none') +
  # annotate("text", label="Cohen's d:", x=0.8, y=8.9, size=fsize) +
  # annotate("text", label="2.51 (large) [2.44,2.58] CI", x=1.5, y=8, color="#19c5ca") +
  # annotate("text", label="9.3e-06 (negligible) [-0.0037, 0.0037] CI", x=2, y=7.5, color="#888888") +
  # annotate("text", label="2.51 [2.44,2.58] CI", x=1, y=7.9, color="#19c5ca", size=fsize) +
  # annotate("text", label="9.3e-06 [-0.0037, 0.0037] CI", x=1.3, y=8.45, color="#888888", size=fsize) +
  
  ylim(0,10)

ggsave(plot=p, 'fig4b_stats_annotated.pdf', height=9, width=11, dpi=300)

# Effect sizes
cohen.d(df.selected[, "+H3.3"], df.selected[, "-"])
cohen.d(df.selected[, "+H3.2"], df.selected[, "-"])
cohen.d(df.selected[, "+H3.3mut"], df.selected[, "-"])

cohen.d(iap.selected[, "+H3.3"], iap.selected[, "-"])
cohen.d(iap.selected[, "+H3.2"], iap.selected[, "-"])
cohen.d(iap.selected[, "+H3.3mut"], iap.selected[, "-"])


  
# Fig 5f
select.cols <- c("Navarro2020.ATAC.H33KO",
                 "Navarro2020.ATAC.H33KO.H33resc", "Navarro2020.ATAC.H33KO.smarcad1KD")

df.selected <- df[, c('chr', 'start','end', select.cols)]
iap.selected <- iap.df[, c('seqnames','start','end', select.cols)]

colnames(df.selected) <- c('chr', 'start','end', '-','+H3.3', 'Smarcad1KD')
colnames(iap.selected) <- c('seqnames','start','end', '-','+H3.3', 'Smarcad1KD')

df.melted <- melt(df.selected, id.vars=c('chr','start','end'))
iap.melted <- melt(iap.selected, id.vars=c('seqnames','start','end'))

df.melted$variable <- factor(df.melted$variable, levels <- c(  '-','+H3.3', 'Smarcad1KD'))
iap.melted$variable <- factor(iap.melted$variable, levels <- c(   '-','+H3.3', 'Smarcad1KD'))


fsize <- 7.5

p <- ggplot(df.melted, aes(x=variable, y=value)) +
  stat_compare_means(data=df.melted, method='wilcox.test', size=fsize, label='p.signif', color='#888888', ref.group='-') +
  geom_violin(fill='#aaaaaa') +
  geom_jitter(data=iap.melted, aes(x=variable, y=value, color=variable), alpha=0.6, size=0.8) +
  ylab('ATAC-seq RPGC') +
  xlab('') +
  ggtitle('ATAC-seq 5kb bins') +
  stat_compare_means(data=iap.melted, method='wilcox.test', label='p.signif', label.y=9.4, size=fsize, color="#eb786e", ref.group='-') +
  theme_classic() +
  theme(legend.position='none') +
  # annotate("text", label="Cohen's d:", x=0.8, y=8.9, size=fsize) +
  # annotate("text", label="2.51 (large) [2.44,2.58] CI", x=1.5, y=8, color="#19c5ca") +
  # annotate("text", label="9.3e-06 (negligible) [-0.0037, 0.0037] CI", x=2, y=7.5, color="#888888") +
  # annotate("text", label="2.51 [2.44,2.58] CI", x=1, y=7.9, color="#19c5ca", size=fsize) +
  # annotate("text", label="9.3e-06 [-0.0037, 0.0037] CI", x=1.3, y=8.45, color="#888888", size=fsize) +
  
  ylim(0,10)

ggsave(plot=p, 'fig5f_stats_annotated.pdf', height=9, width=11, dpi=300)

cohen.d(df.selected[, "+H3.3"], df.selected[, "-"])
cohen.d(df.selected[, "Smarcad1KD"], df.selected[, "-"])

cohen.d(iap.selected[, "+H3.3"], iap.selected[, "-"])
cohen.d(iap.selected[, "Smarcad1KD"], iap.selected[, "-"])
