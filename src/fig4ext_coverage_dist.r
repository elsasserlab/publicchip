library(cowplot)
library(ggplot2)
library(reshape2)
library(rtracklayer)
library(ggpubr)
library(GenomicRanges)

intersect.bw.with.bed <- function(bw.file, bed.file, stat='mean') {
  bed <- import(bed.file)
  bw <- BigWigFile(bw.file)
  # summarize over bed regions using stat
  unlist(summary(bw, bed, type=stat))
}

iapez.file <- '../data/Navarro_2020_IAPEz_consensus.1000bp.flanked.bed'
etn.file <- '../data/Navarro_2020_RepMasker_clusters_2kb.ETn.1000bp.flanked.bed'

files <- c(
  '../intermediate/bw/Shi_2019/ESC1_WT_H3K9me3.ext.uniq.bw',
  '../intermediate/bw/Bunch_2014/Bunch_S2_PolII_WT_rep1.ext.uniq.bw',
  '../intermediate/bw/Bunch_2014/Bunch_Total_PolII_WT_rep1.ext.uniq.bw',
  '../intermediate/bw/Deaton_2016/ES_H33_00h_rep1.ext.uniq.bw'
)

values <- lapply(files, intersect.bw.with.bed, bed.file=etn.file)

etn.granges <- values[[1]]
etn.granges$h3k9m3 <- values[[1]]$score
etn.granges$rpol2.s2 <- values[[2]]$score
etn.granges$rpol2.total <- values[[3]]$score
etn.granges$H33.0h <- values[[4]]$score
etn.granges$score <- NULL

values <- lapply(files, intersect.bw.with.bed, bed.file=iapez.file)

iapez.granges <- values[[1]]
iapez.granges$h3k9m3 <- values[[1]]$score
iapez.granges$rpol2.s2 <- values[[2]]$score
iapez.granges$rpol2.total <- values[[3]]$score
iapez.granges$H33.0h <- values[[4]]$score
iapez.granges$score <- NULL

df.etn <- as.data.frame(mcols(etn.granges))
df.iap <- as.data.frame(mcols(iapez.granges))

p1 <- ggplot(df.etn, aes(x=h3k9m3, y=rpol2.s2, color=H33.0h)) + geom_point(alpha=0.8) + scale_color_gradient(low="#dddddd", high="#c21515") + xlab("H3K9me3") + ylab("RNA PolII S2") + ggtitle("ETN") + theme_classic() + theme(legend.position="bottom")
p2 <- ggplot(df.iap, aes(x=h3k9m3, y=rpol2.s2, color=H33.0h)) + geom_point(alpha=0.8) + scale_color_gradient(low="#dddddd", high="#c21515") + xlab("H3K9me3") + ggtitle("IAP ERV") + ylab('') + theme_classic() + theme(legend.position="bottom")

plot_grid(p1, p2)
ggsave('../figures/ext4_bc_scatter_both.pdf', width=20, height=10, units='in')
