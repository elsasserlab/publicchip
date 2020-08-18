library(ggplot2)
library(reshape2)
library(ggpubr)
library(effsize)
# Note: These plots are generated from 5kb bins files, and intersections of
# bigWig files with bed files provided as supplementary data: IAPEz_consensus,
# ChromHMM17, RepMasker_lt200bp_all and TSS_hi.

# These tables can be downloaded from the extra data folder.

tabledir <- './tables'
figdir <- './figures'

# Figure 6C - Boxplot H3.3 ChIP-seq, H3K9me3 ChIP-seq Smarcad1KD vs WT
# at IAP ERVs and TSS highly expressed.

bin.id <- c('seqnames', 'start', 'end')
df <- read.csv(paste(tabledir,
                     'fig06c_and_ext12b_atac_and_chip_at_iap.csv', sep='/'))

columns <- c('Navarro2020_WT_Input',
             'Navarro2020_WT_H33',
             'Navarro2020_WT_H3K9me3',
             'Navarro2020_Smarcad1KD_Input',
             'Navarro2020_Smarcad1KD_H33',
             'Navarro2020_Smarcad1KD_H3K9me3')

df <- df[, c(bin.id, columns)]
df <- df[df$Navarro2020_WT_Input > 0 & df$Navarro2020_Smarcad1KD_Input > 0, ]

df$WT_H33_norm <- df$Navarro2020_WT_H33/df$Navarro2020_WT_Input
df$Smarcad1KD_H33_norm <- df$Navarro2020_Smarcad1KD_H33/df$Navarro2020_Smarcad1KD_Input
df$WT_H3K9me3_norm <- df$Navarro2020_WT_H3K9me3/df$Navarro2020_WT_Input
df$Smarcad1KD_H3K9me3_norm <- df$Navarro2020_Smarcad1KD_H3K9me3/df$Navarro2020_Smarcad1KD_Input


fig.df <- melt(df[, c(bin.id, 'WT_H33_norm', 'Smarcad1KD_H33_norm', 'WT_H3K9me3_norm', 'Smarcad1KD_H3K9me3_norm')], id.vars=bin.id)

ggplot(fig.df, aes(x=variable, y=value, color=variable)) +
  geom_boxplot(color="black",outlier.alpha = 0,notch = T) +
  geom_jitter(aes(x=variable, y=value, color=variable), alpha=0.15, size=0.5) +
  ylim(0,25) + 
  xlab('') +
  ylab('ChIP RPGC / Input') +
  ggtitle('IAP ERVs') +
  theme(legend.position='none', axis.text.x=element_text(angle=45, hjust=1))

ggsave(paste(figdir, 'fig06c_1_chip_at_IAP.png', sep='/'))

compare_means(value ~ variable, data=fig.df,
              method='t.test', paired=T)

cohen.d(df$Smarcad1KD_H33_norm, df$WT_H33_norm)
cohen.d(df$Smarcad1KD_H3K9me3_norm, df$WT_H3K9me3_norm)


## Same at TSS_hi
bin.id <- c('seqnames', 'start', 'end')
df <- read.csv(paste(tabledir,
                     'fig06c_chip_at_tss.csv', sep='/'))

columns <- c('Navarro2020_WT_Input',
             'Navarro2020_WT_H33',
             'Navarro2020_WT_H3K9me3',
             'Navarro2020_Smarcad1KD_Input',
             'Navarro2020_Smarcad1KD_H33',
             'Navarro2020_Smarcad1KD_H3K9me3')

df <- df[, c(bin.id, columns)]
df <- df[df$Navarro2020_WT_Input > 0 & df$Navarro2020_Smarcad1KD_Input > 0, ]

df$WT_H33_norm <- df$Navarro2020_WT_H33/df$Navarro2020_WT_Input
df$Smarcad1KD_H33_norm <- df$Navarro2020_Smarcad1KD_H33/df$Navarro2020_Smarcad1KD_Input
df$WT_H3K9me3_norm <- df$Navarro2020_WT_H3K9me3/df$Navarro2020_WT_Input
df$Smarcad1KD_H3K9me3_norm <- df$Navarro2020_Smarcad1KD_H3K9me3/df$Navarro2020_Smarcad1KD_Input

fig.df <- melt(df[, c(bin.id, 'WT_H33_norm', 'Smarcad1KD_H33_norm', 'WT_H3K9me3_norm', 'Smarcad1KD_H3K9me3_norm')], id.vars=bin.id)

ggplot(fig.df, aes(x=variable, y=value, color=variable)) +
  geom_boxplot(color="black",outlier.alpha = 0,notch = T) +
  geom_jitter(aes(x=variable, y=value, color=variable), alpha=0.15, size=0.5) +
  ylim(0,10) + 
  xlab('') +
  ylab('ChIP RPGC / Input') +
  ggtitle('TSS highly expressed') +
  theme(legend.position='none', axis.text.x=element_text(angle=45, hjust=1))

ggsave(paste(figdir, 'fig06c_2_chip_at_TSS.png', sep='/'))

compare_means(value ~ variable, data=fig.df,
              method='t.test', paired=T)

cohen.d(df$Smarcad1KD_H33_norm, df$WT_H33_norm)
cohen.d(df$Smarcad1KD_H3K9me3_norm, df$WT_H3K9me3_norm)

