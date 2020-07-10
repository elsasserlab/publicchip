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

# Figure 12 extended b - Boxplot ATAC-seq Smarcad1KD vs WT at IAP ERVs 

bin.id <- c('seqnames', 'start', 'end')
df <- read.csv(paste(tabledir,
                     'fig06c_and_ext12b_atac_and_chip_at_iap.csv', sep='/'))

columns <- colnames(df)[4:10]

df <- df[, c(bin.id, columns)]

fig.df <- melt(df, id.vars=bin.id)

ggplot(fig.df, aes(x=variable, y=value, color=variable)) +
  geom_boxplot(color="black",outlier.alpha = 0,notch = T) +
  geom_jitter(aes(x=variable, y=value, color=variable), alpha=0.15, size=0.5) +
  theme_elsasserlab_screen(base_size = 16) +
  ylim(0,3) + 
  xlab('') +
  ylab('ChIP RPGC / Input') +
  ggtitle('ATAC-seq at IAP ERVs') +
  theme(legend.position='none', axis.text.x=element_text(angle=45, hjust=1))

ggsave(paste(figdir, 'fig12ext_b_atac_seq_at_IAP.png', sep='/'))

compare_means(value ~ variable, data=fig.df,
              method='t.test', paired=T)

cohen.d(df$H33WT_Smarcad1KD, df$H33WT)
cohen.d(df$H33WT_ATRXKD, df$H33WT)
cohen.d(df$H33KO_Smarcad1KD, df$H33KO)
cohen.d(df$H33KO_H33rescue_Smarcad1KD, df$H33KO_H33rescue)

# Shuffled

bin.id <- c('seqnames', 'start', 'end')
df <- read.csv(paste(tabledir,
                     'fig06c_and_ext12b_atac_and_chip_at_iap_shuffled.csv', sep='/'))

columns <- colnames(df)[4:10]

df <- df[, c(bin.id, columns)]

fig.df <- melt(df, id.vars=bin.id)

ggplot(fig.df, aes(x=variable, y=value, color=variable)) +
  geom_boxplot(color="black",outlier.alpha = 0,notch = T) +
  geom_jitter(aes(x=variable, y=value, color=variable), alpha=0.15, size=0.5) +
  theme_elsasserlab_screen(base_size = 16) +
  ylim(0,3) + 
  xlab('') +
  ylab('ChIP RPGC / Input') +
  ggtitle('ATAC-seq at IAP ERVs') +
  theme(legend.position='none', axis.text.x=element_text(angle=45, hjust=1))

ggsave(paste(figdir, 'fig12ext_b_atac_seq_at_IAP_shuffled.png', sep='/'))

compare_means(value ~ variable, data=fig.df,
              method='t.test', paired=T)

cohen.d(df$H33WT_Smarcad1KD, df$H33WT)
cohen.d(df$H33WT_ATRXKD, df$H33WT)
cohen.d(df$H33KO_Smarcad1KD, df$H33KO)
cohen.d(df$H33KO_H33rescue_Smarcad1KD, df$H33KO_H33rescue)



