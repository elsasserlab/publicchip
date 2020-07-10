library(ggplot2)

tabledir <- './tables'
figdir <- './figures'

# Note: These plots are generated from 5kb bins files, and intersections of
# bigWig files with bed files provided as supplementary data: IAPEz_consensus,
# ChromHMM17, RepMasker_lt200bp_all and TSS_hi.

# These tables can be downloaded from the extra data folder.

# Figure 6F - 5kb bins Scatterplot H3.3 wildtype vs Smarcad1KD vs H3K9me3,
# IAP-overlapping bins highlighted

bins.df <- read.csv(paste(tabledir,
                          'fig06_f_H33WT_vs_Smarcad1KD_vs_H3K9me3_5kb_bins.csv',
                          sep='/'))

iap.bins.df <- read.csv(paste(tabledir,
                              'fig06_f_H33WT_vs_Smarcad1KD_vs_H3K9me3_5kb_bins_at_IAP.csv',
                              sep='/'))

# Drop zeros in the input
bins.df <- bins.df[bins.df$Navarro2020_WT_Input > 0 & bins.df$Navarro2020_Smarcad1KD_Input > 0, ]

bins.df$H33wt_norm <- bins.df$Navarro2020_WT_H33 / bins.df$Navarro2020_WT_Input
bins.df$H33smarcad1KD_norm <- bins.df$Navarro2020_Smarcad1KD_H33 / bins.df$Navarro2020_Smarcad1KD_Input

# Drop zeros in smarcad1KD norm, if any
bins.df <- bins.df[bins.df$H33smarcad1KD_norm > 0, ]

bins.df$H33fc_norm <- log2(bins.df$H33wt_norm / bins.df$H33smarcad1KD_norm)

iap.bins.df$H33wt_norm <- iap.bins.df$Navarro2020_WT_H33 / iap.bins.df$Navarro2020_WT_Input
iap.bins.df$H33smarcad1KD_norm <- iap.bins.df$Navarro2020_Smarcad1KD_H33 / iap.bins.df$Navarro2020_Smarcad1KD_Input
iap.bins.df$H33fc_norm <- log2(iap.bins.df$H33wt_norm / iap.bins.df$H33smarcad1KD_norm)

ggplot(bins.df, aes(x=H33fc_norm, y=Shi2019_H3K9me3)) + 
  geom_bin2d(binwidth=c(0.02,0.02)) +
  geom_point(data=iap.bins.df, aes(x=H33fc_norm, y=Shi2019_H3K9me3), color='#32a472', size=0.8, alpha=0.6) +
  scale_fill_gradient(low='#dddddd', high='#333333') + #, high='#bb3215') +
  theme_classic() +
  xlim(-3, 3) +
  ylim(0, 30) +
  xlab("log2(H33 WT/ Smarcad1KD)") +
  ylab("Shi2019 H3K9me3 RPGC") +
  geom_vline(xintercept=0, linetype='dashed') +
  ggtitle("Smarcad1 vs H33 WT / Smarcad1KD (5kb bins)")

ggsave(paste(figdir, 'fig06_f_H33wt_vs_Smarcad1KD_vs_H3K9me3_scatterplot.png', sep='/'))



