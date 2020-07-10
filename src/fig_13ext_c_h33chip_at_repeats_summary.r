library(ggplot2)
library(pheatmap)


# Note: These plots are generated from 5kb bins files, and intersections of
# bigWig files with bed files provided as supplementary data: IAPEz_consensus,
# ChromHMM17, RepMasker_lt200bp_all and TSS_hi.

# These tables can be downloaded from the extra data folder.

# Extended figure 13C - Repeats summary heatmap for H3K9me3, H3.3 and Smarcad1 ChIP.

tabledir <- './tables'
figdir <- './figures'

# Helper plot functions
compute.lims <- function(mat) {
  maxmat <- mat
  maxmat[is.na(maxmat)] <- -Inf
  maxval <- max(maxmat)
  
  minmat <- mat
  minmat[is.na(minmat)] <- Inf
  minval <- min(minmat)
  
  nsteps <- 21.0
  
  breaklim <- ceiling(abs(max(abs(c(minval, maxval)))))
  
  # Compute a reasonable amount of steps?
  stepsize <- 2*breaklim / nsteps
  
  breakslist <- seq(-breaklim, +breaklim, by=stepsize)
  breakslist
}

summary_heatmap <- function(values, title, size=35) {
  bcolor <- "white"
  
  breakslist <- compute.lims(values)
  palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breakslist))
  
  cellsize <- size
  
  pheatmap(values,
           main=title,
           cellwidth=cellsize,
           cluster_rows=F,
           cluster_cols=F,
           cellheight=cellsize,
           border_color=bcolor,
           breaks=breakslist,
           color=palette,
           display_numbers=TRUE)
}

df <- read.csv(paste(tabledir, 'fig12ext_c_repeats_summary.csv', sep='/'), row.names=1)

results_wt <- df[, c("Navarro2020_WT_Input", "Navarro2020_WT_H33", "Navarro2020_WT_H3K9me3")]
results_kd <- df[, c("Navarro2020_Smarcad1KD_Input", "Navarro2020_Smarcad1KD_H33", "Navarro2020_Smarcad1KD_H3K9me3" )]
results_smarcad <- df[, c("Smarcad1_Input_FLAG", "Smarcad1_FLAG")]

results_norm <- cbind(log2(sweep(results_wt, 1, results_wt$Navarro2020_WT_Input, `/`)),
                      log2(sweep(results_kd, 1, results_kd$Navarro2020_Smarcad1KD_Input, `/`)),
                      log2(sweep(results_smarcad, 1, results_smarcad$Smarcad1_Input_FLAG, `/`)))

p <- summary_heatmap(t(results_norm[, c("Navarro2020_WT_H33",
                                   "Navarro2020_Smarcad1KD_H33",
                                   "Navarro2020_WT_H3K9me3",
                                   "Navarro2020_Smarcad1KD_H3K9me3")]),
                title="RepMasker annotation (>200bp elements) (log2fc over input)")

ggsave(plot=p, paste(figdir, 'fig_13ext_repeats_summary_heatmap_a.png', sep='/'), width=20, height=5)

results_delta <- cbind(results_norm$Navarro2020_Smarcad1KD_H33 - results_norm$Navarro2020_WT_H33,
                       results_norm$Navarro2020_Smarcad1KD_H3K9me3 - results_norm$Navarro2020_WT_H3K9me3)

colnames(results_delta) <- c('delta_H33', 'delta_H3K9me3')

p <- summary_heatmap(t(results_delta),
                     title="RepMasker annotation (>200bp elements) (log2fc over input)")

ggsave(plot=p, paste(figdir, 'fig_13ext_repeats_summary_heatmap_b.png', sep='/'), width=20, height=5)


p <- summary_heatmap(t(results_norm[, c("Smarcad1_FLAG")]),
                     title="RepMasker annotation (>200bp elements) (log2fc over input)")
ggsave(plot=p, paste(figdir, 'fig_13ext_repeats_summary_heatmap_c.png', sep='/'), width=20, height=3)


