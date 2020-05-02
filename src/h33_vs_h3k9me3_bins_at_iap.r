library(ggplot2)
library(reshape)
library(effsize)
library(rtracklayer)

compute.effect.size.vs.global <- function(global, subset, col) {
  exclude <- subsetByOverlaps(global, setdiff(global, subset))
  df.global <- mcols(exclude)[, col]
  df.subset <- mcols(subset)[, col]
  
  # Drop non finite values (if log was applied)
  g <- df.global[is.finite(df.global)]
  s <- df.subset[is.finite(df.subset)]
  
  cohen.d(s, g)
}

wilcoxon.category.enriched <- function(global, subset, col) {
  exclude <- subsetByOverlaps(global, setdiff(global, subset))
  df.global <- mcols(exclude)[, col]
  df.subset <- mcols(subset)[, col]
  
  # Drop non finite values (if log was applied)
  g <- df.global[is.finite(df.global)]
  s <- df.subset[is.finite(df.subset)]
  wilcox.test(s, g)
}

args = commandArgs(trailingOnly=TRUE)

in.bins.file <- args[1]
in.iap.file <- args[2]
outfile.prefix <- args[3]

df <- read.table(in.bins.file, sep='\t', header=T)
values <- df[, c('chr', 'start', 'end', 'H33.00h','H33.inp','H3K9me3.WT','H3K9me3.input')]
values$H33logfc <- log2(values$H33.00h/values$H33.inp)
values$H3k9me3.logfc <- log2(values$H3K9me3.WT/values$H3K9me3.input)

iap.ranges <- import(in.iap.file)
bins.grang <- makeGRangesFromDataFrame(values, keep.extra.columns = T)
iap.bins <- subsetByOverlaps(bins.grang, iap.ranges, minoverlap=2500)

iap.bins.df <- data.frame(H33=iap.bins$H33logfc, H3K9me3=iap.bins$H3k9me3.logfc)

# Calculate effect sizes
k9.estimate <- compute.effect.size.vs.global(bins.grang, iap.bins, 'H3k9me3.logfc')
h33.estimate <- compute.effect.size.vs.global(bins.grang, iap.bins, 'H33logfc')

# Without overlay
p1 <- ggplot(values, aes(x=H33logfc,y=H3k9me3.logfc)) +
  geom_bin2d(binwidth=c(0.1,0.1)) +
  ylim(-10, 7) +
  scale_fill_gradient(low="#eeeeee", high="#bb0000") +
  theme_classic()

# With overlay 
p2 <- ggplot(values, aes(x=H33logfc,y=H3k9me3.logfc)) +
  geom_bin2d(binwidth=c(0.1,0.1)) +
  ylim(-10, 7) +
  geom_point(data=iap.bins.df, aes(x=H33,y=H3K9me3), color='black', alpha=0.8, size=0.8) +
  scale_fill_gradient(low="#eeeeee", high="#bb0000") +
  annotate(geom="text", x=-1.2, y=-8.5, label="Wilcoxon rank sum test p-val < 2.2e-16 (both axes)",
           color="black") +
  annotate(geom="text", x=-1.2, y=-9.5,
           label=paste0("Cohen's d effect: ", round(h33.estimate$estimate, digits=2) ,", ", round(k9.estimate$estimate,digits=2)," (H3.3,H3K9me3)"),
           color="black") +
  theme_classic()

ggsave(paste0(outfile.prefix, '_bins_density.png'), plot=p1, dpi=300)
ggsave(paste0(outfile.prefix, '_bins_density_at_iap.png'), plot=p2, dpi=300)

