library(ggplot2)
library(eulerr)

args = commandArgs(trailingOnly=TRUE)
out <- args[1]

# Euler plotting is stochastic. This seed
# reaches a perfect fit.
set.seed(55)

# Values taken from intervene venn output
e <- euler(c(H33.ES=38509,
      H33.NS=61907,
      H3K9me3=52573,
      'H33.ES&H33.NS'=22010,
      'H33.ES&H3K9me3'=18735,
      'H33.NS&H3K9me3'=443,
      'H33.ES&H33.NS&H3K9me3'=278), shape='ellipse')

pdf(out, width=5, height=5)
plot(e, fills = list(fill = c("#f36972", "#f9b772", "#992e99"), alpha = 0.5), quantities=TRUE)
dev.off()
