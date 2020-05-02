#!/usr/bin/env Rscript

library(ggplot2)
library(reshape)

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

infile <- args[1]
outfile.prefix <- args[2]

df <- read.table(infile, sep='\t', header=T)
drops <- c('chr', 'start', 'end')

df <- df[, !(names(df) %in% drops)]

melted.df <- melt(df, id.vars='id')
no.out.df <- remove.outliers.melted(melted.df)

# Plot with outliers
p1 <- ggplot(melted.df, aes(x=variable, y=value, fill=variable)) +
  geom_violin() +
  geom_boxplot(width=0.1, color='#333333', fill='white') +
  xlab('') +
  ylab('') +
  ggtitle('Full bin distribution') +
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(paste(outfile.prefix, 'outliers.pdf', sep='_'), plot=p1, width = 30, height = 22, units = "cm")

# Plot without outliers
p2 <- ggplot(no.out.df, aes(x=variable, y=value, fill=variable)) +
  geom_violin() +
  geom_boxplot(width=0.1, color='#333333', fill='white', outlier.shape=NA) +
  xlab('') +
  ylab('') +
  ggtitle('Bin distribution (no outliers)') +
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(paste(outfile.prefix, 'no_out.pdf', sep='_'), plot=p2, dpi=300, width = 30, height = 22, units = "cm")