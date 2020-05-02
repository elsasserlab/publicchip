library(ggplot2)
library(reshape2)
library(rtracklayer)
library(ggpubr)
library(purrr)
library(GenomicRanges)
library(pheatmap)

source('./code/bwutils.r')
source('./code/heatmap_utils.r')

bw.dir <- './data/bw/'
out.dir <- './final_figures/07_ext/'

repmasker.file <- './data/bed/RepMasker_lt200bp_all.bed'
chromhmm.file <- './data/bed/ChromHMM17.bed'

# Keep the order
levels.order <- c('1_Txn_Elongation',
                  '2_Weak_Txn',
                  '3_Txn_Transition',
                  '4_Poised_Enhancer',
                  '5_Active_Promoter',
                  '6_Strong_Enhancer',
                  '7_Active_Promoter',
                  '8_Strong_Enhancer',
                  '9_Strong_Enhancer',
                  '10_Poised_Promoter',
                  '11_Repressed',
                  '12_Heterochrom',
                  '13_Heterochrom',
                  '14_Heterochrom',
                  '15_Insulator',
                  '16_H3K9_Hetero',
                  '17_H3K9_H33_Hetero')

# Panel 1: ESC H3.3 SNAP time log2 fold change normalized to input replicate 1
files <- paste0(bw.dir, c(
  'Martire_2019/ESC_H33WT_ATAC.bw',
  'Martire_2019/ESC_H33KO_ATAC.bw'
  
))

summaries <- map(files, summarize.by.bed, bed.file <- chromhmm.file, drop.name.col=F)
df.fg <- Reduce(function(...) merge(..., all=TRUE), summaries)
rownames(df.fg) <- df.fg$name
df.fg$name <- NULL

out.file <- paste0(out.dir, 'ATAC_signal_at_chromhmm.pdf')
colnames(df.fg) <- c('H3.3 WT', 'H3.3 KO')
summarized.heatmap(t(df.fg[levels.order, ]), out.file, "ATAC-Seq")

summaries <- map(files, summarize.by.bed, bed.file <- repmasker.file, drop.name.col=F)
df.fg <- Reduce(function(...) merge(..., all=TRUE), summaries)
rownames(df.fg) <- df.fg$name
df.fg$name <- NULL

out.file <- paste0(out.dir, 'ATAC_signal_at_repmasker.pdf')
colnames(df.fg) <- c('H3.3 WT', 'H3.3 KO')
summarized.heatmap(t(df.fg), out.file, "ATAC-Seq")
