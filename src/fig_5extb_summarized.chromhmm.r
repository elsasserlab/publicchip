library(ggplot2)
library(reshape2)
library(rtracklayer)
library(ggpubr)
library(purrr)
library(GenomicRanges)
library(pheatmap)

source('./code/bwutils.r')
source('./code/heatmap_utils.r')

chromhmm.file <- './data/bed/ChromHMM17.bed'
bw.dir <- './data/bw/'
out.dir <- './final_figures/05_ext/'

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
  'Deaton_2016/ES_H33_00h_rep1.bw',
  'Deaton_2016/ES_H33_03h_rep1.bw',
  'Deaton_2016/ES_H33_06h_rep1.bw',
  'Deaton_2016/ES_H33_12h_rep1.bw'
))

files.inp <- paste0(bw.dir, c(
  'Deaton_2016/ES_H33_inp_rep1.bw'
))

summaries <- map(files, summarize.by.bed, bed.file <- chromhmm.file, drop.name.col=F)
df.fg <- Reduce(function(...) merge(..., all=TRUE), summaries)
rownames(df.fg) <- df.fg$name
df.fg$name <- NULL

input <- summarize.by.bed(files.inp[[1]], chromhmm.file)
df.logfc <- log2(sweep(df.fg, 1, input$ES_H33_inp_rep1.bw, `/`))


colnames(df.logfc) <- c('0h','3h', '6h', '12h')
out.file <- paste0(out.dir, 'ES_SNAP_logfc_cov_enrichment_rep1.pdf')
summarized.heatmap(t(df.logfc[levels.order, ]), out.file, "ESC H3.3-SNAP Time-ChIP (Replicate 1)")

# Panel 2: ESC H3.3 SNAP time log2 fold change normalized to input replicate 2
files <- paste0(bw.dir, c(
  'Deaton_2016/ES_H33_00h_rep2.bw',
  'Deaton_2016/ES_H33_03h_rep2.bw',
  'Deaton_2016/ES_H33_06h_rep2.bw',
  'Deaton_2016/ES_H33_12h_rep2.bw'
))

files.inp <- paste0(bw.dir, c(
  'Deaton_2016/ES_H33_inp_rep2.bw'
))

summaries <- map(files, summarize.by.bed, bed.file <- chromhmm.file, drop.name.col=F)
df.fg <- Reduce(function(...) merge(..., all=TRUE), summaries)
rownames(df.fg) <- df.fg$name
df.fg$name <- NULL

input <- summarize.by.bed(files.inp[[1]], chromhmm.file)
df.logfc <- log2(sweep(df.fg, 1, input$ES_H33_inp_rep2.bw, `/`))


colnames(df.logfc) <- c('0h','3h', '6h', '12h')
out.file <- paste0(out.dir, 'ES_SNAP_logfc_cov_enrichment_rep2.pdf')
summarized.heatmap(t(df.logfc[levels.order, ]), out.file, "ESC H3.3-SNAP Time-ChIP (Replicate 2)")


# Panel 2: ESC H3.3 SNAP time log2 fold change normalized to input replicate 3
files <- paste0(bw.dir, c(
  'Deaton_2016/ES_H33_00h_rep3.bw',
  'Deaton_2016/ES_H33_03h_rep3.bw',
  'Deaton_2016/ES_H33_06h_rep3.bw',
  'Deaton_2016/ES_H33_12h_rep3.bw'
))

files.inp <- paste0(bw.dir, c(
  'Deaton_2016/ES_H33_inp_rep2.bw'
))

summaries <- map(files, summarize.by.bed, bed.file <- chromhmm.file, drop.name.col=F)
df.fg <- Reduce(function(...) merge(..., all=TRUE), summaries)
rownames(df.fg) <- df.fg$name
df.fg$name <- NULL

input <- summarize.by.bed(files.inp[[1]], chromhmm.file)
df.logfc <- log2(sweep(df.fg, 1, input$ES_H33_inp_rep2.bw, `/`))


colnames(df.logfc) <- c('0h','3h', '6h', '12h')
out.file <- paste0(out.dir, 'ES_SNAP_logfc_cov_enrichment_rep3.pdf')
summarized.heatmap(t(df.logfc[levels.order, ]), out.file, "ESC H3.3-SNAP Time-ChIP (Replicate 3)")

