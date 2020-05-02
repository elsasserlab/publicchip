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
out.dir <- './final_figures/06_ext/'

repmasker.file <- './data/bed/RepMasker_lt200bp_all.bed'

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

summaries <- map(files, summarize.by.bed, bed.file <- repmasker.file, drop.name.col=F)
df.fg <- Reduce(function(...) merge(..., all=TRUE), summaries)
rownames(df.fg) <- df.fg$name
df.fg$name <- NULL

input <- summarize.by.bed(files.inp[[1]], repmasker.file)
df.logfc <- log2(sweep(df.fg, 1, input$ES_H33_inp_rep1.bw, `/`))

out.file <- paste0(out.dir, 'ES_SNAP_logfc_cov_enrichment_repmasker.pdf')
colnames(df.logfc) <- c('0h','3h','6h','12h')
summarized.heatmap(t(df.logfc), out.file, "ESC H3.3-SNAP Time-ChIP")


# Panel 2: HP1alpha, beta, gamma
files <- paste0(bw.dir, c(
  'Ostapcuk_2018/ESC_HP1a_R1.bw',
  'Ostapcuk_2018/ESC_HP1b_R1.bw',
  'Ostapcuk_2018/ESC_HP1g_R1.bw'
))

summaries <- map(files, summarize.by.bed, bed.file <- repmasker.file, drop.name.col=F)
df.fg <- Reduce(function(...) merge(..., all=TRUE), summaries)
rownames(df.fg) <- df.fg$name
df.fg$name <- NULL

out.file <- paste0(out.dir, 'Hp1_cov_enrichment_repmasker.pdf')
colnames(df.fg) <- c('HP1a', 'HP1b', 'HP1g')
summarized.heatmap(t(df.fg), out.file, "")

# Panel 3: H3K9me3
files <- paste0(bw.dir, c(
  'Shi_2019/2019_Shi_ESC1_WT_H3K9me3.bw'
))

summaries <- map(files, summarize.by.bed, bed.file <- repmasker.file, drop.name.col=F)
df.fg <- Reduce(function(...) merge(..., all=TRUE), summaries)
rownames(df.fg) <- df.fg$name
df.fg$name <- NULL

out.file <- paste0(out.dir, 'H3K9me3_cov_enrichment_repmasker.pdf')
colnames(df.fg) <- c('H3K9me3')
summarized.heatmap(t(df.fg), out.file, "")

