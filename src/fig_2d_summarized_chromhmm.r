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
out.dir <- './final_figures/02/'

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

# Panel 1: ESC H3.3 SNAP time log2 fold change normalized to input
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
out.file <- paste0(out.dir, 'ES_SNAP_logfc_cov_enrichment.pdf')
summarized.heatmap(t(df.logfc[levels.order, ]), out.file, "ESC H3.3-SNAP Time-ChIP")

# Panel 2: ES_H33HA_tet_ON
files <- paste0(bw.dir, c(
  'Yildirim_2014/2014_Rando_H33_P3.mm9.bw',
  'Yildirim_2014/2014_Rando_H33_P6.mm9.bw'
))

summaries <- map(files, summarize.by.bed, bed.file <- chromhmm.file, drop.name.col=F)
df.fg <- Reduce(function(...) merge(..., by='name', all=TRUE), summaries)
rownames(df.fg) <- df.fg$name
df.fg$name <- NULL
out.file <- paste0(out.dir, 'ES_H33HA_tet_ON_cov.pdf')
summarized.heatmap(t(df.fg[levels.order, ]), out.file, "ESC H3.3-HA Tet-ON")

# Panel 3: H3K9me3
files <- paste0(bw.dir, c(
   'Shi_2019/2019_Shi_ESC1_WT_H3K9me3.bw'
))

out.file <- paste0(out.dir, 'ESC_H3K9me3.pdf')
summaries <- summarize.by.bed(files, chromhmm.file)
summarized.heatmap(t(summaries[levels.order, ]), out.file, "ESC H3K9me3")
