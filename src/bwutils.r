
library(rtracklayer)
library(GenomicRanges)

intersect.bw.with.bed <- function(bw.file, bed.file, stat='mean', keep.name=T) {
  bed <- import(bed.file)
  bw <- BigWigFile(bw.file)
  # summarize over bed regions using stat
  scored.gr <- unlist(summary(bw, bed, type=stat))
  if (keep.name == TRUE) {
    scored.gr$name <- bed$name
  }
  scored.gr
}

summarize.by.bed <- function (bw.file, bed.file, stat='mean', drop.name.col=T) {
  values <- intersect.bw.with.bed(bw.file, bed.file, keep.name=T)
  df <- data.frame(mcols(values))
  result <- aggregate(.~name, data=df, mean)
  colnames(result) <- c('name', basename(bw.file))
  if (drop.name.col == TRUE) {
    rownames(result) <- result$name
    result$name <- NULL  
  }
  result
}
