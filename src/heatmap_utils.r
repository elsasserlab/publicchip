library(pheatmap)
library(RColorBrewer)
library(stringr)

compute.lims <- function(mat) {
  maxmat <- mat
  maxmat[is.na(maxmat)] <- -Inf
  maxval <- max(maxmat)
  
  minmat <- mat
  minmat[is.na(minmat)] <- Inf
  minval <- min(minmat)
  
  nsteps <- 21.0
  
  breaklim <- round(abs(max(abs(c(minval, maxval)))), digits=2)
  
  # Compute a reasonable amount of steps?
  stepsize <- 2*breaklim / nsteps

  breakslist <- seq(-breaklim, +breaklim, by=stepsize)
}

summarized.heatmap <- function(values, filename, title, breakslist, clean.names=T, log=T) {
  bcolor <- "white"
  
  breakslist <- compute.lims(values)
  palette <- colorRampPalette(brewer.pal(n = 7, name = "Reds"))(length(breakslist))
  if (log == TRUE) {
    palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breakslist))
  }
  
  cellsize <- 35
  
  rowlabels <- rownames(values)
  collabels <- colnames(values)
  
  if (clean.names == TRUE) {
    rowlabels <- sapply(rowlabels, str_replace_all, '[._]', ' ')
    collabels <- sapply(collabels, str_replace_all, '[._]', ' ')
  }
  
  pheatmap(values,
           filename=filename,
           main=title, cluster_rows=F, cluster_cols=F,
           cellwidth=cellsize, cellheight=cellsize,
           border_color=bcolor, breaks=breakslist,
           labels_row=rowlabels, labels_col=collabels,
           color=palette,
           display_numbers=T)
}
