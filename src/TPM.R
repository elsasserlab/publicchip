
TPM <- read.table("TPM.txt", sep="\t", header = TRUE)

TPM <- TPM[TPM$ESC_WT.bam>10,]

library(ggpubr)

Genes_LG2TPM <- log(TPM[grepl("REF_", TPM$X),]$TPM,2)
Rep_LG2TPM <- log(TPM[grepl("REP_", TPM$X),]$TPM,2)


TPM$REF <- grepl("REF_", TPM$X)
TPM$LG2TPM <- log(TPM$TPM,2)

TPM$LG10TPM <- log(TPM$TPM,10)

ggdensity(TPM, x = "LG10TPM", y="..count..", color="REF", rug=TRUE, add = "mean", fill = "REF",
          palette = c("#00AFBB","#E7B800")) + 
          geom_vline(xintercept = TPM[grepl("REP_IAPEz", TPM$X),]$LG10TPM) +
          geom_vline(xintercept = TPM[grepl("REP_ETnERV-int", TPM$X),]$LG10TPM)
  

abline(v=0)

TPM[grepl("REP_IAP", TPM$X),]

TPM[grepl("REP_ETn", TPM$X),]
