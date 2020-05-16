library(ggpubr)

CountTable <- read.table("../data/CountTable.txt", sep="\t", header=TRUE)

CountTable$Counts <- CountTable$ESC_WT.bam
CountTable$ESC_WT.bam <- NULL
CountTable$Copies <- 1;
CountTable$Length <- CountTable$Stop-CountTable$Start;

CountTable <- CountTable[CountTable$Length<100000,]
CountTable <- CountTable[CountTable$Length>2000,]

CountTable$Stop <- NULL
CountTable$Start <- NULL
CountTable$RPK <- CountTable$Counts / CountTable$Length * 1000
CountTable$TPM <- CountTable$RPK / sum(CountTable$RPK) * 1000000

CountTable$lg10TPM <- log(CountTable$TPM,10)
CountTable$REF <- grepl("REF_", CountTable$Gene)

CountTable$Type <- NULL
CountTable$Type[grepl("REF", CountTable$Gene)] <- "Genes"
CountTable$Type[grepl("REP", CountTable$Gene)] <- "Other Repeats"
CountTable$Type[grepl("REP_IAPEz", CountTable$Gene)] <- "IAPEz"
CountTable$Type[grepl("ETn", CountTable$Gene)] <- "ETn"
CountTable$Type[grepl("REP_ETnERV2-int", CountTable$Gene)] <- "MusD"


ggviolin(CountTable, x = "Type", y = "lg10TPM", color = "Type", alpha = 0.7, draw_quantiles = 0.5, add = "jitter",  add.params = list(size = 0.1, alpha=0.5), ylim = c(-3,5)) +
  grids(axis="y")



Agg_CT <- aggregate(.~Gene, data=CountTable, FUN=sum)
Agg_CT$RPK <- Agg_CT$Counts / Agg_CT$Length * 1000
Agg_CT$RPK <- Agg_CT$Counts / Agg_CT$Copies
Agg_CT$TPM <- Agg_CT$RPK / sum(Agg_CT$RPK) * 1000000

Agg_CT$lg10TPM <- log(Agg_CT$TPM,10)
Agg_CT$REF <- grepl("REF_", Agg_CT$Gene)

write.table(Agg_CT, pipe("pbcopy"), sep="\t")

ggdensity(Agg_CT, x = "lg10TPM", y="..count..", color="REF", rug=TRUE, add = "mean", fill = "REF",
          palette = c("#00AFBB","#E7B800")) + 
  geom_vline(xintercept = Agg_CT[grepl("REP_IAPEz", Agg_CT$Gene),]$lg10TPM, color="black") +
  geom_vline(xintercept = Agg_CT[grepl("REP_ETnERV-int", Agg_CT$Gene),]$lg10TPM, color="green") +
  geom_vline(xintercept = Agg_CT[grepl("REP_ETnERV2-int", Agg_CT$Gene),]$lg10TPM, color="darkgreen") +
  geom_vline(xintercept = Agg_CT[grepl("REP_ETnERV3-int", Agg_CT$Gene),]$lg10TPM, color="turquoise")


CountTable$lg10TPM <- log(CountTable$TPM,10)
CountTable$REF <- grepl("REF_", CountTable$Gene)

CountTable$Type <- NULL
CountTable$Type[grepl("REF", CountTable$Gene)] <- "Genes"
CountTable$Type[grepl("REP", CountTable$Gene)] <- "Other Repeats"
CountTable$Type[grepl("REP_IAPEz", CountTable$Gene)] <- "IAPEz"
CountTable$Type[grepl("ETn", CountTable$Gene)] <- "ETn"
CountTable$Type[grepl("REP_ETnERV2-int", CountTable$Gene)] <- "MusD"

write.table(CountTable, pipe("pbcopy"), sep="\t")

ggviolin(CountTable, x = "Type", y = "lg10TPM", color = "Type", alpha = 0.7, draw_quantiles = 0.5, add = "jitter",  add.params = list(size = 0.1, alpha=0.5), ylim = c(-3,5)) +
  grids(axis="y")
