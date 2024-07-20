### load packages
library(DSS)
library(dmrseq)
library(ChIPseeker)
library(BiocParallel)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(DOSE)
library(GenomicRanges)
library(fdrtool)
library(ArchR)
library(ggplot2)
library(ggfortify)
library(ggsci)
library(methylCC)
library(RColorBrewer)
library(liftOver)
library(rtracklayer)
library(DOSE)
library(DMRichR)
library(viridis)
library(enrichplot)
library(wesanderson)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
mmAnno <- getAnnot("mm10")

### read in data
infile <- list.files(pattern="cov$")

### split up by groups
F0_DDI <- 1:3
F0_FA <- 4:6
F1_DDI <- 7:9
F1_FA <- 10:12
F2_DDI <- 13:15
F2_FA <- 16:18
F3_DDI <- 19:21
F3_FA <- 22:24
F4_Cross <- 25:27
F4_DDI <- 28:30
F4_FA <- 31:33
F4_OCross <- 34:36
F4_RCross <- 37:39

### generate matrices
bs <- read.bismark(files = infile,rmZeroCov=TRUE,strandCollapse=T,verbose=T)
x <- colnames(bs)
x <- gsub(".CpG_report.merged_CpG_evidence.cov","",x)
x <- gsub("Sperm_","",x)
colnames(bs) <- x
pData(bs)$Condition <- rep(c("DDI","FA80","DDI","FA80","DDI","FA80","DDI","FA80","CROSS","DDI","FA80","Outcross","Reverse Outcross"),c(3,3,3,3,3,3,3,3,3,3,3,3,3))
pData(bs)$Generation <- rep(c("F0","F1","F2","F3","F4"),c(6,6,6,6,15))
