#################################
### Goal: Calculate TPM from TMM
### Author: Janay Fox
### R script
###############################

## Set up ##
library(edgeR)

#read in data 
BN.genecounts <- read.table("./data/totalExpr/BN_bf_new_sal.gene.counts.matrix")
BN.genecounts <- round(BN.genecounts)
BA.genecounts <- read.table("./data/totalExpr/BA_bf_new_sal.gene.counts.matrix")
BA.genecounts <- round(BA.genecounts)

BN.lengths <- read.table("./data/totalExpr/BN.Trinity.gene_lengths.txt")
BA.lengths <- read.table("./data/totalExpr/BA.Trinity.gene_lengths.txt")

#read into edger
BN.dge <- DGEList(counts = BN.genecounts, genes = data.frame(Length=BN.lengths))
BA.dge <- DGEList(counts = BA.genecounts, genes = data.frame(Length=BA.lengths))

#calculate TMM
BN.dge.TMM <- calcNormFactors(BN.dge)
BA.dge.TMM <- calcNormFactors(BA.dge)

#calculate rpkm
BN.dge.TMM.rpkm <- rpkm(BN.dge.TMM, BN.dge.TMM$genes$Length.V2)
BA.dge.TMM.rpkm <- rpkm(BA.dge.TMM, BA.dge.TMM$genes$Length.V2)

#calculate tpm
BN.dge.TMM.tpm <- t(t(BN.dge.TMM.rpkm) / colSums(BN.dge.TMM.rpkm)) * 1e6
BA.dge.TMM.tpm <- t(t(BA.dge.TMM.rpkm) / colSums(BA.dge.TMM.rpkm)) * 1e6

#save tpms 
write.table(BN.dge.TMM.tpm, file = "./data/totalExpr/BN.TMM.TPM.matrix")
write.table(BA.dge.TMM.tpm, file = "./data/totalExpr/BA.TMM.TPM.matrix")

write.table(BN.dge.TMM.rpkm, file = "./data/totalExpr/BN.TMM.RPKM.matrix")
write.table(BA.dge.TMM.rpkm, file = "./data/totalExpr/BA.TMM.RPKM.matrix")
