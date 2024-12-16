######################################################################
### Goal: Run goSeq analysis and plot results
### Author: Janay Fox
### R script
#####################################################################

## Set up ##
library(goseq)
library(GO.db)
library(qvalue)
library(ggplot2)
library(dplyr)
library(tidyverse)

#read in list of all DEGs 
BN.allDEG <- read.table("./data/DEG/BN/salmon/diffExpr.P0.01_C2.matrix")
BA.allDEG <- read.table("./data/DEG/BA/salmon/diffExpr.P0.01_C2.matrix")

#read in lists of genes to test 
BN.StFvsSwF.up.factor = read.table("./data/goSeqData/BN_StF.vs.SwF_up_factor.txt", row.names=2, header=F)
colnames(BN.StFvsSwF.up.factor) = c('type')
BN.StFvsSwF.up.genes = rownames(BN.StFvsSwF.up.factor)

BN.StFvsSwF.down.factor = read.table("./data/goSeqData/BN_StF.vs.SwF_down_factor.txt", row.names=2, header=F)
colnames(BN.StFvsSwF.down.factor) = c('type')
BN.StFvsSwF.down.genes = rownames(BN.StFvsSwF.down.factor)

BN.SwFvsSwP.up.factor = read.table("./data/goSeqData/BN_SwF.vs.SwP_up_factor.txt", row.names=2, header=F)
colnames(BN.SwFvsSwP.up.factor) = c('type')
BN.SwFvsSwP.up.genes= rownames(BN.SwFvsSwP.up.factor)

BN.SwFvsSwP.down.factor = read.table("./data/goSeqData/BN_SwF.vs.SwP_down_factor.txt", row.names=2, header=F)
colnames(BN.SwFvsSwP.down.factor) = c('type')
BN.SwFvsSwP.down.genes= rownames(BN.SwFvsSwP.down.factor)

BN.vennGenes.factor = read.table("./data/goSeqData/BN_venn_factor.txt", row.names=2, header=F)
colnames(BN.vennGenes.factor) = c('type')
BN.vennGenes.genes = rownames(BN.vennGenes.factor)

BN.mfuzz.cluster5.factor = read.table("./data/goSeqData/BN_mfuzz_cluster5_factor.txt", row.names=2, header=F)
colnames(BN.mfuzz.cluster5.factor) = c('type')
BN.mfuzz.cluster5.genes = rownames(BN.mfuzz.cluster5.factor)

BN.mfuzz.cluster6.factor = read.table("./data/goSeqData/BN_mfuzz_cluster6_factor.txt", row.names=2, header=F)
colnames(BN.mfuzz.cluster6.factor) = c('type')
BN.mfuzz.cluster6.genes = rownames(BN.mfuzz.cluster6.factor)

BN.k.cluster2.factor = read.table("./data/goSeqData/BN_k_cluster2_factor.txt", row.names=2, header=F)
colnames(BN.k.cluster2.factor) = c('type')
BN.k.cluster2.genes = rownames(BN.k.cluster2.factor)

BN.30perc.cluster7.factor = read.table("./data/goSeqData/BN_30perc_cluster7_factor.txt", row.names=2, header=F)
colnames(BN.30perc.cluster7.factor) = c('type')
BN.30perc.cluster7.genes = rownames(BN.30perc.cluster7.factor)

BN.sameDir.up.factor = read.table("./data/goSeqData/BN_same_direction_up_factor.txt", row.names=2, header=F)
colnames(BN.sameDir.up.factor) = c('type')
BN.sameDir.up.genes = rownames(BN.sameDir.up.factor)

BN.sameDir.down.factor = read.table("./data/goSeqData/BN_same_direction_down_factor.txt", row.names=2, header=F)
colnames(BN.sameDir.down.factor) = c('type')
BN.sameDir.down.genes = rownames(BN.sameDir.down.factor)

BA.StFvsSwF.up.factor = read.table("./data/goSeqData/BA_StF.vs.SwF_up_factor.txt", row.names=2, header=F)
colnames(BA.StFvsSwF.up.factor) = c('type')
BA.StFvsSwF.up.genes = rownames(BA.StFvsSwF.up.factor)

BA.StFvsSwF.down.factor = read.table("./data/goSeqData/BA_StF.vs.SwF_down_factor.txt", row.names=2, header=F)
colnames(BA.StFvsSwF.down.factor) = c('type')
BA.StFvsSwF.down.genes = rownames(BA.StFvsSwF.down.factor)

BA.SwFvsSwP.up.factor = read.table("./data/goSeqData/BA_SwF.vs.SwP_up_factor.txt", row.names=2, header=F)
colnames(BA.SwFvsSwP.up.factor) = c('type')
BA.SwFvsSwP.up.genes = rownames(BA.SwFvsSwP.up.factor)

BA.SwFvsSwP.down.factor = read.table("./data/goSeqData/BA_SwF.vs.SwP_down_factor.txt", row.names=2, header=F)
colnames(BA.SwFvsSwP.down.factor) = c('type')
BA.SwFvsSwP.down.genes = rownames(BA.SwFvsSwP.down.factor)

BA.vennGenes.factor = read.table("./data/goSeqData/BA_venn_factor.txt", row.names=2, header=F)
colnames(BA.vennGenes.factor) = c('type')
BA.vennGenes.genes = rownames(BA.vennGenes.factor)

BA.mfuzz.cluster5.factor = read.table("./data/goSeqData/BA_mfuzz_cluster5_factor.txt", row.names=2, header=F)
colnames(BA.mfuzz.cluster5.factor) = c('type')
BA.mfuzz.cluster5.genes = rownames(BA.mfuzz.cluster5.factor)

BA.k.cluster3.factor = read.table("./data/goSeqData/BA_k_cluster3_factor.txt", row.names=2, header=F)
colnames(BA.k.cluster3.factor) = c('type')
BA.k.cluster3.genes = rownames(BA.k.cluster3.factor)

BA.30perc.cluster5.factor = read.table("./data/goSeqData/BA_30perc_cluster5_factor.txt", row.names=2, header=F)
colnames(BA.30perc.cluster5.factor) = c('type')
BA.30perc.cluster5.genes = rownames(BA.30perc.cluster5.factor)

BA.sameDir.up.factor = read.table("./data/goSeqData/BA_same_direction_up_factor.txt", row.names=2, header=F)
colnames(BA.sameDir.up.factor) = c('type')
BA.sameDir.up.genes = rownames(BA.sameDir.up.factor)

BA.sameDir.down.factor = read.table("./data/goSeqData/BA_same_direction_down_factor.txt", row.names=2, header=F)
colnames(BA.sameDir.down.factor) = c('type')
BA.sameDir.down.genes = rownames(BA.sameDir.down.factor)

BA.snp.factor = read.table("./data/goSeqData/EA_SNPs.txt", header=F)
colnames(BA.snp.factor) = c('type', 'gene')
BA.snp.factor.genes = unique(BA.snp.factor$gene)

BN.snp.factor = read.table("./data/goSeqData/EN_SNPs.txt", header=F)
colnames(BN.snp.factor) = c('type', 'gene')
BN.snp.factor.genes = unique(BN.snp.factor$gene)

#read in gene lengths
BN.gene.lengths = read.table("./data/goSeqData/BN.Trinity.gene_lengths.txt", header=T, row.names=1, com='')
BN.gene.lengths = as.matrix(BN.gene.lengths[,1,drop=F])

BA.gene.lengths = read.table("./data/goSeqData/BA.Trinity.gene_lengths.txt", header=T, row.names=1, com='')
BA.gene.lengths = as.matrix(BA.gene.lengths[,1,drop=F])

#read in background gene lists to test against 
BN.background = read.table("./data/goSeqData/BN_background.txt", header=F, row.names = 1)
BN.background.gene_ids = rownames(BN.background)

BN.background.snps = read.table("./data/goSeqData/EN_background_SNPs.txt", header=F)
BN.background.snps.gene_ids = unique(BN.background.snps$V1)

BA.background = read.table("./data/goSeqData/BA_background.txt", header=F, row.names = 1)
BA.background.gene_ids = rownames(BA.background)

BA.background.snps = read.table("./data/goSeqData/EA_background_SNPs.txt", header=F)
BA.background.snps.gene_ids = unique(BA.background.snps$V1)

# parse GO assignments
BN.GO.info = read.table("./data/goSeqData/BN_go_annotations.txt", header=F, row.names=1,stringsAsFactors=F)
BN.GO.info[is.na(BN.GO.info)] <- "NULL"
BN.GO.info.listed = apply(BN.GO.info, 1, function(x) unlist(strsplit(x,',')))
names(BN.GO.info.listed) = rownames(BN.GO.info)

BA.GO.info = read.table("./data/goSeqData/BA_go_annotations.txt", header=F, row.names=1,stringsAsFactors=F)
BA.GO.info[is.na(BA.GO.info)] <- "NULL"
BA.GO.info.listed = apply(BA.GO.info, 1, function(x) unlist(strsplit(x,',')))
names(BA.GO.info.listed) = rownames(BA.GO.info)

#get gene lengths for background genes only 
BN.background.gene.lengths = BN.gene.lengths[BN.background.gene_ids,]
BA.background.gene.lengths = BA.gene.lengths[BA.background.gene_ids,]

BN.background.snp.gene.lengths = BN.gene.lengths[BN.background.snps.gene_ids,]
BA.background.snp.gene.lengths = BA.gene.lengths[BA.background.snps.gene_ids,]

#get go terms for background genes only 
BN.GO.info.listed.snp = BN.GO.info.listed[ names(BN.GO.info.listed) %in% BN.background.snps.gene_ids ]
BN.GO.info.listed = BN.GO.info.listed[ names(BN.GO.info.listed) %in% BN.background.gene_ids ]

BA.GO.info.listed.snp = BA.GO.info.listed[ names(BA.GO.info.listed) %in% BA.background.snps.gene_ids ]
BA.GO.info.listed = BA.GO.info.listed[ names(BA.GO.info.listed) %in% BA.background.gene_ids ]

## GOseq analysis ##
#convert to vector suitable for use with goseq
BN.allDEG.genes.vec = as.integer(BN.background.gene_ids %in% rownames(BN.allDEG))
BN.StFvsSwF.up.vec = as.integer(BN.background.gene_ids %in% rownames(BN.StFvsSwF.up.factor))
BN.StFvsSwF.down.vec = as.integer(BN.background.gene_ids %in% rownames(BN.StFvsSwF.down.factor))
BN.SwFvsSwP.up.vec = as.integer(BN.background.gene_ids %in% rownames(BN.SwFvsSwP.up.factor))
BN.SwFvsSwP.down.vec = as.integer(BN.background.gene_ids %in% rownames(BN.SwFvsSwP.down.factor))
BN.vennGenes.vec = as.integer(BN.background.gene_ids %in% rownames(BN.vennGenes.factor))
BN.mfuzz.cluster5.vec = as.integer(BN.background.gene_ids %in% rownames(BN.mfuzz.cluster5.factor))
BN.mfuzz.cluster6.vec = as.integer(BN.background.gene_ids %in% rownames(BN.mfuzz.cluster6.factor))
BN.k.cluster2.vec = as.integer(BN.background.gene_ids %in% rownames(BN.k.cluster2.factor))
BN.30perc.cluster7.vec = as.integer(BN.background.gene_ids %in% rownames(BN.30perc.cluster7.factor))
BN.sameDir.up.vec = as.integer(BN.background.gene_ids %in% rownames(BN.sameDir.up.factor))
BN.sameDir.down.vec = as.integer(BN.background.gene_ids %in% rownames(BN.sameDir.down.factor))
BN.snp.vec = as.integer(BN.background.snps.gene_ids %in% BN.snp.factor$gene)

BA.allDEG.genes.vec = as.integer(BA.background.gene_ids %in% rownames(BA.allDEG))
BA.StFvsSwF.up.vec = as.integer(BA.background.gene_ids %in% rownames(BA.StFvsSwF.up.factor))
BA.StFvsSwF.down.vec = as.integer(BA.background.gene_ids %in% rownames(BA.StFvsSwF.down.factor))
BA.SwFvsSwP.up.vec = as.integer(BA.background.gene_ids %in% rownames(BA.SwFvsSwP.up.factor))
BA.SwFvsSwP.down.vec = as.integer(BA.background.gene_ids %in% rownames(BA.SwFvsSwP.down.factor))
BA.vennGenes.vec = as.integer(BA.background.gene_ids %in% rownames(BA.vennGenes.factor))
BA.mfuzz.cluster5.vec = as.integer(BA.background.gene_ids %in% rownames(BA.mfuzz.cluster5.factor))
BA.k.cluster3.vec = as.integer(BA.background.gene_ids %in% rownames(BA.k.cluster3.factor))
BA.30perc.cluster5.vec = as.integer(BA.background.gene_ids %in% rownames(BA.30perc.cluster5.factor))
BA.sameDir.up.vec = as.integer(BA.background.gene_ids %in% rownames(BA.sameDir.up.factor))
BA.sameDir.down.vec = as.integer(BA.background.gene_ids %in% rownames(BA.sameDir.down.factor))
BA.snp.vec = as.integer(BA.background.snps.gene_ids %in% BA.snp.factor$gene)

#fit the probabilty weight function (PWF) on all DE features 
BN.pwf=nullp(BN.allDEG.genes.vec, bias.data=BN.background.gene.lengths)
rownames(BN.pwf) = BN.background.gene_ids

BN.pwf.snp=nullp(BN.snp.vec, bias.data=BN.background.snp.gene.lengths)
rownames(BN.pwf.snp) = BN.background.snps.gene_ids

BA.pwf=nullp(BA.allDEG.genes.vec, bias.data=BA.background.gene.lengths)
rownames(BA.pwf) = BA.background.gene_ids

BA.pwf.snp=nullp(BA.snp.vec, bias.data=BA.background.snp.gene.lengths)
rownames(BA.pwf.snp) = BA.background.snps.gene_ids

#perform functional enrichment testing for each category of DEGs
BN.pwf$DEgenes = BN.StFvsSwF.up.vec
BN.StFvsSwF.up.res = goseq(BN.pwf, gene2cat=BN.GO.info.listed)

BN.pwf$DEgenes = BN.StFvsSwF.down.vec
BN.StFvsSwF.down.res = goseq(BN.pwf, gene2cat=BN.GO.info.listed)

BN.pwf$DEgenes = BN.SwFvsSwP.up.vec
BN.SwFvsSwP.up.res = goseq(BN.pwf, gene2cat=BN.GO.info.listed)

BN.pwf$DEgenes = BN.SwFvsSwP.down.vec
BN.SwFvsSwP.down.res = goseq(BN.pwf, gene2cat=BN.GO.info.listed)

BN.pwf$DEgenes = BN.vennGenes.vec
BN.vennGenes.res = goseq(BN.pwf, gene2cat=BN.GO.info.listed)

BN.pwf$DEgenes = BN.mfuzz.cluster5.vec
BN.mfuzz.cluster5.res = goseq(BN.pwf, gene2cat=BN.GO.info.listed)

BN.pwf$DEgenes = BN.mfuzz.cluster6.vec
BN.mfuzz.cluster6.res = goseq(BN.pwf, gene2cat=BN.GO.info.listed)

BN.pwf$DEgenes = BN.k.cluster2.vec
BN.k.cluster2.res = goseq(BN.pwf, gene2cat=BN.GO.info.listed)

BN.pwf$DEgenes = BN.30perc.cluster7.vec
BN.30perc.cluster7.res = goseq(BN.pwf, gene2cat=BN.GO.info.listed)

BN.pwf$DEgenes = BN.sameDir.up.vec
BN.sameDir.up.res = goseq(BN.pwf, gene2cat=BN.GO.info.listed)

BN.pwf$DEgenes = BN.sameDir.down.vec
BN.sameDir.down.res = goseq(BN.pwf, gene2cat=BN.GO.info.listed)

BA.pwf$DEgenes = BA.StFvsSwF.up.vec
BA.StFvsSwF.up.res = goseq(BA.pwf, gene2cat=BA.GO.info.listed)

BA.pwf$DEgenes = BA.StFvsSwF.down.vec
BA.StFvsSwF.down.res = goseq(BA.pwf, gene2cat=BA.GO.info.listed)

BA.pwf$DEgenes = BA.SwFvsSwP.up.vec
BA.SwFvsSwP.up.res = goseq(BA.pwf, gene2cat=BA.GO.info.listed)

BA.pwf$DEgenes = BA.SwFvsSwP.down.vec
BA.SwFvsSwP.down.res = goseq(BA.pwf, gene2cat=BA.GO.info.listed)

BA.pwf$DEgenes = BA.vennGenes.vec
BA.vennGenes.res = goseq(BA.pwf, gene2cat=BA.GO.info.listed)

BA.pwf$DEgenes = BA.mfuzz.cluster5.vec
BA.mfuzz.cluster5.res = goseq(BA.pwf, gene2cat=BA.GO.info.listed)

BA.pwf$DEgenes = BA.k.cluster3.vec
BA.k.cluster3.res = goseq(BA.pwf, gene2cat=BA.GO.info.listed)

BA.pwf$DEgenes = BA.30perc.cluster5.vec
BA.30perc.cluster5.res = goseq(BA.pwf, gene2cat=BA.GO.info.listed)

BA.pwf$DEgenes = BA.sameDir.up.vec
BA.sameDir.up.res = goseq(BN.pwf, gene2cat=BA.GO.info.listed)

BA.pwf$DEgenes = BA.sameDir.down.vec
BA.sameDir.down.res = goseq(BN.pwf, gene2cat=BA.GO.info.listed)

#for snps 
BN.pwf.snp$DEgenes = BN.snp.vec
BN.snp.res = goseq(BN.pwf.snp, gene2cat = BN.GO.info.listed.snp)

BA.pwf.snp$DEgenes = BA.snp.vec
BA.snp.res = goseq(BA.pwf.snp, gene2cat = BA.GO.info.listed.snp)

#correct p-values to q values 
q.val.correction <- function(go.res) {
                            pvals = go.res$over_represented_pvalue
                            pvals[pvals > 1 - 1e-10] = 1 - 1e-10
                            q = qvalue(pvals)
                            go.res$over_represented_FDR = q$qvalues
                            return(go.res)
                          }

BN.StFvsSwF.up.res <- q.val.correction(BN.StFvsSwF.up.res)
BN.StFvsSwF.down.res <- q.val.correction(BN.StFvsSwF.down.res)
BN.SwFvsSwP.up.res <- q.val.correction(BN.SwFvsSwP.up.res)
BN.SwFvsSwP.down.res <- q.val.correction(BN.SwFvsSwP.down.res)
BN.vennGenes.res <- q.val.correction(BN.vennGenes.res)
BN.mfuzz.cluster5.res <- q.val.correction(BN.mfuzz.cluster5.res)
BN.mfuzz.cluster6.res <- q.val.correction(BN.mfuzz.cluster6.res)
BN.k.cluster2.res <- q.val.correction(BN.k.cluster2.res)
BN.30perc.cluster7.res <- q.val.correction(BN.30perc.cluster7.res)
BN.sameDir.up.res <- q.val.correction(BN.sameDir.up.res)
BN.sameDir.down.res <- q.val.correction(BN.sameDir.down.res)

BA.StFvsSwF.up.res <- q.val.correction(BA.StFvsSwF.up.res)
BA.StFvsSwF.down.res <- q.val.correction(BA.StFvsSwF.down.res)
BA.SwFvsSwP.up.res <- q.val.correction(BA.SwFvsSwP.up.res)
BA.SwFvsSwP.down.res <- q.val.correction(BA.SwFvsSwP.down.res)
BA.vennGenes.res <- q.val.correction(BA.vennGenes.res)
BA.mfuzz.cluster5.res <- q.val.correction(BA.mfuzz.cluster5.res)
BA.k.cluster3.res <- q.val.correction(BA.k.cluster3.res)
BA.30perc.cluster5.res <- q.val.correction(BA.30perc.cluster5.res)
BA.sameDir.up.res <- q.val.correction(BA.sameDir.up.res)
BA.sameDir.down.res <- q.val.correction(BA.sameDir.down.res)

BN.snp.res$over_represented_FDR <- p.adjust(BN.snp.res$over_represented_pvalue, method = "hommel")
BA.snp.res <- q.val.correction(BA.snp.res)

#extract significant results
BN.StFvsSwF.up.res.sig <- BN.StFvsSwF.up.res[BN.StFvsSwF.up.res$over_represented_FDR<=0.05,]
BN.StFvsSwF.down.res.sig <- BN.StFvsSwF.down.res[BN.StFvsSwF.down.res$over_represented_FDR<=0.05,]
BN.SwFvsSwP.up.res.sig <- BN.SwFvsSwP.up.res[BN.SwFvsSwP.up.res$over_represented_FDR<=0.05,]
BN.SwFvsSwP.down.res.sig <- BN.SwFvsSwP.down.res[BN.SwFvsSwP.down.res$over_represented_FDR<=0.05,]
BN.vennGenes.res.sig <- BN.vennGenes.res[BN.vennGenes.res$over_represented_FDR<=0.05,]
BN.mfuzz.cluster5.res.sig <- BN.mfuzz.cluster5.res[BN.mfuzz.cluster5.res$over_represented_FDR<=0.05,]
BN.mfuzz.cluster6.res.sig <- BN.mfuzz.cluster6.res[BN.mfuzz.cluster6.res$over_represented_FDR<=0.05,]
BN.k.cluster2.res.sig <- BN.k.cluster2.res[BN.k.cluster2.res$over_represented_FDR<=0.05,]
BN.30perc.cluster7.res.sig <- BN.30perc.cluster7.res[BN.30perc.cluster7.res$over_represented_FDR<=0.05,]
BN.sameDir.up.res.sig <- BN.sameDir.up.res[BN.sameDir.up.res$over_represented_FDR<=0.05,]
BN.sameDir.down.res.sig <- BN.sameDir.down.res[BN.sameDir.down.res$over_represented_FDR<=0.05,]

BA.StFvsSwF.up.res.sig <- BA.StFvsSwF.up.res[BA.StFvsSwF.up.res$over_represented_FDR<=0.05,]
BA.StFvsSwF.down.res.sig <- BA.StFvsSwF.down.res[BA.StFvsSwF.down.res$over_represented_FDR<=0.05,]
BA.SwFvsSwP.up.res.sig <- BA.SwFvsSwP.up.res[BA.SwFvsSwP.up.res$over_represented_FDR<=0.05,]
BA.SwFvsSwP.down.res.sig <- BA.SwFvsSwP.down.res[BA.SwFvsSwP.down.res$over_represented_FDR<=0.05,]
BA.vennGenes.res.sig <- BA.vennGenes.res[BA.vennGenes.res$over_represented_FDR<=0.05,]
BA.mfuzz.cluster5.res.sig <- BA.mfuzz.cluster5.res[BA.mfuzz.cluster5.res$over_represented_FDR<=0.05,]
BA.k.cluster3.res.sig <- BA.k.cluster3.res[BA.k.cluster3.res$over_represented_FDR<=0.05,]
BA.30perc.cluster5.res.sig <- BA.30perc.cluster5.res[BA.30perc.cluster5.res$over_represented_FDR<=0.05,]
BA.sameDir.up.res.sig <- BA.sameDir.up.res[BA.sameDir.up.res$over_represented_FDR<=0.05,]
BA.sameDir.down.res.sig <- BA.sameDir.down.res[BA.sameDir.down.res$over_represented_FDR<=0.05,]

BN.snp.res.sig <- BN.snp.res[BN.snp.res$over_represented_FDR<=0.05,]
BA.snp.res.sig <- BA.snp.res[BA.snp.res$over_represented_FDR<=0.05,]

#remove obsolete terms
BN.StFvsSwF.up.res.sig <- BN.StFvsSwF.up.res.sig[BN.StFvsSwF.up.res.sig$category != "GO:0044421",]
BN.StFvsSwF.up.res.sig <- BN.StFvsSwF.up.res.sig[BN.StFvsSwF.up.res.sig$category != "GO:0001871",]

BN.SwFvsSwP.down.res.sig <- BN.SwFvsSwP.down.res.sig[BN.SwFvsSwP.down.res.sig$category != "GO:0044421",]
BN.SwFvsSwP.down.res.sig <- BN.SwFvsSwP.down.res.sig[BN.SwFvsSwP.down.res.sig$category != "GO:0055114",]
BN.SwFvsSwP.down.res.sig <- BN.SwFvsSwP.down.res.sig[BN.SwFvsSwP.down.res.sig$category != "GO:0051704",]
BN.SwFvsSwP.down.res.sig <- BN.SwFvsSwP.down.res.sig[BN.SwFvsSwP.down.res.sig$category != "GO:0050880",]
BN.SwFvsSwP.down.res.sig <- BN.SwFvsSwP.down.res.sig[BN.SwFvsSwP.down.res.sig$category != "GO:0044707",]

BN.mfuzz.cluster5.res.sig <- BN.mfuzz.cluster5.res.sig[BN.mfuzz.cluster5.res.sig$category != "GO:0001871",]
BN.mfuzz.cluster5.res.sig <- BN.mfuzz.cluster5.res.sig[BN.mfuzz.cluster5.res.sig$category != "GO:0044421",]

BN.sameDir.up.res.sig <- BN.sameDir.up.res.sig[BN.sameDir.up.res.sig$category != "GO:0050880",]
BN.sameDir.up.res.sig <- BN.sameDir.up.res.sig[BN.sameDir.up.res.sig$category != "GO:0015077",]
BN.sameDir.up.res.sig <- BN.sameDir.up.res.sig[BN.sameDir.up.res.sig$category != "GO:0022892",]
BN.sameDir.up.res.sig <- BN.sameDir.up.res.sig[BN.sameDir.up.res.sig$category != "GO:0071804",]
BN.sameDir.up.res.sig <- BN.sameDir.up.res.sig[BN.sameDir.up.res.sig$category != "GO:0044421",]

BA.StFvsSwF.up.res.sig <- BA.StFvsSwF.up.res.sig[BA.StFvsSwF.up.res.sig$category != "GO:0044421",]

BA.SwFvsSwP.down.res.sig <- BA.SwFvsSwP.down.res.sig[BA.SwFvsSwP.down.res.sig$category != "GO:0030529",]
BA.SwFvsSwP.down.res.sig <- BA.SwFvsSwP.down.res.sig[BA.SwFvsSwP.down.res.sig$category != "GO:0044445",]
BA.SwFvsSwP.down.res.sig <- BA.SwFvsSwP.down.res.sig[BA.SwFvsSwP.down.res.sig$category != "GO:0044267",]

BA.mfuzz.cluster5.res.sig <- BA.mfuzz.cluster5.res.sig[BA.mfuzz.cluster5.res.sig$category != "GO:0044421",]

BA.30perc.cluster5.res.sig <- BA.30perc.cluster5.res.sig[BA.30perc.cluster5.res.sig$category != "GO:0044421",]
BA.30perc.cluster5.res.sig <- BA.30perc.cluster5.res.sig[BA.30perc.cluster5.res.sig$category != "GO:0035587",]
BA.30perc.cluster5.res.sig <- BA.30perc.cluster5.res.sig[BA.30perc.cluster5.res.sig$category != "GO:0051818",]
BA.30perc.cluster5.res.sig <- BA.30perc.cluster5.res.sig[BA.30perc.cluster5.res.sig$category != "GO:0051852",]
BA.30perc.cluster5.res.sig <- BA.30perc.cluster5.res.sig[BA.30perc.cluster5.res.sig$category != "GO:0051883",]
BA.30perc.cluster5.res.sig <- BA.30perc.cluster5.res.sig[BA.30perc.cluster5.res.sig$category != "GO:0051704",]

#organize go id -> list of genes
BN.GO.to.gene.list = list()
for (gene_id in intersect(names(BN.GO.info.listed), BN.background.gene_ids)) {
  BN.go.list = BN.GO.info.listed[[gene_id]]
  for (go_id in BN.go.list) {
    BN.GO.to.gene.list[[go_id]] = c(BN.GO.to.gene.list[[go_id]], gene_id)
  }
}

BA.GO.to.gene.list = list()
for (gene_id in intersect(names(BA.GO.info.listed), BA.background.gene_ids)) {
  BA.go.list = BA.GO.info.listed[[gene_id]]
  for (go_id in BA.go.list) {
    BA.GO.to.gene.list[[go_id]] = c(BA.GO.to.gene.list[[go_id]], gene_id)
  }
}

#add on gene id column to results 
BN.StFvsSwF.up.res.sig$gene_ids = do.call(rbind, lapply(BN.StFvsSwF.up.res.sig$category, function(x) { 
  gene_list = BN.GO.to.gene.list[[x]]
  gene_list = gene_list[gene_list %in% rownames(BN.StFvsSwF.up.factor)]
  paste(gene_list, collapse=', ');
}) )

BN.StFvsSwF.down.res.sig$gene_ids = do.call(rbind, lapply(BN.StFvsSwF.down.res.sig$category, function(x) { 
  gene_list = BN.GO.to.gene.list[[x]]
  gene_list = gene_list[gene_list %in% rownames(BN.StFvsSwF.down.factor)]
  paste(gene_list, collapse=', ');
}) )

BN.SwFvsSwP.up.res.sig$gene_ids = do.call(rbind, lapply(BN.SwFvsSwP.up.res.sig$category, function(x) { 
  gene_list = BN.GO.to.gene.list[[x]]
  gene_list = gene_list[gene_list %in% rownames(BN.SwFvsSwP.up.factor)]
  paste(gene_list, collapse=', ');
}) )

BN.SwFvsSwP.down.res.sig$gene_ids = do.call(rbind, lapply(BN.SwFvsSwP.down.res.sig$category, function(x) { 
  gene_list = BN.GO.to.gene.list[[x]]
  gene_list = gene_list[gene_list %in% rownames(BN.SwFvsSwP.down.factor)]
  paste(gene_list, collapse=', ');
}) )

BN.vennGenes.res.sig$gene_ids = do.call(rbind, lapply(BN.vennGenes.res.sig$category, function(x) { 
  gene_list = BN.GO.to.gene.list[[x]]
  gene_list = gene_list[gene_list %in% rownames(BN.vennGenes.factor)]
  paste(gene_list, collapse=', ');
}) )

BN.mfuzz.cluster5.res.sig$gene_ids = do.call(rbind, lapply(BN.mfuzz.cluster5.res.sig$category, function(x) { 
  gene_list = BN.GO.to.gene.list[[x]]
  gene_list = gene_list[gene_list %in% rownames(BN.mfuzz.cluster5.factor)]
  paste(gene_list, collapse=', ');
}) )

BN.mfuzz.cluster6.res.sig$gene_ids = do.call(rbind, lapply(BN.mfuzz.cluster6.res.sig$category, function(x) { 
  gene_list = BN.GO.to.gene.list[[x]]
  gene_list = gene_list[gene_list %in% rownames(BN.mfuzz.cluster6.factor)]
  paste(gene_list, collapse=', ');
}) )

BN.30perc.cluster7.res.sig$gene_ids = do.call(rbind, lapply(BN.30perc.cluster7.res.sig$category, function(x) { 
  gene_list = BN.GO.to.gene.list[[x]]
  gene_list = gene_list[gene_list %in% rownames(BN.30perc.cluster7.factor)]
  paste(gene_list, collapse=', ');
}) )

BN.sameDir.up.res.sig$gene_ids = do.call(rbind, lapply(BN.sameDir.up.res.sig$category, function(x) { 
  gene_list = BN.GO.to.gene.list[[x]]
  gene_list = gene_list[gene_list %in% rownames(BN.sameDir.up.factor)]
  paste(gene_list, collapse=', ');
}) )

BN.sameDir.down.res.sig$gene_ids = do.call(rbind, lapply(BN.sameDir.down.res.sig$category, function(x) { 
  gene_list = BN.GO.to.gene.list[[x]]
  gene_list = gene_list[gene_list %in% rownames(BN.sameDir.down.factor)]
  paste(gene_list, collapse=', ');
}) )

BA.StFvsSwF.up.res.sig$gene_ids = do.call(rbind, lapply(BA.StFvsSwF.up.res.sig$category, function(x) { 
  gene_list = BA.GO.to.gene.list[[x]]
  gene_list = gene_list[gene_list %in% rownames(BA.StFvsSwF.up.factor)]
  paste(gene_list, collapse=', ');
}) )

BA.StFvsSwF.down.res.sig$gene_ids = do.call(rbind, lapply(BA.StFvsSwF.down.res.sig$category, function(x) { 
  gene_list = BA.GO.to.gene.list[[x]]
  gene_list = gene_list[gene_list %in% rownames(BA.StFvsSwF.down.factor)]
  paste(gene_list, collapse=', ');
}) )

BA.SwFvsSwP.up.res.sig$gene_ids = do.call(rbind, lapply(BA.SwFvsSwP.up.res.sig$category, function(x) { 
  gene_list = BA.GO.to.gene.list[[x]]
  gene_list = gene_list[gene_list %in% rownames(BA.SwFvsSwP.up.factor)]
  paste(gene_list, collapse=', ');
}) )

BA.SwFvsSwP.down.res.sig$gene_ids = do.call(rbind, lapply(BA.SwFvsSwP.down.res.sig$category, function(x) { 
  gene_list = BA.GO.to.gene.list[[x]]
  gene_list = gene_list[gene_list %in% rownames(BA.SwFvsSwP.down.factor)]
  paste(gene_list, collapse=', ');
}) )

BA.mfuzz.cluster5.res.sig$gene_ids = do.call(rbind, lapply(BA.mfuzz.cluster5.res.sig$category, function(x) { 
  gene_list = BA.GO.to.gene.list[[x]]
  gene_list = gene_list[gene_list %in% rownames(BA.mfuzz.cluster5.factor)]
  paste(gene_list, collapse=', ');
}) )

BA.30perc.cluster5.res.sig$gene_ids = do.call(rbind, lapply(BA.30perc.cluster5.res.sig$category, function(x) { 
  gene_list = BA.GO.to.gene.list[[x]]
  gene_list = gene_list[gene_list %in% rownames(BA.30perc.cluster5.factor)]
  paste(gene_list, collapse=', ');
}) )

#save results
write.table(BN.StFvsSwF.up.res.sig, file="BN_StFvsSwF_up_sigGO.tsv", sep='	', quote=F, row.names=F)
write.table(BN.StFvsSwF.down.res.sig, file="BN_StFvsSwF_down_sigGO.tsv", sep='	', quote=F, row.names=F)
write.table(BN.SwFvsSwP.up.res.sig, file="BN_SwFvsSwP_up_sigGO.tsv", sep='	', quote=F, row.names=F)
write.table(BN.SwFvsSwP.down.res.sig, file="BN_SwFvsSwP_down_sigGO.tsv", sep='	', quote=F, row.names=F)
write.table(BN.vennGenes.res.sig, file="BN_vennGenes_sigGO.tsv", sep='	', quote=F, row.names=F)
write.table(BN.mfuzz.cluster5.res.sig, file="BN_mfuzz_cluster5_sigGO.tsv", sep='	', quote=F, row.names=F)
write.table(BN.mfuzz.cluster6.res.sig, file="BN_mfuzz_cluster6_sigGO.tsv", sep='	', quote=F, row.names=F)
#write.table(BN.k.cluster2.res.sig, file="BN_k_cluster2_sigGO.tsv", sep='	', quote=F, row.names=F) #empty
write.table(BN.30perc.cluster7.res.sig, file="BN_30perc.cluster7_sigGO.tsv", sep='	', quote=F, row.names=F)
write.table(BN.sameDir.up.res.sig, file="BN_sameDir_up_sigGO.tsv", sep='	', quote=F, row.names=F)
write.table(BN.sameDir.down.res.sig, file="BN_sameDir_down_sigGO.tsv", sep='	', quote=F, row.names=F)

write.table(BA.StFvsSwF.up.res.sig, file="BA_StFvsSwF_up_sigGO.tsv", sep='	', quote=F, row.names=F)
write.table(BA.StFvsSwF.down.res.sig, file="BA_StFvsSwF_down_sigGO.tsv", sep='	', quote=F, row.names=F)
write.table(BA.SwFvsSwP.up.res.sig, file="BA_SwFvsSwP_up_sigGO.tsv", sep='	', quote=F, row.names=F)
write.table(BA.SwFvsSwP.down.res.sig, file="BA_SwFvsSwP_down_sigGO.tsv", sep='	', quote=F, row.names=F)
#write.table(BA.vennGenes.res.sig, file="BA_vennGenes_sigGO.tsv", sep='	', quote=F, row.names=F)
write.table(BA.mfuzz.cluster5.res.sig, file="BA_mfuzz_cluster5_sigGO.tsv", sep='	', quote=F, row.names=F)
#write.table(BA.k.cluster3.res.sig, file="BA_k_cluster3_sigGO.tsv", sep='	', quote=F, row.names=F)
write.table(BA.30perc.cluster5.res.sig, file="BA_30perc.cluster5_sigGO.tsv", sep='	', quote=F, row.names=F)

#save a list for REVIGO 
write.table(BN.StFvsSwF.up.res.sig[,c(1,8)], file="BN_StFvsSwF_up_sigGO_revigo.tsv", sep='	', quote=F, row.names=F)
write.table(BN.StFvsSwF.down.res.sig[,c(1,8)], file="BN_StFvsSwF_down_sigGO_revigo.tsv", sep='	', quote=F, row.names=F)
write.table(BN.SwFvsSwP.up.res.sig[,c(1,8)], file="BN_SwFvsSwP_up_sigGO_revigo.tsv", sep='	', quote=F, row.names=F)
write.table(BN.SwFvsSwP.down.res.sig[,c(1,8)], file="BN_SwFvsSwP_down_sigGO_revigo.tsv", sep='	', quote=F, row.names=F)
write.table(BN.vennGenes.res.sig[,c(1,8)], file="BN_vennGenes_sigGO_revigo.tsv", sep='	', quote=F, row.names=F)
write.table(BN.mfuzz.cluster5.res.sig[,c(1,8)], file="BN_mfuzz_cluster5_sigGO_revigo.tsv", sep='	', quote=F, row.names=F)
write.table(BN.mfuzz.cluster6.res.sig[,c(1,8)], file="BN_mfuzz_cluster6_sigGO_revigo.tsv", sep='	', quote=F, row.names=F)
#write.table(BN.k.cluster2.res.sig[,c(1,8)], file="BN_k_cluster2_sigGO_revigo.tsv", sep='	', quote=F, row.names=F) #empty
write.table(BN.30perc.cluster7.res.sig[,c(1,8)], file="BN_30perc.cluster7_sigGO_revigo.tsv", sep='	', quote=F, row.names=F)
write.table(BN.sameDir.up.res.sig[,c(1,8)], file="BN_sameDir_up_sigGO_revigo.tsv", sep='	', quote=F, row.names=F)
write.table(BN.sameDir.down.res.sig[,c(1,8)], file="BN_sameDir_down_sigGO_revigo.tsv", sep='	', quote=F, row.names=F)

write.table(BA.StFvsSwF.up.res.sig[,c(1,8)], file="BA_StFvsSwF_up_sigGO_revigo.tsv", sep='	', quote=F, row.names=F)
write.table(BA.StFvsSwF.down.res.sig[,c(1,8)], file="BA_StFvsSwF_down_sigGO_revigo.tsv", sep='	', quote=F, row.names=F)
write.table(BA.SwFvsSwP.up.res.sig[,c(1,8)], file="BA_SwFvsSwP_up_sigGO_revigo.tsv", sep='	', quote=F, row.names=F)
write.table(BA.SwFvsSwP.down.res.sig[,c(1,8)], file="BA_SwFvsSwP_down_sigGO_revigo.tsv", sep='	', quote=F, row.names=F)
#write.table(BA.vennGenes.res.sig[,c(1,8)], file="BA_vennGenes_sigGO_revigo.tsv", sep='	', quote=F, row.names=F)
write.table(BA.mfuzz.cluster5.res.sig[,c(1,8)], file="BA_mfuzz_cluster5_sigGO_revigo.tsv", sep='	', quote=F, row.names=F)
#write.table(BA.k.cluster3.res.sig[,c(1,8)], file="BA_k_cluster3_sigGO_revigo.tsv", sep='	', quote=F, row.names=F)
write.table(BA.30perc.cluster5.res.sig[,c(1,8)], file="BA_30perc.cluster5_sigGO_revigo.tsv", sep='	', quote=F, row.names=F)

#extract the top 10 results and calculate hits percentage
BN.StFvsSwF.up.top <- BN.StFvsSwF.up.res.sig %>% slice_min(over_represented_FDR, n = 10)
BN.StFvsSwF.up.top <- BN.StFvsSwF.up.top %>% mutate(hitsPerc = numDEInCat * 100/numInCat)

BN.StFvsSwF.down.top <- BN.StFvsSwF.down.res.sig %>% slice_min(over_represented_FDR, n = 10)
BN.StFvsSwF.down.top <- BN.StFvsSwF.down.top %>% mutate(hitsPerc = numDEInCat * 100/numInCat)

BN.SwFvsSwP.up.top <- BN.SwFvsSwP.up.res.sig %>% slice_min(over_represented_FDR, n = 10)
BN.SwFvsSwP.up.top <- BN.SwFvsSwP.up.top %>% mutate(hitsPerc = numDEInCat * 100/numInCat)

BN.SwFvsSwP.down.top <- BN.SwFvsSwP.down.res.sig %>% slice_min(over_represented_FDR, n = 10)
BN.SwFvsSwP.down.top <- BN.SwFvsSwP.down.top %>% mutate(hitsPerc = numDEInCat * 100/numInCat)

BN.vennGenes.top <- BN.vennGenes.res.sig %>% slice_min(over_represented_FDR, n = 10)
BN.vennGenes.top <- BN.vennGenes.top %>% mutate(hitsPerc = numDEInCat * 100/numInCat)

BN.mfuzz.cluster5.top <- BN.mfuzz.cluster5.res.sig %>% slice_min(over_represented_FDR, n = 10)
BN.mfuzz.cluster5.top <- BN.mfuzz.cluster5.top %>% mutate(hitsPerc = numDEInCat * 100/numInCat)

BN.mfuzz.cluster6.top <- BN.mfuzz.cluster6.res.sig %>% slice_min(over_represented_FDR, n = 10)
BN.mfuzz.cluster6.top <- BN.mfuzz.cluster6.top %>% mutate(hitsPerc = numDEInCat * 100/numInCat)

BN.30perc.cluster7.top <- BN.30perc.cluster7.res.sig %>% slice_min(over_represented_FDR, n = 10)
BN.30perc.cluster7.top <- BN.30perc.cluster7.top %>% mutate(hitsPerc = numDEInCat * 100/numInCat)

BN.sameDir.up.top <- BN.sameDir.up.res.sig %>% slice_min(over_represented_FDR, n = 10)
BN.sameDir.up.top <- BN.sameDir.up.top %>% mutate(hitsPerc = numDEInCat * 100/numInCat)

BN.sameDir.down.top <- BN.sameDir.down.res.sig %>% slice_min(over_represented_FDR, n = 10)
BN.sameDir.down.top <- BN.sameDir.down.top %>% mutate(hitsPerc = numDEInCat * 100/numInCat)

BA.StFvsSwF.up.top <- BA.StFvsSwF.up.res.sig %>% slice_min(over_represented_FDR, n = 10)
BA.StFvsSwF.up.top <- BA.StFvsSwF.up.top %>% mutate(hitsPerc = numDEInCat * 100/numInCat)

BA.StFvsSwF.down.top <- BA.StFvsSwF.down.res.sig %>% slice_min(over_represented_FDR, n = 10)
BA.StFvsSwF.down.top <- BA.StFvsSwF.down.top %>% mutate(hitsPerc = numDEInCat * 100/numInCat)

BA.SwFvsSwP.up.top <- BA.SwFvsSwP.up.res.sig %>% slice_min(over_represented_FDR, n = 10)
BA.SwFvsSwP.up.top <- BA.SwFvsSwP.up.top %>% mutate(hitsPerc = numDEInCat * 100/numInCat)

BA.SwFvsSwP.down.top <- BA.SwFvsSwP.down.res.sig %>% slice_min(over_represented_FDR, n = 10)
BA.SwFvsSwP.down.top <- BA.SwFvsSwP.down.top %>% mutate(hitsPerc = numDEInCat * 100/numInCat)

BA.mfuzz.cluster5.top <- BA.mfuzz.cluster5.res.sig %>% slice_min(over_represented_FDR, n = 10)
BA.mfuzz.cluster5.top <- BA.mfuzz.cluster5.top %>% mutate(hitsPerc = numDEInCat * 100/numInCat)

BA.30perc.cluster5.top <- BA.30perc.cluster5.res.sig %>% slice_min(over_represented_FDR, n = 10)
BA.30perc.cluster5.top <- BA.30perc.cluster5.top %>% mutate(hitsPerc = numDEInCat * 100/numInCat)

#plot the top 10 results
plot.go <- function(go.res){
  plot.top <- ggplot(data = go.res, aes(x = hitsPerc, y = term, color = over_represented_FDR, size = numDEInCat)) +
    theme_bw() + geom_point() + expand_limits(x = 0) + labs(x = "Hits (%)", y ="GO term", color = "q value", size = "Count") +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), legend.text = element_text(size=12), 
                      legend.title = element_text(size = 13))
  
  return(plot.top)
}

BN.StFvsSwF.up.plot <- plot.go(BN.StFvsSwF.up.top)
BN.StFvsSwF.up.plot

BN.StFvsSwF.down.plot <- plot.go(BN.StFvsSwF.down.top)
BN.StFvsSwF.down.plot

BN.SwFvsSwP.up.plot <- plot.go(BN.SwFvsSwP.up.top)
BN.SwFvsSwP.up.plot

BN.SwFvsSwP.down.plot <- plot.go(BN.SwFvsSwP.down.top)
BN.SwFvsSwP.down.plot

BN.vennGenes.plot <- plot.go(BN.vennGenes.top)
BN.vennGenes.plot

BN.mfuzz.cluster5.plot <- plot.go(BN.mfuzz.cluster5.top)
BN.mfuzz.cluster5.plot

BN.mfuzz.cluster6.plot <- plot.go(BN.mfuzz.cluster6.top)
BN.mfuzz.cluster6.plot

BN.30perc.cluster7.plot <- plot.go(BN.30perc.cluster7.top)
BN.30perc.cluster7.plot

BN.sameDir.up.plot <- plot.go(BN.sameDir.up.top)
BN.sameDir.up.plot

BA.StFvsSwF.up.plot <- plot.go(BA.StFvsSwF.up.top)
BA.StFvsSwF.up.plot

BA.StFvsSwF.down.plot <- plot.go(BA.StFvsSwF.down.top)
BA.StFvsSwF.down.plot

BA.SwFvsSwP.up.plot <- plot.go(BA.SwFvsSwP.up.top)
BA.SwFvsSwP.up.plot

BA.SwFvsSwP.down.plot <- plot.go(BA.SwFvsSwP.down.top)
BA.SwFvsSwP.down.plot

BA.mfuzz.cluster5.plot <- plot.go(BA.mfuzz.cluster5.top)
BA.mfuzz.cluster5.plot

BA.30perc.cluster5.plot <- plot.go(BA.30perc.cluster5.top)
BA.30perc.cluster5.plot


