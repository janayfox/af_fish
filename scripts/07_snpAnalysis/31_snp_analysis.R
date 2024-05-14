######################################################################
### Goal: Extract results from Trinotate
### Author: Janay Fox
### R script
#####################################################################

## Set up ##
library(OutFLANK)
library(vcfR)
library(hierfstat)
library(pcadapt)
library(radiator)
library(qvalue)
library(mmod)
library(dplyr)
library(topGO)

#read pop info
pop.code.EA <- read.table("./data/snp/pop_EA.txt", header = TRUE)
pop.code.EN <- read.table("./data/snp/pop_EN.txt", header = TRUE)

colnames(pop.code.EA) <- c("INDIVIDUALS", "STRATA")
colnames(pop.code.EN) <- c("INDIVIDUALS", "STRATA")

#read in vcf files and convert format 
data.EA <- radiator::genomic_converter(
  data = "./data/snp/BA.freebayes.results.filtered.hwe.ld.norm.edited.vcf",
  strata = pop.code.EA,
  output = "hierfstat")

data.EA.evolved <- radiator::genomic_converter(
  data = "./data/snp/BA.freebayes.evolved.vcf",
  strata = pop.code.EA,
  output = "hierfstat")

data.EN <- radiator::genomic_converter(
  data = "./data/snp/BN.freebayes.results.filtered.hwe.ld.norm.edited.vcf",
  strata = pop.code.EN,
  output = "hierfstat")

data.EN.evolved <- radiator::genomic_converter(
  data = "./data/snp/BN.freebayes.evolved.vcf",
  strata = pop.code.EN,
  output = "hierfstat")

#read in lists of genesets 
EA.evolved <- read.table("./data/DEG/BA_evolved.txt")
EA.plastic <- read.table("./data/DEG/BA_plastic.tsv")

EN.evolved <- read.table("./data/DEG/BN_evolved.tsv")
EN.plastic <- read.table("./data/DEG/BN_plastic.tsv")

## Calculate FST ##
#calculate basic stats 
stat.EA <- basic.stats(data.EA$hierfstat, diploid = TRUE)
stat.EN <- basic.stats(data.EN$hierfstat, diploid = TRUE)

basic.stats(data.EA.full$hierfstat, diploid = TRUE)
basic.stats(data.EN.full$hierfstat, diploid = TRUE)

#calculate population stats 
wc(data.EA$hierfstat, diploid = TRUE)
wc(data.EA.plastic$hierfstat, diploid = TRUE)
wc(data.EA.evolved$hierfstat, diploid = TRUE)

wc(data.EA.full$hierfstat, diploid = TRUE)
wc(data.EA.plastic.full$hierfstat, diploid = TRUE)
wc(data.EA.evolved.full$hierfstat, diploid = TRUE)

wc(data.EN$hierfstat, diploid = TRUE)
wc(data.EN.plastic$hierfstat, diploid = TRUE)
wc(data.EN.evolved$hierfstat, diploid = TRUE)

wc(data.EN.full$hierfstat, diploid = TRUE)
wc(data.EN.plastic.full$hierfstat, diploid = TRUE)
wc(data.EN.evolved.full$hierfstat, diploid = TRUE)

#calculate bootstrap confidence intervals
boot.ppfst(data.EA$hierfstat, nboot=1000, diploid = TRUE)
boot.ppfst(data.EA.plastic$hierfstat, nboot=1000, diploid = TRUE)
boot.ppfst(data.EA.evolved$hierfstat, nboot=1000, diploid = TRUE)

boot.ppfst(data.EA.full$hierfstat, nboot=1000, diploid = TRUE)
boot.ppfst(data.EA.plastic.full$hierfstat, nboot=1000, diploid = TRUE)
boot.ppfst(data.EA.evolved.full$hierfstat, nboot=1000, diploid = TRUE)

boot.ppfst(data.EN$hierfstat, nboot=1000, diploid = TRUE)
boot.ppfst(data.EN.plastic$hierfstat, nboot=1000, diploid = TRUE)
boot.ppfst(data.EN.evolved$hierfstat, nboot=1000, diploid = TRUE)

boot.ppfst(data.EN.full$hierfstat, nboot=1000, diploid = TRUE)
boot.ppfst(data.EN.plastic.full$hierfstat, nboot=1000, diploid = TRUE)
boot.ppfst(data.EN.evolved.full$hierfstat, nboot=1000, diploid = TRUE)

## PCADAPT outlier analysis ##
#read in vcf files 
vcf.EA <- read.vcfR("./data/snp/BA.freebayes.results.filtered.hwe.ld.norm.vcf")
vcf.EN <- read.vcfR("./data/snp/BN.freebayes.results.filtered.hwe.ld.norm.vcf")

#extract genotypes 
geno.EA <- extract.gt(vcf.EA)
geno.EN <- extract.gt(vcf.EN)

#convert genotypes to proper format 
convert.geno <- function(geno){
  G <- geno  #if we mess up we want to be able to go back to geno
  
  G[geno %in% c("0/0")] <- 0
  G[geno  %in% c("0/1")] <- 1
  G[geno %in% c("1/1")] <- 2
  G[is.na(G)] <- 9
  tG <- t(G)
  print(dim(tG))
  return(tG)
}

G.EA <- convert.geno(geno.EA)
G.EN <- convert.geno(geno.EN)

#add on pop info 
rownames(pop.code.EA) <- pop.code.EA$INDIVIDUALS
rownames(pop.code.EN) <- pop.code.EN$INDIVIDUALS

#make sure that IDs are in the same order 
G.EA <- G.EA[order(rownames(G.EA)), , drop = FALSE]
pop.code.EA <- pop.code.EA[order(pop.code.EA$INDIVIDUALS), ]

G.EN <- G.EN[order(rownames(G.EN)), , drop = FALSE]
pop.code.EN <- pop.code.EN[order(pop.code.EN$INDIVIDUALS), ]

#check
identical(rownames(G.EA), pop.code.EA$INDIVIDUALS)
identical(rownames(G.EN), pop.code.EN$INDIVIDUALS)

#read into pcadapt 
data.pcadapt.EA.read <- read.pcadapt(t(G.EA), type = "pcadapt")
data.pcadapt.EN.read <- read.pcadapt(t(G.EN), type = "pcadapt")

#run pcadapt 
pca.EA <- pcadapt(data.pcadapt.EA.read,K=25)
pca.EN <- pcadapt(data.pcadapt.EN.read,K=25)

#determine optimum K
plot(pca.EA, option = "screeplot") # k = 5
plot(pca.EN, option = "screeplot") # k = 8

#rerun pcas 
pca.EA <- pcadapt(data.pcadapt.EA.read,K=5)
pca.EN <- pcadapt(data.pcadapt.EN.read,K=8)

#visualize to see what PC separates populations 
plot(pca.EA, option = "scores", pop = pop.code.EA$STRATA)
plot(pca.EA, option = "scores", pop = pop.code.EA$STRATA, i = 3, j = 4)
plot(pca.EA, option = "scores", pop = pop.code.EA$STRATA, i = 4, j = 5)
#PC 3 and 4 separates populations

plot(pca.EN, option = "scores", pop = pop.code.EN$STRATA)
plot(pca.EN, option = "scores", pop = pop.code.EN$STRATA, i = 3, j = 4)
plot(pca.EN, option = "scores", pop = pop.code.EN$STRATA, i = 5, j = 6)
plot(pca.EN, option = "scores", pop = pop.code.EN$STRATA, i = 7, j = 8)
#PC3 and 5

plot(pca.EA, option = "manhattan")
plot(pca.EN, option = "manhattan")

hist(pca.EA$pvalues, breaks = 50)
hist(pca.EN$pvalues, breaks = 50)

#find outliers 
qval.EA <- qvalue(pca.EA$pvalues)$qvalues
outliers.EA <- which(qval.EA<0.05)
length(outliers.EA)

qval.EN <- qvalue(pca.EN$pvalues)$qvalues
outliers.EN <- which(qval.EN<0.05)
length(outliers.EN)

#association between PC and outliers
snp.pc.EA <- get.pc(pca.EA, outliers.EA)
snp.pc.EN <- get.pc(pca.EN, outliers.EN)

#get outliers that are associated with specific PCs 
snp.pc.EA.keep <- subset(snp.pc.EA, PC == 3 | PC == 4, select = SNP)
snp.pc.EN.keep <- subset(snp.pc.EN, PC == 3 | PC == 5, select = SNP)

#find overlap 
EA.snps <- intersect(outliers.EA, snp.pc.EA.keep$SNP)
length(EA.snps)

EN.snps <- intersect(outliers.EN, snp.pc.EN.keep$SNP)
length(EN.snps)

#get names of SNPs 
snp.names.EA <- colnames(G.EA)[EA.snps]
snp.names.EN <- colnames(G.EN)[EN.snps]

#make dataframe 
snps.res.EA <- data.frame(snp_name = snp.names.EA)
snps.res.EN <- data.frame(snp_name = snp.names.EN)

#get location 
snps.res.EA$location <- sub("_[^_]*$", "", snps.res.EA$snp_name)
snps.res.EN$location <- sub("_[^_]*$", "", snps.res.EN$snp_name)

#find intersection with evolved and plastic 
intersect(snps.res.EA$location, EA.evolved$V1)
intersect(snps.res.EA$location, EA.plastic$V1)

intersect(snps.res.EN$location, EN.evolved$V1)
intersect(snps.res.EN$location, EN.plastic$V1)

#get names of all SNPs 
background.snps.EA <- colnames(G.EA)
background.snps.EA <- sub("_[^_]*$", "", background.snps.EA)

background.snps.EN <- colnames(G.EN)
background.snps.EN <- sub("_[^_]*$", "", background.snps.EN)

#save gene lists for go seq 
snp.goseq.EA <- data.frame(factor = "snp_EA",
                           name = snps.res.EA$location)

snp.goseq.EN <- data.frame(factor = "snp_EN",
                           name = snps.res.EN$location)

write.table(snp.goseq.EA, file = "./EA_SNPs.txt",
            quote = FALSE, col.names = FALSE,
            row.names = FALSE, sep = "\t")

write.table(snp.goseq.EN, file = "./EN_SNPs.txt",
            quote = FALSE, col.names = FALSE,
            row.names = FALSE, sep = "\t")

write.table(background.snps.EA, file = "./EA_background_SNPs.txt",
            quote = FALSE, col.names = FALSE,
            row.names = FALSE)

write.table(background.snps.EN, file = "./EN_background_SNPs.txt",
            quote = FALSE, col.names = FALSE,
            row.names = FALSE)

#run chi square to test for differences in proportions 
chi.data <- data.frame("sig" = c(330,1434), "not_sig" = c(118355,226268))
rownames(chi.data) <- c("BA", "BN")
chi.data <- t(chi.data)

chi.res <- chisq.test(chi.data, correct = TRUE)
chi.res


## topGO ##
#create gene universe 
gene.uni.EA <- data.frame(snp = background.snps.EA,
                          q_value = qval.EA)

gene.uni.EN <- data.frame(snp = background.snps.EN,
                          q_value = qval.EN)


#add on q values for sig snps 
qval.EA.outlier <- qval.EA[EA.snps]
snps.res.EA$q_value <- qval.EA.outlier

qval.EN.outlier <- qval.EN[EN.snps]
snps.res.EN$q_value <- qval.EN.outlier

#remove snp_name column 
go.data.EA <- snps.res.EA[, -1]
go.data.EN <- snps.res.EN[, -1]

#retain only min q value for each gene 
go.data.EA <- go.data.EA %>% group_by(location) %>%
  filter(q_value == min(q_value)) %>%
  ungroup()

go.data.EN <- go.data.EN %>% group_by(location) %>%
  filter(q_value == min(q_value)) %>%
  ungroup()

go.data.EA.uni <- gene.uni.EA %>% group_by(snp) %>%
  filter(q_value == min(q_value)) %>%
  ungroup() %>%
  distinct(snp, .keep_all = TRUE)

go.data.EN.uni <- gene.uni.EN %>% group_by(snp) %>%
  filter(q_value == min(q_value)) %>%
  ungroup() %>%
  distinct(snp, .keep_all = TRUE)

#convert to named vector 
EA.genes.vec <- setNames(go.data.EA.uni$q_value, go.data.EA.uni$snp)
EN.genes.vec <- setNames(go.data.EN.uni$q_value, go.data.EN.uni$snp)

#generate gene2go
BA.gene2GO <- readMappings(file = "./data/goSeqData/BA_go_annotations.txt")
BN.gene2GO <- readMappings(file = "./data/goSeqData/BN_go_annotations.txt")

#get gene of interest lists 
geneNames.EA <- names(BA.gene2GO)
geneList.EA <- factor(as.integer(geneNames.EA %in% go.data.EA$location)) #filter to only be sig genes
names(geneList.EA) <- geneNames.EA #add names back on

geneNames.EN <- names(BN.gene2GO)
geneList.EN <- factor(as.integer(geneNames.EN %in% go.data.EN$location)) #filter to only be sig genes
names(geneList.EN) <- geneNames.EN #add names back on

#construct data 
EA.topgo.data <- new("topGOdata",
                     description = "EN",
                     ontology = "BP",
                     allGenes = geneList.EA, 
                     nodeSize = 5, annot = annFUN.gene2GO, gene2GO = BA.gene2GO)


EN.topgo.data <- new("topGOdata",
                     description = "EN",
                     ontology = "BP",
                     allGenes = geneList.EN,
                     nodeSize = 10, annot = annFUN.gene2GO, gene2GO = BN.gene2GO)

#run test
results.BA <- runTest(EA.topgo.data, algorithm = "weight01", statistic = "fisher")
results.BN <- runTest(EN.topgo.data, algorithm = "weight01", statistic = "fisher")

#generate table of results 
allGO.BA <- usedGO(EA.topgo.data)
all.res.BA <- GenTable(EA.topgo.data, weightFisher=results.BA, orderBy='weightFisher', topNodes=length(allGO.BA))

allGO.BN <- usedGO(EN.topgo.data)
all.res.BN <- GenTable(EN.topgo.data, weightFisher=results.BN, orderBy='weightFisher', topNodes=length(allGO.BN))

#extract only significant results 
all.res.BA <- all.res.BA[all.res.BA$weightFisher <= 0.05,]
all.res.BN <- all.res.BN[all.res.BN$weightFisher <= 0.05,]

#save results 
write.table(all.res.BA, file = "./SNP_GO_enrichRes_EA.tsv", sep = "\t", row.names = FALSE)
write.table(all.res.BN, file = "./SNP_GO_enrichRes_EN.tsv", sep = "\t", row.names = FALSE)

#format for revigo 
write.table(all.res.BA[,c(1,6)], file = "./SNP_GO_enrichRes_EA_REVIGO.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(all.res.BN[,c(1,6)], file = "./SNP_GO_enrichRes_EN_REVIGO.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# ## Outflank FST outlier analysis ##
# #calculate FST 
# FST.EA <- MakeDiploidFSTMat(G.EA, locusNames = 1:ncol(G.EA), popNames = pop.code.EA$STRATA)
# FST.EN <- MakeDiploidFSTMat(G.EN, locusNames = 1:ncol(G.EN), popNames = pop.code.EN$STRATA)
# 
# #summarize
# hist(FST.EA$FST, breaks = 50)
# hist(FST.EN$FST, breaks = 50)
# 
# summary(FST.EA$FST)
# summary(FST.EN$FST)
# 
# #data checks 
# plot(FST.EA$He, FST.EA$FST)
# plot(FST.EN$He, FST.EN$FST)
# 
# plot(FST.EA$FST, FST.EA$FSTNoCorr)
# plot(FST.EN$FST, FST.EN$FSTNoCorr)
# 
# #look for outliers 
# of.EA <- OutFLANK(FST.EA, RightTrimFraction=0.01, NumberOfSamples=2, qthreshold=0.05)
# 
# OutFLANKResultsPlotter(of.EA,withOutliers=T,
#                        NoCorr=T,Hmin=0.1,binwidth=0.005,
#                        Zoom=F,RightZoomFraction=0.05,titletext=NULL)
# 
# P.EA <- pOutlierFinderChiSqNoCorr(FST.EA,Fstbar=of.EA$FSTNoCorrbar,
#                                 dfInferred=of.EA$dfInferred,qthreshold=0.05,Hmin=0.1)
# of.outliers.EA <- P.EA$OutlierFlag
# table(of.outliers.EA)
# 
# of.outliers.EA <- of.EA$results$LocusName[of.EA$results$OutlierFlag == TRUE]
# 
# of.EN <- OutFLANK(FST.EN, RightTrimFraction=0.01, NumberOfSamples=2, qthreshold=0.05)
# 
# OutFLANKResultsPlotter(of.EN,withOutliers=T,
#                        NoCorr=T,Hmin=0.1,binwidth=0.005,
#                        Zoom=F,RightZoomFraction=0.05,titletext=NULL)
# 
# P.EN <- pOutlierFinderChiSqNoCorr(FST.EN,Fstbar=of.EN$FSTNoCorrbar,
#                                   dfInferred=of.EN$dfInferred,qthreshold=0.05,Hmin=0.1)
# of.outliers.EN <- P.EN$OutlierFlag==TRUE
# table(of.outliers.EN)
# of.outliers.EN <- of.EN$results$LocusName[of.EN$results$OutlierFlag == TRUE]
# 
# ## write files for bayescan ##
# radiator::genomic_converter(
#   data = "./data/snp/BA.freebayes.results.filtered.hwe.ld.norm.edited.vcf",
#   strata = pop.code.EA,
#   output = "bayescan",
#   filename = "BA_bayescan")
# 
# radiator::genomic_converter(
#   data = "./data/snp/BN.freebayes.results.filtered.hwe.ld.norm.edited.vcf",
#   strata = pop.code.EN,
#   output = "bayescan",
#   filename = "BN_bayescan")
# 
# radiator::genomic_converter(
#   data = "./data/snp/BN.freebayes.results.filtered.ld.norm.edited.vcf",
#   strata = pop.code.EN,
#   output = "bayescan",
#   filename = "BN_bayescan_full")
# 
# radiator::genomic_converter(
#   data = "./data/snp/BA.freebayes.results.filtered.ld.norm.edited.vcf",
#   strata = pop.code.EA,
#   output = "bayescan",
#   filename = "BA_bayescan_full")
# 
# #read in bayescan results 
# bayescan.filtered.BA <- read.table("./data/snp/BA_bayescan_filtered_fst.txt")
# bayescan.filtered.BN <- read.table("./data/snp/BN_bayescan_filtered_fst.txt")
# bayescan.full.BA <- read.table("./data/snp/BA_bayescan_full_fst.txt")
# bayescan.full.BN <- read.table("./data/snp/BN_bayescan_full_fst.txt")
# 
# ## Try mmod for Nei's GST and Jost's D ##
# geneind.EA <- radiator::genomic_converter(
#   data = "./data/snp/BA.freebayes.results.filtered.hwe.ld.norm.edited.vcf",
#   strata = pop.code.EA,
#   output = "genind")
# 
# geneind.EN <- radiator::genomic_converter(
#   data = "./data/snp/BN.freebayes.results.filtered.hwe.ld.norm.edited.vcf",
#   strata = pop.code.EN,
#   output = "genind")
# 
# #calculate stats
# EA.stat <- diff_stats(geneind.EA$genind, phi_st = TRUE)
# EN.stat <- diff_stats(geneind.EN$genind, phi_st = TRUE)
# 
# Nei.EA <- pairwise_Gst_Nei(geneind.EA$genind)
# Nei.EN <- pairwise_Gst_Nei(geneind.EN$genind)
# 
# Phi.EA <- Phi_st_Meirmans(geneind.EA$genind)
# Phi.EN <- Phi_st_Meirmans(geneind.EN$genind)
# 
# D.EA <- pairwise_D(geneind.EA$genind)
# D.EN <- pairwise_D(geneind.EN$genind)
# 
# bs.EA <- chao_bootstrap(geneind.EA$genind, nreps = 100)
# summarise_bootstrap(bs.EA, Gst_Nei)     # for Nei's Gst
# 
# bs.EN <- chao_bootstrap(geneind.EN$genind, nreps = 100)
# summarise_bootstrap(bs.EN, Gst_Nei)     # for Nei's Gst
# 
# ## prepare SNP lists for PANTHER ##
# #read in snps
# BN_snps <- read.table("./data/goSeqData/EN_snps.txt")
# BA_snps <- read.table("./data/goSeqData/EA_snps.txt")
# 
# #get list of unique genes with snps 
# BN_snps <- unique(BN_snps$V2)
# BA_snps <- unique(BA_snps$V2)
# 
# #read in annotations
# BN_trinotate <- read.table("./data/trinotate/BN_trinotate.tsv", sep = "\t", fill = TRUE, quote = "", header = TRUE)
# BA_trinotate <- read.table("./data/trinotate/BA_trinotate.tsv", sep = "\t", fill = TRUE, quote = "", header = TRUE)
# 
# #reduce to only uniprot
# BN_trinotate <- BN_trinotate[c("gene_id", "sprot_Top_BLASTP_hit")]
# BA_trinotate <- BA_trinotate[c("gene_id", "sprot_Top_BLASTP_hit")]
# 
# #filter out rows with no annotation 
# BN_trinotate <- BN_trinotate %>% filter(sprot_Top_BLASTP_hit != ".")
# BA_trinotate <- BA_trinotate %>% filter(sprot_Top_BLASTP_hit != ".")
# 
# #get data for snps 
# BN_snp_anno <- BN_trinotate[BN_trinotate$gene_id %in% BN_snps, ]
# BA_snp_anno <- BA_trinotate[BA_trinotate$gene_id %in% BA_snps, ]
# 
# #edit uniprot 
# BN_snp_anno$sprot_Top_BLASTP_hit <- sapply(strsplit(as.character(BN_snp_anno$sprot_Top_BLASTP_hit), "\\^"), `[`, 1)
# BA_snp_anno$sprot_Top_BLASTP_hit <- sapply(strsplit(as.character(BA_snp_anno$sprot_Top_BLASTP_hit), "\\^"), `[`, 1)
# 
# BN_trinotate$sprot_Top_BLASTP_hit <- sapply(strsplit(as.character(BN_trinotate$sprot_Top_BLASTP_hit), "\\^"), `[`, 1)
# BA_trinotate$sprot_Top_BLASTP_hit <- sapply(strsplit(as.character(BA_trinotate$sprot_Top_BLASTP_hit), "\\^"), `[`, 1)
# 
# #save 
# write.table(unique(BN_snp_anno$sprot_Top_BLASTP_hit), "./BN_snp_PANTHER.txt",
#           quote = FALSE, col.names = FALSE,
#           row.names = FALSE, sep = "\t")
# 
# write.table(unique(BA_snp_anno$sprot_Top_BLASTP_hit), "./BA_snp_PANTHER.txt",
#             quote = FALSE, col.names = FALSE,
#             row.names = FALSE, sep = "\t")
# 
# write.table(unique(BN_trinotate$sprot_Top_BLASTP_hit), "./BN_background_PANTHER.txt",
#             quote = FALSE, col.names = FALSE,
#             row.names = FALSE, sep = "\t")
# 
# write.table(unique(BA_trinotate$sprot_Top_BLASTP_hit), "./BA_background_PANTHER.txt",
#             quote = FALSE, col.names = FALSE,
#             row.names = FALSE, sep = "\t")
# 
# #get go annotation lists 
# BN_go <- read.table("./data/goSeqData/BN_go_annotations.txt")
# BA_go <- read.table("./data/goSeqData/BA_go_annotations.txt")
# 
# #get data for snps 
# BN_snp_anno_go <- BN_go[BN_go$V1 %in% BN_snps, ]
# BA_snp_anno_go <- BA_go[BA_go$V1 %in% BA_snps, ]
# 
# #write file 
# file_conn <- file("go_snp_BN.txt", "w")
# all_go_ids <- unlist(strsplit(BN_snp_anno_go$V2, ","))
# writeLines(all_go_ids, file_conn)
# close(file_conn)
