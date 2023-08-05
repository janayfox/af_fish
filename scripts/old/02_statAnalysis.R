######################################################################
### Goal: PCA and cluster analysis on DEG and total gene expression
### Author: Janay Fox
### R script
#####################################################################

## Set up ##
library(pheatmap)
library(ggplot2)
library(factoextra)
library(data.table)
library(RColorBrewer)
library(tibble)
library(stringr)
library(dendextend)
library(tidyverse)
library(cluster)
library(ggupset)
library(VennDiagram)
library(RColorBrewer)
library(GGally)
library(Mfuzz)
library(dplyr)
library(ggpubr)

#read in data 
degBN <- read.table("./data/DEG/BN/salmon/diffExpr.P0.01_C2.matrix")
degBA <- read.table("./data/DEG/BA/salmon/diffExpr.P0.01_C2.matrix")

allExprBN <- read.table("./data/totalExpr/BN_bf_new_sal.gene.TMM.EXPR.matrix")
allExprBA <- read.table("./data/totalExpr/BA_bf_new_sal.gene.TMM.EXPR.matrix")

BN.StF.vs.StP.DEG <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.stream_field_vs_stream_pond.edgeR.DE_results.P0.01_C2.DE.subset")
BN.StF.vs.SwF.DEG <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.stream_field_vs_swamp_field.edgeR.DE_results.P0.01_C2.DE.subset")
BN.StP.vs.SwP.DEG <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.stream_pond_vs_swamp_pond.edgeR.DE_results.P0.01_C2.DE.subset")
BN.SwF.vs.SwP.DEG <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.swamp_field_vs_swamp_pond.edgeR.DE_results.P0.01_C2.DE.subset")

BA.StF.vs.StP.DEG <- read.table("./data/DEG/BA/salmon/BA_bf_new_sal.gene.counts.matrix.stream_field_vs_stream_pond.edgeR.DE_results.P0.01_C2.DE.subset")
BA.StF.vs.SwF.DEG <- read.table("./data/DEG/BA/salmon/BA_bf_new_sal.gene.counts.matrix.stream_field_vs_swamp_field.edgeR.DE_results.P0.01_C2.DE.subset")
BA.StP.vs.SwP.DEG <- read.table("./data/DEG/BA/salmon/BA_bf_new_sal.gene.counts.matrix.stream_pond_vs_swamp_pond.edgeR.DE_results.P0.01_C2.DE.subset")
BA.SwF.vs.SwP.DEG <- read.table("./data/DEG/BA/salmon/BA_bf_new_sal.gene.counts.matrix.swamp_field_vs_swamp_pond.edgeR.DE_results.P0.01_C2.DE.subset")

BN.StF.vs.StP.result <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.stream_field_vs_stream_pond.edgeR.DE_results")
BN.StF.vs.SwF.result <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.stream_field_vs_swamp_field.edgeR.DE_results")
BN.StP.vs.SwP.result <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.stream_pond_vs_swamp_pond.edgeR.DE_results")
BN.SwF.vs.SwP.result <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.swamp_field_vs_swamp_pond.edgeR.DE_results")

BA.StF.vs.StP.result <- read.table("./data/DEG/BA/salmon/BA_bf_new_sal.gene.counts.matrix.stream_field_vs_stream_pond.edgeR.DE_results")
BA.StF.vs.SwF.result <- read.table("./data/DEG/BA/salmon/BA_bf_new_sal.gene.counts.matrix.stream_field_vs_swamp_field.edgeR.DE_results")
BA.StP.vs.SwP.result <- read.table("./data/DEG/BA/salmon/BA_bf_new_sal.gene.counts.matrix.stream_pond_vs_swamp_pond.edgeR.DE_results")
BA.SwF.vs.SwP.result <- read.table("./data/DEG/BA/salmon/BA_bf_new_sal.gene.counts.matrix.swamp_field_vs_swamp_pond.edgeR.DE_results")

BN.pond.data <- read.csv("./data/pondData/pond_data_BN.csv")
rownames(BN.pond.data) <- BN.pond.data[,1] #convert first column to row names 
BN.pond.data[,1] <- NULL

BA.pond.data <- read.csv("./data/pondData/pond_data_BA.csv")
rownames(BA.pond.data) <- BA.pond.data[,1] 
BA.pond.data[,1] <- NULL

##  Make UpSet plots and venn diagrams ## 
#get dataframe in proper format 
BN.StF.vs.StP.DEG.list <- data.frame(gene = rownames(BN.StF.vs.StP.DEG), comparison = "H-I vs H-A")
BN.StF.vs.SwF.DEG.list <- data.frame(gene = rownames(BN.StF.vs.SwF.DEG), comparison = "H-I vs L-I")
BN.StP.vs.SwP.DEG.list <- data.frame(gene = rownames(BN.StP.vs.SwP.DEG), comparison = "H-A vs L-A")
BN.SwF.vs.SwP.DEG.list <- data.frame(gene = rownames(BN.SwF.vs.SwP.DEG), comparison = "L-I vs L-A")

BN.DEG.upset <- rbind(BN.StF.vs.StP.DEG.list, BN.StF.vs.SwF.DEG.list, BN.StP.vs.SwP.DEG.list, BN.SwF.vs.SwP.DEG.list)
BN.DEG.upset <- BN.DEG.upset %>% group_by(gene) %>% summarize(comparisons = list(comparison))

BA.StF.vs.StP.DEG.list <- data.frame(gene = rownames(BA.StF.vs.StP.DEG), comparison = "H-I vs H-A")
BA.StF.vs.SwF.DEG.list <- data.frame(gene = rownames(BA.StF.vs.SwF.DEG), comparison = "H-I vs L-I")
BA.StP.vs.SwP.DEG.list <- data.frame(gene = rownames(BA.StP.vs.SwP.DEG), comparison = "H-A vs L-A")
BA.SwF.vs.SwP.DEG.list <- data.frame(gene = rownames(BA.SwF.vs.SwP.DEG), comparison = "L-I vs L-A")

BA.DEG.upset <- rbind(BA.StF.vs.StP.DEG.list, BA.StF.vs.SwF.DEG.list, BA.StP.vs.SwP.DEG.list, BA.SwF.vs.SwP.DEG.list)
BA.DEG.upset <- BA.DEG.upset %>% group_by(gene) %>% summarize(comparisons = list(comparison))

#Plot upset 
tiff("BN_upset.tiff", units="in", width = 7, height = 8, res = 600)
BN.upset <- ggplot(data = BN.DEG.upset, aes(x = comparisons)) + geom_bar(fill = "#27187E") + scale_x_upset() + theme_bw() + ylab("Intersection Size") + 
                  xlab("") + geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5, color = "black", size = 3) + 
                    theme_combmatrix(combmatrix.panel.point.color.fill = "#2a9d8f",
                                     combmatrix.panel.line.color = "#2a9d8f",
                                     combmatrix.label.extra_spacing	= 10,
                                     combmatrix.label.text = element_text(size = 12, color = "black")) + 
                    theme(axis.title = element_text(size = 15, color = "black"), axis.text.y = element_text(size = 11, color = "black")) 

BN.upset
dev.off()

tiff("BA_upset.tiff", units="in", width = 7, height = 8, res = 600)
BA.upset <- ggplot(data = BA.DEG.upset, aes(x = comparisons)) + geom_bar(fill = "#27187E") + scale_x_upset() + theme_bw() + ylab("Intersection Size") + 
  xlab("") + geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5, color = "black", size = 3) + 
  theme_combmatrix(combmatrix.panel.point.color.fill = "#2a9d8f",
                   combmatrix.panel.line.color = "#2a9d8f",
                   combmatrix.label.extra_spacing	= 10,
                   combmatrix.label.text = element_text(size = 12, color = "black")) + 
  theme(axis.title = element_text(size = 15, color = "black"), axis.text.y = element_text(size = 11, color = "black")) 

BA.upset
dev.off()

#Plot venn diagram 
#function to display venn diagram 
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

display_venn(BN.DEG.list.upset , lwd = 2, col = c("#e76f51", "#0081a7", "#FFBC42", "#a7c957"),
             fill = c(alpha("#e76f51", 0.8), alpha("#0081a7", 0.8), alpha("#FFBC42",0.8), alpha("#a7c957",0.8)),
             fontface = "bold", cex = 1.5, cat.cex = 1.5)

display_venn(BA.DEG.list.upset, lwd = 2, col = c("#e76f51", "#0081a7", "#FFBC42", "#a7c957"),
             fill = c(alpha("#e76f51", 0.8), alpha("#0081a7", 0.8), alpha("#FFBC42",0.8), alpha("#a7c957",0.8)),
             fontface = "bold", cex = 1.5, cat.cex = 1.5)

## Run and Plot PCA on DEG##
#make list of DEG from all four important comparisons 
BN.DEG.list <- c(BN.StF.vs.StP.DEG.list$gene, BN.StF.vs.SwF.DEG.list$gene, BN.SwF.vs.SwP.DEG.list$gene, BN.StP.vs.SwP.DEG.list$gene)
BN.DEG.list <- unique(BN.DEG.list) #remove duplicates

BA.DEG.list <- c(BA.StF.vs.StP.DEG.list$gene, BA.StF.vs.SwF.DEG.list$gene, BA.SwF.vs.SwP.DEG.list$gene, BA.StP.vs.SwP.DEG.list$gene)
BA.DEG.list <- unique(BA.DEG.list) #remove duplicates

#extract DEG
extractDEG <- function(DEG.df, DEG.list){
  return(DEG.df[DEG.list,])
}

BN.DEG <- extractDEG(allExprBN, BN.DEG.list)
BA.DEG <- extractDEG(allExprBA, BA.DEG.list)

#log2 transform data 
log.trans <- function(df){
  return(log((df + 1), 2))
}

BN.DEG.log <- log.trans(BN.DEG)
BA.DEG.log <- log.trans(BA.DEG)

#median center data 
med.center <- function(df){
  rowMed <- apply(df, 1, median)
  return(df-rowMed)
}

BN.DEG.log.center <- med.center(BN.DEG.log)
BA.DEG.log.center <- med.center(BA.DEG.log)

#transpose data 
transpose.names <- function(data){
  t.data <- t(data)
  rownames(t.data) <- colnames(data)
  colnames(t.data) <- rownames(data)
  return(as.data.frame(t.data))
}

t.BN.DEG.log.center <- transpose.names(BN.DEG.log.center)
t.BA.DEG.log.center <- transpose.names(BA.DEG.log.center)

#run PCA
pca.BN.DEG <- prcomp(t.BN.DEG.log.center)
pca.BA.DEG <- prcomp(t.BA.DEG.log.center)

#view results 
summary(pca.BN.DEG)
summary(pca.BA.DEG)

#extract first 4 PCs
t.BN.DEG.log.center <- cbind(t.BN.DEG.log.center, pca.BN.DEG$x[,1:4])
t.BA.DEG.log.center <- cbind(t.BA.DEG.log.center, pca.BA.DEG$x[,1:4])

#add on site and collection data 
add.site.collection.BN <- function(data.BN) {
  data.BN.new <- tibble::rownames_to_column(data.BN, "group")
  data.BN.new$group <- str_sub(data.BN.new$group, end = -13)
  
  data.BN.new$site <- data.BN.new$group
  data.BN.new$site <- str_sub(data.BN.new$site, end = -6)
  data.BN.new$site <- gsub("_", "", data.BN.new$site)
  
  data.BN.new$collection <- data.BN.new$group
  data.BN.new$collection <- str_sub(data.BN.new$collection, start = 7)
  data.BN.new$collection <- gsub("_", "", data.BN.new$collection)
  
  rownames(data.BN.new) <- rownames(data.BN)
  return(data.BN.new)
}

add.site.collection.BA <- function(data.BA) {
  data.BA.new <- tibble::rownames_to_column(data.BA, "group")
  data.BA.new$group <- str_sub(data.BA.new$group, end = -10)
  data.BA.new$group <- gsub("d_", "d", data.BA.new$group)
  
  data.BA.new$site <- data.BA.new$group
  data.BA.new$site <- str_sub(data.BA.new$site, end = -6)
  data.BA.new$site <- gsub("_", "", data.BA.new$site)
  
  data.BA.new$collection <- data.BA.new$group
  data.BA.new$collection <- str_sub(data.BA.new$collection, start = 7)
  data.BA.new$collection <- gsub("_", "", data.BA.new$collection)
  
  rownames(data.BA.new) <- rownames(data.BA)
  
  return(data.BA.new)
}

t.BN.DEG.log.center <- add.site.collection.BN(t.BN.DEG.log.center)
t.BA.DEG.log.center <- add.site.collection.BA(t.BA.DEG.log.center)

#plot PCA
plot.pca <- function(pca.data, pca.x, pca.y, xlabel, ylabel){
  ggplot(data = pca.data, aes(x = pca.x, y = pca.y, color = group, shape = group)) + theme_bw() + geom_point() +
    scale_color_manual(name = "Source and Collection",
                      values = c("#27187E", "#82DCDA", "#ad2e24", "#ff9b54"), 
                       labels = c("H-I", "H-A", 
                                  "L-I", "L-A")) + 
    scale_shape_manual(name = "Source and Collection",
                       labels = c("H-I", "H-A", 
                                  "L-I", "L-A"),
                       values = c(15, 17, 15, 17)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
    stat_ellipse(aes(group = group)) + xlab(xlabel) + ylab(ylabel) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), legend.text = element_text(size=12), 
          legend.title = element_text(size = 13))
}

tiff("BN_pca12_DEG.tiff", units="in", width = 8, height = 6, res = 600)
BN.DEG.pca12.plot <- plot.pca(t.BN.DEG.log.center, t.BN.DEG.log.center$PC1, 
                                  t.BN.DEG.log.center$PC2, "PC1 (43.4%)", "PC2 (13.2%)")
BN.DEG.pca12.plot
dev.off()

tiff("BN_pca23_DEG.tiff", units="in", width = 8, height = 6, res = 600)
BN.DEG.pca23.plot <- plot.pca(t.BN.DEG.log.center, t.BN.DEG.log.center$PC2, 
                                  t.BN.DEG.log.center$PC3, "PC2 (13.2%)", "PC3 (9.3%)")
BN.DEG.pca23.plot
dev.off()

tiff("BA_pca12_DEG.tiff", units="in", width = 8, height = 6, res = 600)
BA.DEG.pca12.plot <- plot.pca(t.BA.DEG.log.center, t.BA.DEG.log.center$PC1, 
                                  t.BA.DEG.log.center$PC2, "PC1 (61.5%)", "PC2 (13.1%)")
BA.DEG.pca12.plot
dev.off()

tiff("BA_pca23_DEG.tiff", units="in", width = 8, height = 6, res = 600)
BA.DEG.pca23.plot <- plot.pca(t.BA.DEG.log.center, t.BA.DEG.log.center$PC2, 
                                  t.BA.DEG.log.center$PC3, "PC2 (13.1%)", "PC3 (7.6%)")
BA.DEG.pca23.plot
dev.off()

#plot pond information 
#add on pond information
t.BN.DEG.log.center <- merge(t.BN.DEG.log.center, BN.pond.data, by = "row.names", all = TRUE)
rownames(t.BN.DEG.log.center) <- t.BN.DEG.log.center[,1]
t.BN.DEG.log.center[,1] <- NULL

t.BA.DEG.log.center <- merge(t.BA.DEG.log.center, BA.pond.data, by = "row.names", all = TRUE)
rownames(t.BA.DEG.log.center) <- t.BA.DEG.log.center[,1]
t.BA.DEG.log.center[,1] <- NULL

#plot PCA with pond information 
plot.pond.pca <- function(pca.data, pca.x, pca.y, xlabel, ylabel){
  ggplot(data = pca.data, aes(x = pca.x, y = pca.y, color = pond, shape = group)) + theme_bw() + geom_point() +
    scale_color_manual(values = c("#ad2e24", "#27187E", "#82DCDA")) + 
    scale_shape_manual(values = c(16, 15, 17, 4), labels = c("H-I", "H-A", 
                                                             "L-I", "L-A")) +
    labs(color = "Pond", shape = "Group") + geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") + stat_ellipse(aes(group = pond)) + xlab(xlabel) + ylab(ylabel) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), legend.text = element_text(size=11))
}

tiff("BN_pca12_DEG_pond.tiff", units="in", width = 8, height = 6, res = 600)
BN.DEG.pca12.pond.plot <- plot.pond.pca(t.BN.DEG.log.center, t.BN.DEG.log.center$PC1, t.BN.DEG.log.center$PC2, "PC1 (43.4%)", "PC2 (13.2%)")
BN.DEG.pca12.pond.plot
dev.off()

tiff("BN_pca23_DEG_pond.tiff", units="in", width = 8, height = 6, res = 600)
BN.DEG.pca23.pond.plot <- plot.pond.pca(t.BN.DEG.log.center, t.BN.DEG.log.center$PC2, t.BN.DEG.log.center$PC3, "PC2 (13.2%)", "PC3 (9.3%)")
BN.DEG.pca23.pond.plot
dev.off()

tiff("BA_pca12_DEG_pond.tiff", units="in", width = 8, height = 6, res = 600)
BA.DEG.pca12.pond.plot <- plot.pond.pca(t.BA.DEG.log.center, t.BA.DEG.log.center$PC1, t.BA.DEG.log.center$PC2, "PC1 (61.5%)", "PC2 (13.1%)")
BA.DEG.pca12.pond.plot
dev.off()

tiff("BA_pca23_DEG_pond.tiff", units="in", width = 8, height = 6, res = 600)
BA.DEG.pca23.pond.plot <- plot.pond.pca(t.BA.DEG.log.center, t.BA.DEG.log.center$PC2, t.BA.DEG.log.center$PC3, "PC2 (13.1%)", "PC3 (7.6%)")
BA.DEG.pca23.pond.plot
dev.off()

##Plot PCA on total gene expression##
#transform data 
allExprBN.log <- log.trans(allExprBN)
allExprBA.log <- log.trans(allExprBA)

allExprBN.log.center <- med.center(allExprBN.log)
allExprBA.log.center <- med.center(allExprBA.log)

#transpose data 
t.allExprBN.log.center <- transpose.names(allExprBN.log.center)
t.allExprBA.log.center <- transpose.names(allExprBA.log.center)

#run PCA
pca.allExpr.BN <- prcomp(t.allExprBN.log.center)
pca.allExpr.BA <- prcomp(t.allExprBA.log.center)

#add on site and collection data 
t.allExprBN.log.center <- add.site.collection.BN(t.allExprBN.log.center)
t.allExprBA.log.center <- add.site.collection.BA(t.allExprBA.log.center)

#view results 
summary(pca.allExpr.BN)
summary(pca.allExpr.BA)

#extract first 4 PCs
t.allExprBN.log.center <- cbind(t.allExprBN.log.center, pca.allExpr.BN$x[,1:4])
t.allExprBA.log.center <- cbind(t.allExprBA.log.center, pca.allExpr.BA$x[,1:4])

#plot PCA
tiff("BN_pca12_allExpr.tiff", units="in", width = 8, height = 6, res = 600)
allExprBN.pca12.plot <- plot.pca(t.allExprBN.log.center, t.allExprBN.log.center$PC1, t.allExprBN.log.center$PC2, "PC1 (9%)", "PC2 (5.8%)")
allExprBN.pca12.plot
dev.off()

tiff("BN_pca23_allExpr.tiff", units="in", width = 8, height = 6, res = 600)
allExprBN.pca23.plot <- plot.pca(t.allExprBN.log.center, t.allExprBN.log.center$PC2, t.allExprBN.log.center$PC3, "PC2 (5.8%)", "PC3 (4.8%)")
allExprBN.pca23.plot
dev.off()

tiff("BA_pca12_allExpr.tiff", units="in", width = 8, height = 6, res = 600)
allExprBA.pca12.plot <- plot.pca(t.allExprBA.log.center, t.allExprBA.log.center$PC1, t.allExprBA.log.center$PC2, "PC1 (14.4%)", "PC2 (7.3%)")
allExprBA.pca12.plot
dev.off()

tiff("BA_pca23_allExpr.tiff", units="in", width = 8, height = 6, res = 600)
allExprBA.pca23.plot <- plot.pca(t.allExprBA.log.center, t.allExprBA.log.center$PC2, t.allExprBA.log.center$PC3, "PC2 (7.3%)", "PC3 (6.2%)")
allExprBA.pca23.plot
dev.off()

#plot pond information 
#add on pond information
t.allExprBN.log.center <- merge(t.allExprBN.log.center, BN.pond.data, by = "row.names", all = TRUE)
rownames(t.allExprBN.log.center) <- t.allExprBN.log.center[,1]
t.allExprBN.log.center[,1] <- NULL

t.allExprBA.log.center <- merge(t.allExprBA.log.center, BA.pond.data, by = "row.names", all = TRUE)
rownames(t.allExprBA.log.center) <- t.allExprBA.log.center[,1]
t.allExprBA.log.center[,1] <- NULL

#plot pca with pond information
tiff("BN_pca12_allExpr_pond.tiff", units="in", width = 8, height = 6, res = 600)
allExprBN.pca12.pond.plot <- plot.pond.pca(t.allExprBN.log.center, t.allExprBN.log.center$PC1, t.allExprBN.log.center$PC2, "PC1 (9%)", "PC2 (5.8%)")
allExprBN.pca12.pond.plot
dev.off()

tiff("BN_pca23_allExpr_pond.tiff", units="in", width = 8, height = 6, res = 600)
allExprBN.pca23.pond.plot <- plot.pond.pca(t.allExprBN.log.center, t.allExprBN.log.center$PC2, t.allExprBN.log.center$PC3, "PC2 (5.8%)", "PC3 (4.8%)")
allExprBN.pca23.pond.plot
dev.off()

tiff("BA_pca12_allExpr_pond.tiff", units="in", width = 8, height = 6, res = 600)
allExprBA.pca12.pond.plot <- plot.pond.pca(t.allExprBA.log.center, t.allExprBA.log.center$PC1, t.allExprBA.log.center$PC2, "PC1 (14.4%)", "PC2 (7.3%)")
allExprBA.pca12.pond.plot
dev.off()

tiff("BA_pca23_allExpr_pond.tiff", units="in", width = 8, height = 6, res = 600)
allExprBA.pca23.pond.plot <- plot.pond.pca(t.allExprBA.log.center, t.allExprBA.log.center$PC2, t.allExprBA.log.center$PC3, "PC2 (7.3%)", "PC3 (6.2%)")
allExprBA.pca23.pond.plot
dev.off()

## Run cluster analysis and plot heatmap on DEG##
#do cluster analysis on samplesand genes
BN.DEG.sample.hclust <- hclust(dist(t(BN.DEG.log.center), method = "euclidean"), method = "ward.D2")
BA.DEG.sample.hclust <- hclust(dist(t(BA.DEG.log.center), method = "euclidean"), method = "ward.D2")

BN.DEG.gene.hclust <- hclust(dist(BN.DEG.log.center, method = "euclidean"), method = "ward.D2")
BA.DEG.gene.hclust <- hclust(dist(BA.DEG.log.center, method = "euclidean"), method = "ward.D2")

#make dataframe for sample annotations 
sample.info.BN <- data.frame(Group = t.BN.DEG.log.center$group)
rownames(sample.info.BN) <- row.names(t.BN.DEG.log.center)
sample.info.BN$Group <- gsub("stream_field", "H-I", sample.info.BN$Group)
sample.info.BN$Group <- gsub("stream_pond", "H-A", sample.info.BN$Group)
sample.info.BN$Group <- gsub("swamp_field", "L-I", sample.info.BN$Group)
sample.info.BN$Group <- gsub("swamp_pond", "L-A", sample.info.BN$Group)

sample.info.BA <- data.frame(Group = t.BA.DEG.log.center$group)
rownames(sample.info.BA) <- row.names(t.BA.DEG.log.center)
sample.info.BA$Group <- gsub("stream_field", "H-I", sample.info.BA$Group)
sample.info.BA$Group <- gsub("stream_pond", "H-A", sample.info.BA$Group)
sample.info.BA$Group <- gsub("swamp_field", "L-I", sample.info.BA$Group)
sample.info.BA$Group <- gsub("swamp_pond", "L-A", sample.info.BA$Group)

#find gene clusters by cutting tree at % height 
find.gene.clust <- function(hclust.dat, perc.val){
  gene.clust <- cutree(tree = as.dendrogram(hclust.dat), h = (perc.val*max(hclust.dat$height))) 
  print(max(gene.clust))
  return(gene.clust)
}

BN.gene.20perc.clust.dat <- find.gene.clust(BN.DEG.gene.hclust, 0.2)
BA.gene.20perc.clust.dat <- find.gene.clust(BA.DEG.gene.hclust, 0.2)

BN.gene.30perc.clust.dat <- find.gene.clust(BN.DEG.gene.hclust, 0.3)
BA.gene.30perc.clust.dat <- find.gene.clust(BA.DEG.gene.hclust, 0.3)

#make cluster info into dataframe for annotations
gene.clust.BN.20perc <- data.frame(Cluster = BN.gene.20perc.clust.dat)
gene.clust.BN.20perc$Cluster <- as.factor(gene.clust.BN.20perc$Cluster)
gene.clust.BN.20perc$Cluster <- gsub("\\<1\\>", "Cluster1", gene.clust.BN.20perc$Cluster)
gene.clust.BN.20perc$Cluster <- gsub("\\<2\\>", "Cluster2", gene.clust.BN.20perc$Cluster)
gene.clust.BN.20perc$Cluster <- gsub("\\<3\\>", "Cluster3", gene.clust.BN.20perc$Cluster)
gene.clust.BN.20perc$Cluster <- gsub("\\<4\\>", "Cluster4", gene.clust.BN.20perc$Cluster)
gene.clust.BN.20perc$Cluster <- gsub("\\<5\\>", "Cluster5", gene.clust.BN.20perc$Cluster)
gene.clust.BN.20perc$Cluster <- gsub("\\<6\\>", "Cluster6", gene.clust.BN.20perc$Cluster)
gene.clust.BN.20perc$Cluster <- gsub("\\<7\\>", "Cluster7", gene.clust.BN.20perc$Cluster)
gene.clust.BN.20perc$Cluster <- gsub("\\<8\\>", "Cluster8", gene.clust.BN.20perc$Cluster)
gene.clust.BN.20perc$Cluster <- gsub("\\<9\\>", "Cluster9", gene.clust.BN.20perc$Cluster)
gene.clust.BN.20perc$Cluster <- gsub("\\<10\\>", "Cluster10", gene.clust.BN.20perc$Cluster)
gene.clust.BN.20perc$Cluster <- gsub("\\<11\\>", "Cluster11", gene.clust.BN.20perc$Cluster)
gene.clust.BN.20perc$Cluster <- gsub("\\<12\\>", "Cluster12", gene.clust.BN.20perc$Cluster)
gene.clust.BN.20perc$Cluster <- gsub("\\<13\\>", "Cluster13", gene.clust.BN.20perc$Cluster)
gene.clust.BN.20perc$Cluster <- gsub("\\<14\\>", "Cluster14", gene.clust.BN.20perc$Cluster)
gene.clust.BN.20perc$Cluster <- gsub("\\<15\\>", "Cluster15", gene.clust.BN.20perc$Cluster)

gene.clust.BA.20perc <- data.frame(Cluster = BA.gene.20perc.clust.dat)
gene.clust.BA.20perc$Cluster <- as.factor(gene.clust.BA.20perc$Cluster)
gene.clust.BA.20perc$Cluster <- gsub("\\<1\\>", "Cluster1", gene.clust.BA.20perc$Cluster)
gene.clust.BA.20perc$Cluster <- gsub("\\<2\\>", "Cluster2", gene.clust.BA.20perc$Cluster)
gene.clust.BA.20perc$Cluster <- gsub("\\<3\\>", "Cluster3", gene.clust.BA.20perc$Cluster)
gene.clust.BA.20perc$Cluster <- gsub("\\<4\\>", "Cluster4", gene.clust.BA.20perc$Cluster)
gene.clust.BA.20perc$Cluster <- gsub("\\<5\\>", "Cluster5", gene.clust.BA.20perc$Cluster)
gene.clust.BA.20perc$Cluster <- gsub("\\<6\\>", "Cluster6", gene.clust.BA.20perc$Cluster)
gene.clust.BA.20perc$Cluster <- gsub("\\<7\\>", "Cluster7", gene.clust.BA.20perc$Cluster)
gene.clust.BA.20perc$Cluster <- gsub("\\<8\\>", "Cluster8", gene.clust.BA.20perc$Cluster)
gene.clust.BA.20perc$Cluster <- gsub("\\<9\\>", "Cluster9", gene.clust.BA.20perc$Cluster)
gene.clust.BA.20perc$Cluster <- gsub("\\<10\\>", "Cluster10", gene.clust.BA.20perc$Cluster)
gene.clust.BA.20perc$Cluster <- gsub("\\<11\\>", "Cluster11", gene.clust.BA.20perc$Cluster)

gene.clust.BN.30perc <- data.frame(Cluster = BN.gene.30perc.clust.dat)
gene.clust.BN.30perc$Cluster <- as.factor(gene.clust.BN.30perc$Cluster)
gene.clust.BN.30perc$Cluster <- gsub("\\<1\\>", "Cluster1", gene.clust.BN.30perc$Cluster)
gene.clust.BN.30perc$Cluster <- gsub("\\<2\\>", "Cluster2", gene.clust.BN.30perc$Cluster)
gene.clust.BN.30perc$Cluster <- gsub("\\<3\\>", "Cluster3", gene.clust.BN.30perc$Cluster)
gene.clust.BN.30perc$Cluster <- gsub("\\<4\\>", "Cluster4", gene.clust.BN.30perc$Cluster)
gene.clust.BN.30perc$Cluster <- gsub("\\<5\\>", "Cluster5", gene.clust.BN.30perc$Cluster)
gene.clust.BN.30perc$Cluster <- gsub("\\<6\\>", "Cluster6", gene.clust.BN.30perc$Cluster)
gene.clust.BN.30perc$Cluster <- gsub("\\<7\\>", "Cluster7", gene.clust.BN.30perc$Cluster)
gene.clust.BN.30perc$Cluster <- gsub("\\<8\\>", "Cluster8", gene.clust.BN.30perc$Cluster)
gene.clust.BN.30perc$Cluster <- gsub("\\<9\\>", "Cluster9", gene.clust.BN.30perc$Cluster)

gene.clust.BA.30perc <- data.frame(Cluster = BA.gene.30perc.clust.dat)
gene.clust.BA.30perc$Cluster <- as.factor(gene.clust.BA.30perc$Cluster)
gene.clust.BA.30perc$Cluster <- gsub("\\<1\\>", "Cluster1", gene.clust.BA.30perc$Cluster)
gene.clust.BA.30perc$Cluster <- gsub("\\<2\\>", "Cluster2", gene.clust.BA.30perc$Cluster)
gene.clust.BA.30perc$Cluster <- gsub("\\<3\\>", "Cluster3", gene.clust.BA.30perc$Cluster)
gene.clust.BA.30perc$Cluster <- gsub("\\<4\\>", "Cluster4", gene.clust.BA.30perc$Cluster)
gene.clust.BA.30perc$Cluster <- gsub("\\<5\\>", "Cluster5", gene.clust.BA.30perc$Cluster)
gene.clust.BA.30perc$Cluster <- gsub("\\<6\\>", "Cluster6", gene.clust.BA.30perc$Cluster)
gene.clust.BA.30perc$Cluster <- gsub("\\<7\\>", "Cluster7", gene.clust.BA.30perc$Cluster)
gene.clust.BA.30perc$Cluster <- gsub("\\<8\\>", "Cluster8", gene.clust.BA.30perc$Cluster)
gene.clust.BA.30perc$Cluster <- gsub("\\<9\\>", "Cluster9", gene.clust.BA.30perc$Cluster)

#set colors for annotations
BN.clust.20perc.anno.colors <- list(
  Group = c("H-DO,Imm." = "#2a9d8f", "H-DO,Acc." = "#e9c46a",
            "L-DO,Imm." = "#f4a261", "L-DO,Acc." = "#e76f51"),
  Cluster = c(Cluster1 = "#780000", Cluster2 = "#ea5545", Cluster3 = "#f46a9b", Cluster4 = "#fb5607", Cluster5 =  "#ef9b20",
              Cluster6 = "#ffd500", Cluster7 = "#b5e48c", Cluster8 =  "#87bc45", Cluster9 = "#006400", 
              Cluster10 = "#80ced7", Cluster11 = "#27aeef", Cluster12 = "#00509d", Cluster13 = "#b33dc6", 
              Cluster14 = "#9d4edd", Cluster15 = "#5a189a"))

BA.clust.20perc.anno.colors <- list(
  Group = c("H-DO,Imm." = "#2a9d8f", "H-DO,Acc." = "#e9c46a",
            "L-DO,Imm." = "#f4a261", "L-DO,Acc." = "#e76f51"),
  Cluster = c(Cluster1 = "#780000", Cluster2 = "#ea5545", Cluster3 = "#f46a9b", Cluster4 = "#fb5607", Cluster5 =  "#ef9b20",
              Cluster6 = "#ffd500", Cluster7 = "#b5e48c", Cluster8 =  "#87bc45", Cluster9 = "#006400", 
              Cluster10 = "#80ced7", Cluster11 = "#27aeef"))

clust.30perc.anno.colors <- list(
  Group = c("H-DO,Imm." = "#2a9d8f", "H-DO,Acc." = "#e9c46a",
            "L-DO,Imm." = "#f4a261", "L-DO,Acc." = "#e76f51"),
  Cluster = c(Cluster1 = "#ea5545", Cluster2 = "#f46a9b", Cluster3 = "#fb5607", Cluster4 = "#ffd500", Cluster5 =  "#b5e48c",
              Cluster6 = "#87bc45", Cluster7 = "#27aeef", Cluster8 =  "#b33dc6", Cluster9 = "#5a189a"))

#generate color pallette 
cor.color <- colorRampPalette(c("#FFFF00", "black", "#0000FF"))(60) 

#calculate breaks 
breaks = c(seq(-6, 0, length.out = 30),seq(0.01, 6, length.out = 30))

set.breaks <- function(expr.data){
  new.breaks <- breaks
  new.breaks[length(breaks)] <- max(max(expr.data),max(breaks))
  new.breaks[1] <- min(min(expr.data),min(breaks))
  return(new.breaks)
}

BN.DEG.breaks <- set.breaks(BN.DEG.log.center)
BA.DEG.breaks <- set.breaks(BA.DEG.log.center)

#plot heatmaps
make.heatmap <- function(DEG, sample.anno, gene.anno, sample.hclust, gene.hclust, breaks.samp, anno.colors){
  pheatmap(DEG, cluster_rows = gene.hclust, cluster_cols = sample.hclust, show_colnames = FALSE, show_rownames = FALSE,
           color = cor.color, breaks = breaks.samp, annotation_row = gene.anno,
           annotation_col = sample.anno, annotation_colors = anno.colors,
           annotation_names_col = FALSE, annotation_names_row = FALSE, legend_breaks = c(-6,-3,0,3,6), border_color = NA)
}

tiff("BN_heatmap_DEG_20perc.tiff", units="in", width = 7, height = 8, res = 600)
BN.DEG.heatmap.20perc <- make.heatmap(BN.DEG.log.center, sample.info.BN, gene.clust.BN.20perc, BN.DEG.sample.hclust, BN.DEG.gene.hclust, BN.DEG.breaks, BN.clust.20perc.anno.colors)
BN.DEG.heatmap.20perc
dev.off()

tiff("BA_heatmap_DEG_20perc.tiff", units="in", width = 7, height = 8, res = 600)
BA.DEG.heatmap.20perc <- make.heatmap(BA.DEG.log.center, sample.info.BA, gene.clust.BA.20perc, BA.DEG.sample.hclust, BA.DEG.gene.hclust, BA.DEG.breaks, BA.clust.20perc.anno.colors)
BA.DEG.heatmap.20perc
dev.off()

tiff("BN_heatmap_DEG_30perc.tiff", units="in", width = 7, height = 8, res = 600)
BN.DEG.heatmap.30perc <- make.heatmap(BN.DEG.log.center, sample.info.BN, gene.clust.BN.30perc, BN.DEG.sample.hclust, BN.DEG.gene.hclust, BN.DEG.breaks, clust.30perc.anno.colors)
BN.DEG.heatmap.30perc
dev.off()

tiff("BA_heatmap_DEG_30perc.tiff", units="in", width = 7, height = 8, res = 600)
BA.DEG.heatmap.TMM.30perc <- make.heatmap(BA.DEG.log.center, sample.info.BA, gene.clust.BA.30perc, BA.DEG.sample.hclust, BA.DEG.gene.hclust, BA.DEG.breaks, clust.30perc.anno.colors)
BA.DEG.heatmap.TMM.30perc
dev.off()

#make list of genes in clusters of interest 
BN.30perc.cluster7 <- rownames(subset(gene.clust.BN.30perc, Cluster == "Cluster7"))
BA.30perc.cluster5 <- rownames(subset(gene.clust.BA.30perc, Cluster == "Cluster5"))

#try cutting tree at k means 5-15 (or until there are repeated patterns)
find.gene.clust.k <- function(hclust.dat, k.val, cluster.name){
  gene.clust <- cutree(tree = as.dendrogram(hclust.dat), k = k.val) 
  gene.clust <- enframe(gene.clust, name = "gene", value = cluster.name)
  return(gene.clust)
}

BN.gene.clust.dat.k5 <- find.gene.clust.k(BN.DEG.gene.hclust, 5, "cluster.k5")
BN.gene.clust.dat.k6 <- find.gene.clust.k(BN.DEG.gene.hclust, 6, "cluster.k6")
BN.gene.clust.dat.k7 <- find.gene.clust.k(BN.DEG.gene.hclust, 7, "cluster.k7")
BN.gene.clust.dat.k8 <- find.gene.clust.k(BN.DEG.gene.hclust, 8, "cluster.k8")
BN.gene.clust.dat.k9 <- find.gene.clust.k(BN.DEG.gene.hclust, 9, "cluster.k9")
BN.gene.clust.dat.k10 <- find.gene.clust.k(BN.DEG.gene.hclust, 10, "cluster.k10")
BN.gene.clust.dat.k11 <- find.gene.clust.k(BN.DEG.gene.hclust, 11, "cluster.k11")
BN.gene.clust.dat.k12 <- find.gene.clust.k(BN.DEG.gene.hclust, 12, "cluster.k12")
BN.gene.clust.dat.k13 <- find.gene.clust.k(BN.DEG.gene.hclust, 13, "cluster.k13")
BN.gene.clust.dat.k14 <- find.gene.clust.k(BN.DEG.gene.hclust, 14, "cluster.k14")
BN.gene.clust.dat.k15 <- find.gene.clust.k(BN.DEG.gene.hclust, 15, "cluster.k15")

BA.gene.clust.dat.k5 <- find.gene.clust.k(BA.DEG.gene.hclust, 5, "cluster.k5")
BA.gene.clust.dat.k6 <- find.gene.clust.k(BA.DEG.gene.hclust, 6, "cluster.k6")
BA.gene.clust.dat.k7 <- find.gene.clust.k(BA.DEG.gene.hclust, 7, "cluster.k7")
BA.gene.clust.dat.k8 <- find.gene.clust.k(BA.DEG.gene.hclust, 8, "cluster.k8")
BA.gene.clust.dat.k9 <- find.gene.clust.k(BA.DEG.gene.hclust, 9, "cluster.k9")
BA.gene.clust.dat.k10 <- find.gene.clust.k(BA.DEG.gene.hclust, 10, "cluster.k10")
BA.gene.clust.dat.k11 <- find.gene.clust.k(BA.DEG.gene.hclust, 11, "cluster.k11")
BA.gene.clust.dat.k12 <- find.gene.clust.k(BA.DEG.gene.hclust, 12, "cluster.k12")
BA.gene.clust.dat.k13 <- find.gene.clust.k(BA.DEG.gene.hclust, 13, "cluster.k13")
BA.gene.clust.dat.k14 <- find.gene.clust.k(BA.DEG.gene.hclust, 14, "cluster.k14")
BA.gene.clust.dat.k15 <- find.gene.clust.k(BA.DEG.gene.hclust, 15, "cluster.k15")

#make version of gene expression with row names as column and in right order for plotting 
BN.DEG.log.center.ordered <- BN.DEG.log.center[,c(6,7,14,15,16,17,18,27,28,1,2,12,13,22,24,26,30,4,8,11,20,21,31,32,3,5,9,10,19,23,25,29)]
BA.DEG.log.center.ordered <- BA.DEG.log.center[,c(7,11,15,17,18,23,25,32,34,5,6,20,21,27,30,31,1,2,3,4,13,14,22,28,29,36,37,8,9,10,12,16,19,24,26,33,35)]

#change names of columns 
colnames(BN.DEG.log.center.ordered) <- c("HDO-I.1", "HDO-I.2", "HDO-I.3", "HDO-I.4",
                                         "HDO-I.5", "HDO-I.6", "HDO-I.7", "HDO-I.8",
                                         "HDO-I.9","HDO-A.1", "HDO-A.2", "HDO-A.3",
                                         "HDO-A.4", "HDO-A.5", "HDO-A.6", "HDO-A.7",
                                         "HDO-A.8", "LDO-A.1", "LDO-A.2", "LDO-A.3",
                                         "LDO-A.4", "LDO-A.5", "LDO-A.6", "LDO-A.7",
                                         "LDO-I.1", "LDO-I.2", "LDO-I.3", "LDO-I.4",
                                         "LDO-I.5", "LDO-I.6", "LDO-I.7", "LDO-I.8")

colnames(BA.DEG.log.center.ordered) <- c("HDO-I.1", "HDO-I.2", "HDO-I.3", "HDO-I.4",
                                        "HDO-I.5", "HDO-I.6", "HDO-I.7", "HDO-I.8",
                                        "HDO-I.9","HDO-A.1", "HDO-A.2", "HDO-A.3",
                                        "HDO-A.4", "HDO-A.5", "HDO-A.6", "HDO-A.7",
                                        "LDO-A.1", "LDO-A.2", "LDO-A.3", "LDO-A.4",
                                        "LDO-A.5", "LDO-A.6", "LDO-A.7", "LDO-A.8",
                                        "LDO-A.9", "LDO-A.10", "LDO-A.11", "LDO-I.1",
                                        "LDO-I.2", "LDO-I.3", "LDO-I.4", "LDO-I.5",
                                        "LDO-I.6", "LDO-I.7", "LDO-I.8", "LDO-I.9",
                                        "LDO-I.10")

#change rownames to column 
BN.DEG.log.center.ordered <- rownames_to_column(BN.DEG.log.center.ordered, "gene")
BA.DEG.log.center.ordered <- rownames_to_column(BA.DEG.log.center.ordered, "gene")

#make list of all kmeans clusters 
BN.gene.clust.list <- list(BN.DEG.log.center.ordered, BN.gene.clust.dat.k5, BN.gene.clust.dat.k6, BN.gene.clust.dat.k7, BN.gene.clust.dat.k8,
                       BN.gene.clust.dat.k9, BN.gene.clust.dat.k10, BN.gene.clust.dat.k11, BN.gene.clust.dat.k12,
                       BN.gene.clust.dat.k13, BN.gene.clust.dat.k14, BN.gene.clust.dat.k15)

BA.gene.clust.list <- list(BA.DEG.log.center.ordered, BA.gene.clust.dat.k5, BA.gene.clust.dat.k6, BA.gene.clust.dat.k7, BA.gene.clust.dat.k8,
                       BA.gene.clust.dat.k9, BA.gene.clust.dat.k10, BA.gene.clust.dat.k11, BA.gene.clust.dat.k12,
                       BA.gene.clust.dat.k13, BA.gene.clust.dat.k14, BA.gene.clust.dat.k15)

#merge gene expression data with all kmeans cluster data
BN.DEG.log.center.kclusters <- BN.gene.clust.list %>% reduce(full_join, by = "gene")
BA.DEG.log.center.kclusters <- BA.gene.clust.list %>% reduce(full_join, by = "gene")

#convert to long format 
BN.DEG.log.center.k.long <- pivot_longer(BN.DEG.log.center.kclusters, 2:33, names_to = "sample", values_to = "relativeExpression")
BA.DEG.log.center.k.long <- pivot_longer(BA.DEG.log.center.kclusters, 2:37, names_to = "sample", values_to = "relativeExpression")

#visualize expression per cluster 
plot.kClusters <- function(k.data, cluster.k){
  ggplot(k.data, aes_(substitute(sample), substitute(relativeExpression))) + geom_line(aes_(group = substitute(gene)), alpha = 0.3) + facet_grid(rows = vars({{cluster.k}})) +
    theme_bw() + xlab("Sample") + ylab("Relative Gene Expression") + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11),
                                                                           axis.text.x = element_text(angle = 60, vjust = 0.5))   
} 

tiff("BN_5K_cluster.tiff", units="in", width = 12, height = 8, res = 600)
BN.5k <- plot.kClusters(BN.DEG.log.center.k.long, BN.DEG.log.center.k.long$cluster.k5)
BN.5k
dev.off()

tiff("BN_6K_cluster.tiff", units="in", width = 12, height = 8, res = 600)
BN.6k <- plot.kClusters(BN.DEG.log.center.k.long, cluster.k6)
BN.6k
dev.off()

tiff("BN_7K_cluster.tiff", units="in", width = 12, height = 8, res = 600)
BN.7k <- plot.kClusters(BN.DEG.log.center.k.long, cluster.k7)
BN.7k
dev.off()

tiff("BN_8K_cluster.tiff", units="in", width = 12, height = 10, res = 600)
BN.8k <- plot.kClusters(BN.DEG.log.center.k.long, cluster.k8)
BN.8k
dev.off()

tiff("BN_9K_cluster.tiff", units="in", width = 12, height = 12, res = 600)
BN.9k <- plot.kClusters(BN.DEG.log.center.k.long, cluster.k9)
BN.9k
dev.off()

tiff("BN_10K_cluster.tiff", units="in", width = 12, height = 12, res = 600)
BN.10k <- plot.kClusters(BN.DEG.log.center.k.long, cluster.k10)
BN.10k
dev.off()

tiff("BN_11K_cluster.tiff", units="in", width = 12, height = 12, res = 600)
BN.11k <- plot.kClusters(BN.DEG.log.center.k.long, cluster.k11)
BN.11k
dev.off()

tiff("BN_12K_cluster.tiff", units="in", width = 12, height = 14, res = 600)
BN.12k <- plot.kClusters(BN.DEG.log.center.k.long, cluster.k12)
BN.12k
dev.off()

tiff("BN_13K_cluster.tiff", units="in", width = 12, height = 14, res = 600)
BN.13k <- plot.kClusters(BN.DEG.log.center.k.long, cluster.k13)
BN.13k
dev.off()

tiff("BN_14K_cluster.tiff", units="in", width = 12, height = 16, res = 600)
BN.14k <- plot.kClusters(BN.DEG.log.center.k.long, cluster.k14)
BN.14k
dev.off()

tiff("BN_15K_cluster.tiff", units="in", width = 12, height = 16, res = 600)
BN.15k <- plot.kClusters(BN.DEG.log.center.k.long, cluster.k15)
BN.15k
dev.off()

tiff("BA_5K_cluster.tiff", units="in", width = 12, height = 8, res = 600)
BA.5k <- plot.kClusters(BA.DEG.log.center.k.long, BA.DEG.log.center.k.long$cluster.k5)
BA.5k
dev.off()

tiff("BA_6K_cluster.tiff", units="in", width = 12, height = 8, res = 600)
BA.6k <- plot.kClusters(BA.DEG.log.center.k.long, cluster.k6)
BA.6k
dev.off()

tiff("BA_7K_cluster.tiff", units="in", width = 12, height = 8, res = 600)
BA.7k <- plot.kClusters(BA.DEG.log.center.k.long, cluster.k7)
BA.7k
dev.off()

tiff("BA_8K_cluster.tiff", units="in", width = 12, height = 10, res = 600)
BA.8k <- plot.kClusters(BA.DEG.log.center.k.long, cluster.k8)
BA.8k
dev.off()

tiff("BA_9K_cluster.tiff", units="in", width = 12, height = 12, res = 600)
BA.9k <- plot.kClusters(BA.DEG.log.center.k.long, cluster.k9)
BA.9k
dev.off()

tiff("BA_10K_cluster.tiff", units="in", width = 12, height = 12, res = 600)
BA.10k <- plot.kClusters(BA.DEG.log.center.k.long, cluster.k10)
BA.10k
dev.off()

tiff("BA_11K_cluster.tiff", units="in", width = 12, height = 12, res = 600)
BA.11k <- plot.kClusters(BA.DEG.log.center.k.long, cluster.k11)
BA.11k
dev.off()

tiff("BA_12K_cluster.tiff", units="in", width = 12, height = 14, res = 600)
BA.12k <- plot.kClusters(BA.DEG.log.center.k.long, cluster.k12)
BA.12k
dev.off()

tiff("BA_13K_cluster.tiff", units="in", width = 12, height = 14, res = 600)
BA.13k <- plot.kClusters(BA.DEG.log.center.k.long, cluster.k13)
BA.13k
dev.off()

tiff("BA_14K_cluster.tiff", units="in", width = 12, height = 16, res = 600)
BA.14k <- plot.kClusters(BA.DEG.log.center.k.long, cluster.k14)
BA.14k
dev.off()

tiff("BA_15K_cluster.tiff", units="in", width = 12, height = 16, res = 600)
BA.15k <- plot.kClusters(BA.DEG.log.center.k.long, cluster.k15)
BA.15k
dev.off()

#plot heatmap with cluster results 
#change genes to rownames
BN.gene.clust.dat.k8 <- as.data.frame(BN.gene.clust.dat.k8)
rownames(BN.gene.clust.dat.k8) <- BN.gene.clust.dat.k8[,1]
BN.gene.clust.dat.k8[,1] <- NULL

BA.gene.clust.dat.k8 <- as.data.frame(BA.gene.clust.dat.k8)
rownames(BA.gene.clust.dat.k8) <- BA.gene.clust.dat.k8[,1]
BA.gene.clust.dat.k8[,1] <- NULL

#change name of cluster column 
colnames(BN.gene.clust.dat.k8) <- "Cluster"
colnames(BA.gene.clust.dat.k8) <- "Cluster"

#change names of clusters
BN.gene.clust.dat.k8$Cluster <- as.factor(BN.gene.clust.dat.k8$Cluster)
BN.gene.clust.dat.k8$Cluster <- gsub("\\<1\\>", "Cluster1", BN.gene.clust.dat.k8$Cluster)
BN.gene.clust.dat.k8$Cluster <- gsub("\\<2\\>", "Cluster2", BN.gene.clust.dat.k8$Cluster)
BN.gene.clust.dat.k8$Cluster <- gsub("\\<3\\>", "Cluster3", BN.gene.clust.dat.k8$Cluster)
BN.gene.clust.dat.k8$Cluster <- gsub("\\<4\\>", "Cluster4", BN.gene.clust.dat.k8$Cluster)
BN.gene.clust.dat.k8$Cluster <- gsub("\\<5\\>", "Cluster5", BN.gene.clust.dat.k8$Cluster)
BN.gene.clust.dat.k8$Cluster <- gsub("\\<6\\>", "Cluster6", BN.gene.clust.dat.k8$Cluster)
BN.gene.clust.dat.k8$Cluster <- gsub("\\<7\\>", "Cluster7", BN.gene.clust.dat.k8$Cluster)
BN.gene.clust.dat.k8$Cluster <- gsub("\\<8\\>", "Cluster8", BN.gene.clust.dat.k8$Cluster)

BA.gene.clust.dat.k8$Cluster <- as.factor(BA.gene.clust.dat.k8$Cluster)
BA.gene.clust.dat.k8$Cluster <- gsub("\\<1\\>", "Cluster1", BA.gene.clust.dat.k8$Cluster)
BA.gene.clust.dat.k8$Cluster <- gsub("\\<2\\>", "Cluster2", BA.gene.clust.dat.k8$Cluster)
BA.gene.clust.dat.k8$Cluster <- gsub("\\<3\\>", "Cluster3", BA.gene.clust.dat.k8$Cluster)
BA.gene.clust.dat.k8$Cluster <- gsub("\\<4\\>", "Cluster4", BA.gene.clust.dat.k8$Cluster)
BA.gene.clust.dat.k8$Cluster <- gsub("\\<5\\>", "Cluster5", BA.gene.clust.dat.k8$Cluster)
BA.gene.clust.dat.k8$Cluster <- gsub("\\<6\\>", "Cluster6", BA.gene.clust.dat.k8$Cluster)
BA.gene.clust.dat.k8$Cluster <- gsub("\\<7\\>", "Cluster7", BA.gene.clust.dat.k8$Cluster)
BA.gene.clust.dat.k8$Cluster <- gsub("\\<8\\>", "Cluster8", BA.gene.clust.dat.k8$Cluster)

#set annotation colors
k.clust.anno.colors <- list(
  Group = c("H-DO,Imm." = "#2a9d8f", "H-DO,Acc." = "#e9c46a",
            "L-DO,Imm." = "#f4a261", "L-DO,Acc." = "#e76f51"),
  Cluster = c(Cluster1 = "#ea5545", Cluster2 = "#f46a9b", Cluster3 = "#ffd500", Cluster4 = "#87bc45", Cluster5 =  "#27aeef",
              Cluster6 = "#00509d", Cluster7 = "#b33dc6", Cluster8 =  "#5a189a"))

#plot
tiff("BN_heatmap_DEG_kclust.tiff", units="in", width = 7, height = 8, res = 600)
BN.DEG.heatmap.kclust <- make.heatmap(BN.DEG.log.center, sample.info.BN, BN.gene.clust.dat.k8, BN.DEG.sample.hclust, BN.DEG.gene.hclust, BN.DEG.breaks, k.clust.anno.colors)
BN.DEG.heatmap.kclust
dev.off()

tiff("BA_heatmap_DEG_kclust.tiff", units="in", width = 7, height = 8, res = 600)
BA.DEG.heatmap.kclust <- make.heatmap(BA.DEG.log.center, sample.info.BA, BA.gene.clust.dat.k8, BA.DEG.sample.hclust, BA.DEG.gene.hclust, BA.DEG.breaks, k.clust.anno.colors)
BA.DEG.heatmap.kclust
dev.off()

#make list of genes in clusters of interest 
BN.k.cluster2 <- rownames(subset(BN.gene.clust.dat.k8, Cluster == "Cluster2"))
BA.k.cluster3 <- rownames(subset(BA.gene.clust.dat.k8, Cluster == "Cluster3"))

## Plot Bar Graphs comparing numbers of DEGS ## 
#construct dataframe 
comparison <- c("H-DO,Imm. vs H-DO,Acc.", "H-DO,Imm. vs L-DO,Imm.", "L-DO,Imm. vs L-DO,Acc.", "H-DO,Acc. vs L-DO,Acc.",
                    "H-DO,Imm. vs H-DO,Acc.", "H-DO,Imm. vs L-DO,Imm.", "L-DO,Imm. vs L-DO,Acc.", "H-DO,Acc. vs L-DO,Acc.")

species <- c("EN", "EN", "EN", "EN", "EA", "EA", "EA", "EA")

numDEG <- c(nrow(BN.StF.vs.StP.DEG), nrow(BN.StF.vs.SwF.DEG), nrow(BN.SwF.vs.SwP.DEG), nrow(BN.StP.vs.SwP.DEG), 
                nrow(BA.StF.vs.StP.DEG), nrow(BA.StF.vs.SwF.DEG), nrow(BA.SwF.vs.SwP.DEG), nrow(BA.StP.vs.SwP.DEG))

compareDEG <- data.frame(numDEG, species, comparison)

#order factor in order I want in barplot
compareDEG$comparison <- factor(compareDEG$comparison, levels = c("H-DO,Imm. vs H-DO,Acc.", "L-DO,Imm. vs L-DO,Acc.", 
                                                                                  "H-DO,Imm. vs L-DO,Imm.", "H-DO,Acc. vs L-DO,Acc."))
#plot
plot.compareDEG<- ggplot(data = compareDEG, aes(x = comparison, y = numDEG, fill = species)) + geom_bar(stat = "identity", position = position_dodge()) + 
  geom_text(aes(label=numDEG), vjust=1.6, color="white",position = position_dodge(0.9), size=3.5)+
  theme_bw() + scale_fill_manual(values = c("#2a9d8f", "#27187E"), labels = c("E. apleurogramma", "E. neumayeri")) + 
  xlab("Comparison Type") + ylab ("Number of DEGs") + labs(fill = "Species") +
  theme(legend.position = "bottom", axis.title = element_text(size = 14), axis.text = element_text(size = 11),
        legend.text = element_text(size=11)) + scale_x_discrete(labels = function(x) str_wrap(x, width = 5))

tiff("barplot.tiff", units="in", width = 6, height = 7, res = 600)
plot.compareDEG
dev.off()

#plot barplots separately 
compareDEG.BN <- compareDEG[compareDEG$species == "EN",]
compareDEG.BA <- compareDEG[compareDEG$species == "EA",]

plot.compareDEG.BN <- ggplot(data = compareDEG.BN, aes(x = comparison, y = numDEG)) + geom_bar(stat = "identity", position = position_dodge(), fill = "#2a9d8f") + 
                              geom_text(aes(label=numDEG), vjust=1.6, color="white",position = position_dodge(0.9), size=3.5) +
                              theme_bw()  + xlab("Comparison Type") + ylab ("Number of DEGs") +
                              theme(legend.position = "bottom", axis.title = element_text(size = 14), axis.text = element_text(size = 11)) + 
                              scale_x_discrete(labels = function(x) str_wrap(x, width = 5))

plot.compareDEG.BN

plot.compareDEG.BA <- ggplot(data = compareDEG.BA, aes(x = comparison, y = numDEG)) + geom_bar(stat = "identity", position = position_dodge(), fill = "#2a9d8f") + 
  geom_text(aes(label=numDEG), vjust=1.6, color="white",position = position_dodge(0.9), size=3.5) +
  theme_bw()  + xlab("Comparison Type") + ylab ("Number of DEGs") +
  theme(legend.position = "bottom", axis.title = element_text(size = 14), axis.text = element_text(size = 11)) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 5))

plot.compareDEG.BA


## Partition genes into clusters using soft clustering ## 
#try mfuzz method 
#reorganize data 
BN.DEG.log <- BN.DEG.log[,c(6,7,14,15,16,17,18,27,28,1,2,12,13,22,24,26,30,4,8,11,20,21,31,32,3,5,9,10,19,23,25,29)]
BA.DEG.log <- BA.DEG.log[,c(7,11,15,17,18,23,25,32,34,5,6,20,21,27,30,31,1,2,3,4,13,14,22,28,29,36,37,8,9,10,12,16,19,24,26,33,35)]

#convert to eset
BN.DEG.log.matrix <- as.matrix(BN.DEG.log)
BN.DEG.log.eset <- new("ExpressionSet", exprs = BN.DEG.log.matrix)

BA.DEG.log.matrix <- as.matrix(BA.DEG.log)
BA.DEG.log.eset <- new("ExpressionSet", exprs = BA.DEG.log.matrix)

#standardize
BN.DEG.log.eset.std <- standardise(BN.DEG.log.eset)
BA.DEG.log.eset.std <- standardise(BA.DEG.log.eset)

#estimate fuzzifier parameter
m.BN <- mestimate(BN.DEG.log.eset.std)
m.BN
m.BA <- mestimate(BA.DEG.log.eset.std)
m.BA

#estimate optimum number of clusters
cselect.BN <- cselection(BN.DEG.log.eset.std, m = m.BN, crange = seq(2,50,2), repeats = 5)
dmin.BN <- Dmin(BN.DEG.log.eset.std, m = m.BN, crange = seq(2,50,2), repeats = 5)

cselect.BA <- cselection(BA.DEG.log.eset.std, m = m.BA, crange = seq(2,50,2), repeats = 5, visu = TRUE)
dmin.BA <- Dmin(BA.DEG.log.eset.std, m = m.BA, crange = seq(2,50,2), repeats = 5)

#run clustering to visualize clusters
mfuzz.BN.2 <- mfuzz(BN.DEG.log.eset.std, m = m.BN, c = 2)
mfuzz.BN.3 <- mfuzz(BN.DEG.log.eset.std, m = m.BN, c = 3)
mfuzz.BN.4 <- mfuzz(BN.DEG.log.eset.std, m = m.BN, c = 4)
mfuzz.BN.5 <- mfuzz(BN.DEG.log.eset.std, m = m.BN, c = 5)
mfuzz.BN.6 <- mfuzz(BN.DEG.log.eset.std, m = m.BN, c = 6)
mfuzz.BN.7 <- mfuzz(BN.DEG.log.eset.std, m = m.BN, c = 7)
mfuzz.BN.8 <- mfuzz(BN.DEG.log.eset.std, m = m.BN, c = 8)
mfuzz.BN.9 <- mfuzz(BN.DEG.log.eset.std, m = m.BN, c = 9)
mfuzz.BN.10 <- mfuzz(BN.DEG.log.eset.std, m = m.BN, c = 10)
mfuzz.BN.11 <- mfuzz(BN.DEG.log.eset.std, m = m.BN, c = 11)
mfuzz.BN.12 <- mfuzz(BN.DEG.log.eset.std, m = m.BN, c = 12)
mfuzz.BN.13 <- mfuzz(BN.DEG.log.eset.std, m = m.BN, c = 13)
mfuzz.BN.14 <- mfuzz(BN.DEG.log.eset.std, m = m.BN, c = 14)
mfuzz.BN.15 <- mfuzz(BN.DEG.log.eset.std, m = m.BN, c = 15)
mfuzz.BN.16 <- mfuzz(BN.DEG.log.eset.std, m = m.BN, c = 16)
mfuzz.BN.17 <- mfuzz(BN.DEG.log.eset.std, m = m.BN, c = 17)
mfuzz.BN.18 <- mfuzz(BN.DEG.log.eset.std, m = m.BN, c = 18)
mfuzz.BN.19 <- mfuzz(BN.DEG.log.eset.std, m = m.BN, c = 19)
mfuzz.BN.20 <- mfuzz(BN.DEG.log.eset.std, m = m.BN, c = 20)

mfuzz.BA.2 <- mfuzz(BA.DEG.log.eset.std, m = m.BA, c = 2)
mfuzz.BA.3 <- mfuzz(BA.DEG.log.eset.std, m = m.BA, c = 3)
mfuzz.BA.4 <- mfuzz(BA.DEG.log.eset.std, m = m.BA, c = 4)
mfuzz.BA.5 <- mfuzz(BA.DEG.log.eset.std, m = m.BA, c = 5)
mfuzz.BA.6 <- mfuzz(BA.DEG.log.eset.std, m = m.BA, c = 6)
mfuzz.BA.7 <- mfuzz(BA.DEG.log.eset.std, m = m.BA, c = 7)
mfuzz.BA.8 <- mfuzz(BA.DEG.log.eset.std, m = m.BA, c = 8)
mfuzz.BA.9 <- mfuzz(BA.DEG.log.eset.std, m = m.BA, c = 9)
mfuzz.BA.10 <- mfuzz(BA.DEG.log.eset.std, m = m.BA, c = 10)
mfuzz.BA.11 <- mfuzz(BA.DEG.log.eset.std, m = m.BA, c = 11)
mfuzz.BA.12 <- mfuzz(BA.DEG.log.eset.std, m = m.BA, c = 12)
mfuzz.BA.13 <- mfuzz(BA.DEG.log.eset.std, m = m.BA, c = 13)
mfuzz.BA.14 <- mfuzz(BA.DEG.log.eset.std, m = m.BA, c = 14)
mfuzz.BA.15 <- mfuzz(BA.DEG.log.eset.std, m = m.BA, c = 15)
mfuzz.BA.16 <- mfuzz(BA.DEG.log.eset.std, m = m.BA, c = 16)
mfuzz.BA.17 <- mfuzz(BA.DEG.log.eset.std, m = m.BA, c = 17)
mfuzz.BA.18 <- mfuzz(BA.DEG.log.eset.std, m = m.BA, c = 18)
mfuzz.BA.19 <- mfuzz(BA.DEG.log.eset.std, m = m.BA, c = 19)
mfuzz.BA.20 <- mfuzz(BA.DEG.log.eset.std, m = m.BA, c = 20)

#plot clusters
mfuzz.plot.BN <- function(filename, mfuzz.cl, r, c){
  tiff(filename, units="in", width = 14, height = 14, res = 600)
  par(mar=c(5,6,4,1) + 0.1)
  mfuzz.plot2(BN.DEG.log.eset.std, cl = mfuzz.cl, min.mem = 0.7, 
              time.labels = c("HDO.I", "HDO.I", "HDO.I", "HDO.I",
                              "HDO.I", "HDO.I", "HDO.I", "HDO.I",
                              "HDO.I","HDO.A", "HDO.A", "HDO.A",
                              "HDO.A", "HDO.A", "HDO.A", "HDO.A",
                              "HDO.A", "LDO.A", "LDO.A", "LDO.A",
                              "LDO.A", "LDO.A", "LDO.A", "LDO.A",
                              "LDO.I", "LDO.I", "LDO.I", "LDO.I",
                              "LDO.I", "LDO.I", "LDO.I", "LDO.I"),
              xlab = "Source & Collection", ylab = "Relative Expression", 
              mfrow = c(r,c), x11 = FALSE, cex.lab = 2.5, cex.main = 3, cex.axis = 1.75)
  dev.off()
  }

mfuzz.plot.BA <- function(filename, mfuzz.cl, r, c){
  tiff(filename, units="in", width = 12, height = 12, res = 600)
  par(mar=c(5,6,4,1) + 0.1)
  mfuzz.plot2(BA.DEG.log.eset.std, cl = mfuzz.cl, min.mem = 0.7, 
              time.labels = c("HDO.I", "HDO.I", "HDO.I", "HDO.I",
                            "HDO.I", "HDO.I", "HDO.I", "HDO.I",
                            "HDO.I","HDO.A", "HDO.A", "HDO.A",
                            "HDO.A", "HDO.A", "HDO.A", "HDO.A",
                            "LDO.A", "LDO.A", "LDO.A", "LDO.A",
                            "LDO.A", "LDO.A", "LDO.A", "LDO.A",
                            "LDO.A", "LDO.A", "LDO.A", "LDO.I",
                            "LDO.I", "LDO.I", "LDO.I", "LDO.I",
                            "LDO.I", "LDO.I", "LDO.I", "LDO.I",
                            "LDO.I"),
              xlab = "Source & Collection", ylab = "Relative Expression", 
              mfrow = c(r,c), x11 = FALSE, cex.lab = 2.5, cex.main = 3, cex.axis = 1.75)
  dev.off()
}

mfuzz.plot.BN("BN2_mfuzz_plot.tiff", mfuzz.BN.2, 2, 1)
mfuzz.plot.BN("BN3_mfuzz_plot.tif", mfuzz.BN.3, 2, 2)
mfuzz.plot.BN("BN4_mfuzz_plot.tiff", mfuzz.BN.4, 2, 2)
mfuzz.plot.BN("BN5_mfuzz_plot.tiff", mfuzz.BN.5, 3, 2)
mfuzz.plot.BN("BN6_mfuzz_plot.tiff", mfuzz.BN.6, 3, 2)
mfuzz.plot.BN("BN7_mfuzz_plot.tiff", mfuzz.BN.7, 3, 3)
mfuzz.plot.BN("BN8_mfuzz_plot.tiff", mfuzz.BN.8, 3, 3)
mfuzz.plot.BN("BN9_mfuzz_plot.tiff", mfuzz.BN.9, 3, 3)
mfuzz.plot.BN("BN10_mfuzz_plot.tiff", mfuzz.BN.10, 4, 3)
mfuzz.plot.BN("BN11_mfuzz_plot.tiff", mfuzz.BN.11, 4, 3)
mfuzz.plot.BN("BN12_mfuzz_plot.tiff", mfuzz.BN.12, 4, 3)
mfuzz.plot.BN("BN13_mfuzz_plot.tiff", mfuzz.BN.13, 4, 4)
mfuzz.plot.BN("BN14_mfuzz_plot.tiff", mfuzz.BN.14, 4, 4)
mfuzz.plot.BN("BN15_mfuzz_plot.tiff", mfuzz.BN.15, 4, 4)
mfuzz.plot.BN("BN16_mfuzz_plot.tiff", mfuzz.BN.16, 4, 4)
mfuzz.plot.BN("BN17_mfuzz_plot.tiff", mfuzz.BN.17, 5, 4)
mfuzz.plot.BN("BN18_mfuzz_plot.tiff", mfuzz.BN.18, 5, 4)
mfuzz.plot.BN("BN19_mfuzz_plot.tiff", mfuzz.BN.19, 5, 4)
mfuzz.plot.BN("BN20_mfuzz_plot.tiff", mfuzz.BN.20, 5, 4)

mfuzz.plot.BA("BA2_mfuzz_plot.tiff", mfuzz.BA.2, 2, 1)
mfuzz.plot.BA("BA3_mfuzz_plot.tiff", mfuzz.BA.3, 2, 2)
mfuzz.plot.BA("BA4_mfuzz_plot.tiff", mfuzz.BA.4, 2, 2)
mfuzz.plot.BA("BA5_mfuzz_plot.tiff", mfuzz.BA.5, 3, 2)
mfuzz.plot.BA("BA6_mfuzz_plot.tiff", mfuzz.BA.6, 3, 2)
mfuzz.plot.BA("BA7_mfuzz_plot.tiff", mfuzz.BA.7, 3, 3)
mfuzz.plot.BA("BA8_mfuzz_plot.tiff", mfuzz.BA.8, 3, 3)
mfuzz.plot.BA("BA9_mfuzz_plot.tiff", mfuzz.BA.9, 3, 3)
mfuzz.plot.BA("BA10_mfuzz_plot.tiff", mfuzz.BA.10, 4, 3)
mfuzz.plot.BA("BA11_mfuzz_plot.tiff", mfuzz.BA.11, 4, 3)
mfuzz.plot.BA("BA12_mfuzz_plot.tiff", mfuzz.BA.12, 4, 3)
mfuzz.plot.BA("BA13_mfuzz_plot.tiff", mfuzz.BA.13, 4, 4)
mfuzz.plot.BA("BA14_mfuzz_plot.tiff", mfuzz.BA.14, 4, 4)
mfuzz.plot.BA("BA15_mfuzz_plot.tiff", mfuzz.BA.15, 4, 4)
mfuzz.plot.BA("BA16_mfuzz_plot.tiff", mfuzz.BA.16, 4, 4)
mfuzz.plot.BA("BA17_mfuzz_plot.tiff", mfuzz.BA.17, 5, 4)
mfuzz.plot.BA("BA18_mfuzz_plot.tiff", mfuzz.BA.18, 5, 4)
mfuzz.plot.BA("BA19_mfuzz_plot.tiff", mfuzz.BA.19, 5, 4)
mfuzz.plot.BA("BA20_mfuzz_plot.tiff", mfuzz.BA.20, 5, 4)

#extract cluster info  
BN.gene.clusters <- as.data.frame(mfuzz.BN.9$cluster)
colnames(BN.gene.clusters) <- "cluster"
BN.membership <- as.data.frame(mfuzz.BN.9$membership)
  
BA.gene.clusters <- as.data.frame(mfuzz.BA.7$cluster)
colnames(BA.gene.clusters) <- "cluster"
BA.membership <- as.data.frame(mfuzz.BA.7$membership)

#merge datasets
BN.gene.cluster.info <- merge(BN.gene.clusters, BN.membership, by = "row.names", all = TRUE)
rownames(BN.gene.cluster.info) <- BN.gene.cluster.info[,1]
BN.gene.cluster.info[,1] <- NULL

BA.gene.cluster.info <- merge(BA.gene.clusters, BA.membership, by = "row.names", all = TRUE)
rownames(BA.gene.cluster.info) <- BA.gene.cluster.info[,1]
BA.gene.cluster.info[,1] <- NULL

#save cluster info and data for later use
saveRDS(mfuzz.BN.9, file = "mfuzzBN9.RDS")
saveRDS(mfuzz.BA.7, file = "mfuzzBA7.RDS")
saveRDS(BN.gene.cluster.info, file = "mfuzzClusterInfoBN.RDS")
saveRDS(BA.gene.cluster.info, file = "mfuzzClusterInfoBA.RDS")

#make lists of genes in each interesting cluster 
BN.cluster5 <- rownames(subset(BN.gene.cluster.info, cluster == 5))
BN.cluster6 <- rownames(subset(BN.gene.cluster.info, cluster == 6))
BA.cluster5 <- rownames(subset(BA.gene.cluster.info, cluster == 5))

#extract expression data for clusters
BN.DEG.log.center.cluster5 <- subset(BN.DEG.log.center, rownames(BN.DEG.log.center) %in% BN.cluster5)
BN.DEG.log.center.cluster6 <- subset(BN.DEG.log.center, rownames(BN.DEG.log.center) %in% BN.cluster6)
BA.DEG.log.center.cluster5 <- subset(BA.DEG.log.center, rownames(BA.DEG.log.center) %in% BA.cluster5)

#plot heatmaps of these genes
make.heatmap.intGenes <- function(DEG, sample.anno, sample.hclust, gene.hclust, breaks.samp, anno.colors){
  pheatmap(DEG, cluster_rows = gene.hclust, cluster_cols = sample.hclust, show_colnames = FALSE, show_rownames = FALSE,
           color = cor.color, breaks = breaks.samp, 
           annotation_col = sample.anno, annotation_colors = anno.colors,
           annotation_names_col = FALSE, annotation_names_row = FALSE, legend_breaks = c(-6,-3,0,3,6), border_color = NA)
}

sample.anno.colors <- list(
  Group = c("H-DO,Imm." = "#2a9d8f", "H-DO,Acc." = "#e9c46a",
            "L-DO,Imm." = "#f4a261", "L-DO,Acc." = "#e76f51"))
  
BN.cluster5.breaks <- set.breaks(BN.DEG.log.center.cluster5)
BN.cluster6.breaks <- set.breaks(BN.DEG.log.center.cluster6)
BA.cluster5.breaks <- set.breaks(BA.DEG.log.center.cluster5)

tiff("BN_heatmap_cluster5.tiff", units="in", width = 7, height = 8, res = 600)
BN.DEG.heatmap.cluster5 <- make.heatmap.intGenes(BN.DEG.log.center.cluster5, sample.info.BN, BN.DEG.sample.hclust, FALSE, BN.cluster5.breaks, sample.anno.colors)
BN.DEG.heatmap.cluster5
dev.off()

tiff("BN_heatmap_cluster6.tiff", units="in", width = 7, height = 8, res = 600)
BN.DEG.heatmap.cluster6 <- make.heatmap.intGenes(BN.DEG.log.center.cluster6, sample.info.BN, BN.DEG.sample.hclust, FALSE, BN.cluster6.breaks, sample.anno.colors)
BN.DEG.heatmap.cluster6
dev.off()

tiff("BA_heatmap_cluster5.tiff", units="in", width = 7, height = 8, res = 600)
BA.DEG.heatmap.cluster5 <- make.heatmap.intGenes(BA.DEG.log.center.cluster5, sample.info.BA, BA.DEG.sample.hclust, FALSE, BA.cluster5.breaks, sample.anno.colors)
BA.DEG.heatmap.cluster5
dev.off()

#plot heatmaps of all DEGs with mfuzz cluster info 
#extract cluster info
BN.mfuzz.clusters <- BN.gene.cluster.info[,1, drop=FALSE]
BA.mfuzz.clusters <- BA.gene.cluster.info[,1, drop=FALSE]

#change name of cluster column 
colnames(BN.mfuzz.clusters) <- "Cluster"
colnames(BA.mfuzz.clusters) <- "Cluster"

#change names of clusters
BN.mfuzz.clusters$Cluster <- as.factor(BN.mfuzz.clusters$Cluster)
BN.mfuzz.clusters$Cluster <- gsub("\\<1\\>", "Cluster1", BN.mfuzz.clusters$Cluster)
BN.mfuzz.clusters$Cluster <- gsub("\\<2\\>", "Cluster2", BN.mfuzz.clusters$Cluster)
BN.mfuzz.clusters$Cluster <- gsub("\\<3\\>", "Cluster3", BN.mfuzz.clusters$Cluster)
BN.mfuzz.clusters$Cluster <- gsub("\\<4\\>", "Cluster4", BN.mfuzz.clusters$Cluster)
BN.mfuzz.clusters$Cluster <- gsub("\\<5\\>", "Cluster5", BN.mfuzz.clusters$Cluster)
BN.mfuzz.clusters$Cluster <- gsub("\\<6\\>", "Cluster6", BN.mfuzz.clusters$Cluster)
BN.mfuzz.clusters$Cluster <- gsub("\\<7\\>", "Cluster7", BN.mfuzz.clusters$Cluster)
BN.mfuzz.clusters$Cluster <- gsub("\\<8\\>", "Cluster8", BN.mfuzz.clusters$Cluster)
BN.mfuzz.clusters$Cluster <- gsub("\\<9\\>", "Cluster9", BN.mfuzz.clusters$Cluster)

BA.mfuzz.clusters$Cluster <- as.factor(BA.mfuzz.clusters$Cluster)
BA.mfuzz.clusters$Cluster <- gsub("\\<1\\>", "Cluster1", BA.mfuzz.clusters$Cluster)
BA.mfuzz.clusters$Cluster <- gsub("\\<2\\>", "Cluster2", BA.mfuzz.clusters$Cluster)
BA.mfuzz.clusters$Cluster <- gsub("\\<3\\>", "Cluster3", BA.mfuzz.clusters$Cluster)
BA.mfuzz.clusters$Cluster <- gsub("\\<4\\>", "Cluster4", BA.mfuzz.clusters$Cluster)
BA.mfuzz.clusters$Cluster <- gsub("\\<5\\>", "Cluster5", BA.mfuzz.clusters$Cluster)
BA.mfuzz.clusters$Cluster <- gsub("\\<6\\>", "Cluster6", BA.mfuzz.clusters$Cluster)
BA.mfuzz.clusters$Cluster <- gsub("\\<7\\>", "Cluster7", BA.mfuzz.clusters$Cluster)

#set annotation colours
BN.mfuzz.anno.colors <- list(
  Group = c("H-DO,Imm." = "#2a9d8f", "H-DO,Acc." = "#e9c46a",
            "L-DO,Imm." = "#f4a261", "L-DO,Acc." = "#e76f51"),
  Cluster = c(Cluster1 = "#ea5545", Cluster2 = "#f46a9b", Cluster3 = "#ffd500", Cluster4 = "#87bc45", Cluster5 =  "#006400",
              Cluster6 = "#27aeef", Cluster7 = "#00509d", Cluster8 = "#b33dc6", Cluster9 = "#5a189a"))

BA.mfuzz.anno.colors <- list(
  Group = c("H-DO,Imm." = "#2a9d8f", "H-DO,Acc." = "#e9c46a",
            "L-DO,Imm." = "#f4a261", "L-DO,Acc." = "#e76f51"),
  Cluster = c(Cluster1 = "#ea5545", Cluster2 = "#f46a9b", Cluster3 = "#ffd500", Cluster4 = "#87bc45", Cluster5 =  "#27aeef",
              Cluster6 = "#00509d", Cluster7 = "#b33dc6"))

#reorder expression data to be in order of clusters 
#add cluster info on
BN.DEG.log.center.mfuzz <- merge(BN.DEG.log.center, BN.mfuzz.clusters, by = "row.names", all = TRUE)
rownames(BN.DEG.log.center.mfuzz) <- BN.DEG.log.center.mfuzz[,1]
BN.DEG.log.center.mfuzz[,1] <- NULL

BA.DEG.log.center.mfuzz <- merge(BA.DEG.log.center, BA.mfuzz.clusters, by = "row.names", all = TRUE)
rownames(BA.DEG.log.center.mfuzz) <- BA.DEG.log.center.mfuzz[,1]
BA.DEG.log.center.mfuzz[,1] <- NULL

#order by cluster info 
BN.DEG.log.center.mfuzz <- arrange(BN.DEG.log.center.mfuzz, Cluster)
BA.DEG.log.center.mfuzz <- arrange(BA.DEG.log.center.mfuzz, Cluster)

#remove cluster column 
BN.DEG.log.center.mfuzz <- BN.DEG.log.center.mfuzz[,-33]
BA.DEG.log.center.mfuzz <- BA.DEG.log.center.mfuzz[,-38]

#plot
tiff("BN_heatmap_DEG_mfuzz.tiff", units="in", width = 7, height = 8, res = 600)
BN.DEG.heatmap.mfuzz <- make.heatmap(BN.DEG.log.center.mfuzz, sample.info.BN, BN.mfuzz.clusters, BN.DEG.sample.hclust, FALSE, BN.DEG.breaks, BN.mfuzz.anno.colors)
BN.DEG.heatmap.mfuzz
dev.off()

tiff("BA_heatmap_DEG_mfuxx.tiff", units="in", width = 7, height = 8, res = 600)
BA.DEG.heatmap.mfuzz <- make.heatmap(BA.DEG.log.center.mfuzz, sample.info.BA, BA.mfuzz.clusters, BA.DEG.sample.hclust, FALSE, BA.DEG.breaks, BA.mfuzz.anno.colors)
BA.DEG.heatmap.mfuzz
dev.off()

## plot heatmap of genes identified as interesting in venn diagram ##
#get genes of interest 
BN.int.gene <- Reduce(intersect, list(BN.StF.vs.SwF.DEG.list$gene,BN.SwF.vs.SwP.DEG.list$gene))
BN.int.gene <- setdiff(BN.int.gene, BN.StP.vs.SwP.DEG.list$gene)
BN.int.gene <- setdiff(BN.int.gene, BN.StF.vs.StP.DEG.list$gene)

BA.int.gene <- Reduce(intersect, list(BA.StF.vs.SwF.DEG.list$gene,BA.SwF.vs.SwP.DEG.list$gene))
BA.int.gene <- setdiff(BA.int.gene, BA.StP.vs.SwP.DEG.list$gene)
BA.int.gene <- setdiff(BA.int.gene, BA.StF.vs.StP.DEG.list$gene)

#extract expression data for clusters
BN.DEG.log.center.int.gene <- subset(BN.DEG.log.center, rownames(BN.DEG.log.center) %in% BN.int.gene)
BA.DEG.log.center.int.gene <- subset(BA.DEG.log.center, rownames(BA.DEG.log.center) %in% BA.int.gene)

#plot 
BN.intgene.breaks <- set.breaks(BN.DEG.log.center.int.gene)
BA.intgene.breaks <- set.breaks(BA.DEG.log.center.int.gene)

tiff("BN_heatmap_vennGens.tiff", units="in", width = 7, height = 8, res = 600)
BN.DEG.heatmap.int.gene <- make.heatmap.intGenes(BN.DEG.log.center.int.gene, sample.info.BN, BN.DEG.sample.hclust, TRUE, BN.intgene.breaks, sample.anno.colors)
BN.DEG.heatmap.int.gene
dev.off()

tiff("BA_heatmap_vennGens.tiff", units="in", width = 7, height = 8, res = 600)
BA.DEG.heatmap.int.gene <- make.heatmap.intGenes(BA.DEG.log.center.int.gene, sample.info.BA, BA.DEG.sample.hclust, TRUE, BA.intgene.breaks, sample.anno.colors)
BA.DEG.heatmap.int.gene
dev.off()

## Plot volcano plots ##
## add in information on up/down regulation ##
add.diffExpr <- function(data){
  data$diffExpr <- "NO"
  data$diffExpr[data$logFC >= 2 & data$FDR <= 0.01] <- "UP"
  data$diffExpr[data$logFC <= -2 & data$FDR <= 0.01] <- "DOWN"
  return(data)
}

BN.StF.vs.StP.result <- add.diffExpr(BN.StF.vs.StP.result)
BN.StF.vs.SwF.result <- add.diffExpr(BN.StF.vs.SwF.result)
BN.StP.vs.SwP.result <- add.diffExpr(BN.StP.vs.SwP.result)
BN.SwF.vs.SwP.result <- add.diffExpr(BN.SwF.vs.SwP.result)

BA.StF.vs.StP.result <- add.diffExpr(BA.StF.vs.StP.result)
BA.StF.vs.SwF.result <- add.diffExpr(BA.StF.vs.SwF.result)
BA.StP.vs.SwP.result <- add.diffExpr(BA.StP.vs.SwP.result)
BA.SwF.vs.SwP.result <- add.diffExpr(BA.SwF.vs.SwP.result)

# plot
plot.volcano <- function(data){
  ggplot(data = data, aes(x = logFC, y = -log2(FDR), col = diffExpr)) + geom_point(alpha = 0.7) + 
    theme_bw() + scale_color_manual(values =  c("#e76f51", "black", "#82DCDA")) + xlab("log2(FC)") + ylab("-log2(FDR)") +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), legend.position = "none")
}

tiff("BN_StF.vs.StP_volcano.tiff", units="in", width = 6, height = 6, res = 600)
BN.StF.vs.StP.volcano <- plot.volcano(BN.StF.vs.StP.result)
BN.StF.vs.StP.volcano
dev.off()

tiff("BN_StF.vs.SwF_volcano.tiff", units="in", width = 6, height = 6, res = 600)
BN.StF.vs.SwF.volcano <- plot.volcano(BN.StF.vs.SwF.result)
BN.StF.vs.SwF.volcano
dev.off()

tiff("BN_StP.vs.SwP_volcano.tiff", units="in", width = 6, height = 6, res = 600)
BN.StP.vs.SwP.volcano <- plot.volcano(BN.StP.vs.SwP.result)
BN.StP.vs.SwP.volcano
dev.off()

tiff("BN_SwF.vs.SwP_volcano.tiff", units="in", width = 6, height = 6, res = 600)
BN.SwF.vs.SwP.volcano <- plot.volcano(BN.SwF.vs.SwP.result)
BN.SwF.vs.SwP.volcano
dev.off()

tiff("BA_StF.vs.StP_volcano.tiff", units="in", width = 6, height = 6, res = 600)
BA.StF.vs.StP.volcano <- plot.volcano(BA.StF.vs.StP.result)
BA.StF.vs.StP.volcano
dev.off()

tiff("BA_StF.vs.SwF_volcano.tiff", units="in", width = 6, height = 6, res = 600)
BA.StF.vs.SwF.volcano <- plot.volcano(BA.StF.vs.SwF.result)
BA.StF.vs.SwF.volcano
dev.off()

tiff("BA_StP.vs.SwP_volcano.tiff", units="in", width = 6, height = 6, res = 600)
BA.StP.vs.SwP.volcano <- plot.volcano(BA.StP.vs.SwP.result)
BA.StP.vs.SwP.volcano
dev.off()

tiff("BA_SwF.vs.SwP_volcano.tiff", units="in", width = 6, height = 6, res = 600)
BA.SwF.vs.SwP.volcano <- plot.volcano(BA.SwF.vs.SwP.result)
BA.SwF.vs.SwP.volcano
dev.off()

## Make paneled figures ## 
tiff("PCA_panel.tiff", units="in", width = 10, height = 10, res = 600)
PCA.panel <- ggarrange(allExprBN.pca12.plot, allExprBA.pca12.plot, BN.DEG.pca12.plot, BA.DEG.pca12.plot,
                       labels = c("A", "B", "C", "D"),
                       ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

PCA.panel
dev.off()

tiff("PCA_pond_panel.tiff", units="in", width = 10, height = 10, res = 600)
pond.PCA.panel <- ggarrange(allExprBN.pca12.pond.plot, allExprBA.pca12.pond.plot, BN.DEG.pca12.pond.plot, BA.DEG.pca12.pond.plot,
                           labels = c("A", "B", "C", "D"),
                           ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
pond.PCA.panel
dev.off()

tiff("BN_volcano_panel.tiff", units="in", width = 10, height = 10, res = 600)
BN.volcano.plot <- ggarrange(BN.StF.vs.SwF.volcano, BN.StP.vs.SwP.volcano, BN.StF.vs.StP.volcano, BN.SwF.vs.SwP.volcano,
                             labels = c("A", "B", "C", "D"),
                             ncol = 2, nrow = 2, common.legend = FALSE)

BN.volcano.plot
dev.off()

tiff("BA_volcano_panel.tiff", units="in", width = 10, height = 10, res = 600)
BA.volcano.plot <- ggarrange(BA.StF.vs.SwF.volcano, BA.StP.vs.SwP.volcano, BA.StF.vs.StP.volcano, BA.SwF.vs.SwP.volcano,
                             labels = c("A", "B", "C", "D"),
                             ncol = 2, nrow = 2, common.legend = FALSE)

BA.volcano.plot
dev.off()

tiff("compareDEG_panel.tiff", units="in", width = 10, height = 10, res = 600)
num.DEG.plot <- ggarrange(plot.compareDEG, 
                          ggarrange(BN.upset, BA.upset, labels = c("B", "C"), ncol = 2), 
                          nrow = 2, labels = "A")

num.DEG.plot
dev.off()


tiff("compareDEG_panel_sep.tiff", units="in", width = 11, height = 10, res = 600)
num.DEG.plot <- ggarrange(plot.compareDEG.BN, plot.compareDEG.BA,
                          BN.upset, BA.upset, 
                          labels = c("A", "B", "C", "D"),
                          ncol = 2, nrow = 2)

num.DEG.plot
dev.off()


## Make lists of genes to analyze for GO analysis ## 
#make factor lists (genes to analyze)
BN.venn.genes <- data.frame(factor = "BN_vennGenes", gene_id = BN.int.gene)
BN.k.cluster2.genes <- data.frame(factor = "BN_k_cluster2", gene_id = BN.k.cluster2)
BN.30perc.cluster7.genes <- data.frame(factor = "BN_30perc_cluster7", gene_id = BN.30perc.cluster7)
BN.mfuzz.cluster5.genes <- data.frame(factor = "BN_mfuzz_cluster5", gene_id = BN.cluster5)
BN.mfuzz.cluster6.genes <- data.frame(factor = "BN_mfuzz_cluster6", gene_id = BN.cluster6)
BN.SwF.vs.SwP.up.genes <- data.frame(factor = "BN_SwFvsSwP_DEG_up", gene_id = rownames(BN.SwF.vs.SwP.DEG[BN.SwF.vs.SwP.DEG$logFC > 0,]))
BN.SwF.vs.SwP.down.genes <- data.frame(factor = "BN_SwFvsSwP_DEG_down", gene_id = rownames(BN.SwF.vs.SwP.DEG[BN.SwF.vs.SwP.DEG$logFC < 0,]))
BN.StF.vs.SwF.up.genes <- data.frame(factor = "BN_StFvsSwF_DEG_up", gene_id = rownames(BN.StF.vs.SwF.DEG[BN.StF.vs.SwF.DEG$logFC > 0,]))
BN.StF.vs.SwF.down.genes <- data.frame(factor = "BN_StFvsSwF_DEG_down", gene_id = rownames(BN.StF.vs.SwF.DEG[BN.StF.vs.SwF.DEG$logFC < 0,]))


BA.venn.genes <- data.frame(factor = "BA_vennGenes", gene_id = BA.int.gene)
BA.k.cluster3.genes <- data.frame(factor = "BA_k_cluster3", gene_id = BA.k.cluster3)
BA.30perc.cluster5.genes <- data.frame(factor = "BA_30perc_cluster5", gene_id = BA.30perc.cluster5)
BA.mfuzz.cluster5.genes <- data.frame(factor = "BA_mfuzz_cluster5", gene_id = BA.cluster5)
BA.SwF.vs.SwP.up.genes <- data.frame(factor = "BA_SwFvsSwP_DEG_up", gene_id = rownames(BA.SwF.vs.SwP.DEG[BA.SwF.vs.SwP.DEG$logFC > 0,]))
BA.SwF.vs.SwP.down.genes <- data.frame(factor = "BA_SwFvsSwP_DEG_down", gene_id = rownames(BA.SwF.vs.SwP.DEG[BA.SwF.vs.SwP.DEG$logFC < 0,]))
BA.StF.vs.SwF.up.genes <- data.frame(factor = "BA_StFvsSwF_DEG_up", gene_id = rownames(BA.StF.vs.SwF.DEG[BA.StF.vs.SwF.DEG$logFC > 0,]))
BA.StF.vs.SwF.down.genes <- data.frame(factor = "BA_StFvsSwF_DEG_down", gene_id = rownames(BA.StF.vs.SwF.DEG[BA.StF.vs.SwF.DEG$logFC < 0,]))

#make background lists 
BN.StF.vs.StP.background <- rownames(BN.StF.vs.StP.result) 
BN.StF.vs.SwF.background <- rownames(BN.StF.vs.SwF.result) 
BN.StP.vs.SwP.background <- rownames(BN.StP.vs.SwP.result) 
BN.SwF.vs.SwP.background <- rownames(BN.SwF.vs.SwP.result)

BN.background <- c(BN.StF.vs.StP.background, BN.StF.vs.SwF.background, BN.StP.vs.SwP.background, BN.SwF.vs.SwP.background)
BN.background <- unique(BN.background)

BA.StF.vs.StP.background <- rownames(BA.StF.vs.StP.result) 
BA.StF.vs.SwF.background <- rownames(BA.StF.vs.SwF.result) 
BA.StP.vs.SwP.background <- rownames(BA.StP.vs.SwP.result) 
BA.SwF.vs.SwP.background <- rownames(BA.SwF.vs.SwP.result) 

BA.background <- c(BA.StF.vs.StP.background, BA.StF.vs.SwF.background, BA.StP.vs.SwP.background, BA.SwF.vs.SwP.background)
BA.background <- unique(BA.background)

#save gene lists 
write.table(BN.venn.genes, "./data/goSeqData/BN_venn_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BN.30perc.cluster7.genes, "./data/goSeqData/BN_30perc_cluster7_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BN.k.cluster2.genes, "./data/goSeqData/BN_k_cluster2_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BN.mfuzz.cluster5.genes, "./data/goSeqData/BN_mfuzz_cluster5_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BN.mfuzz.cluster6.genes, "./data/goSeqData/BN_mfuzz_cluster6_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BN.SwF.vs.SwP.up.genes, "./data/goSeqData/BN_SwF.vs.SwP_up_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BN.SwF.vs.SwP.down.genes, "./data/goSeqData/BN_SwF.vs.SwP_down_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BN.StF.vs.SwF.up.genes, "./data/goSeqData/BN_StF.vs.SwF_up_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BN.StF.vs.SwF.down.genes, "./data/goSeqData/BN_StF.vs.SwF_down_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(BA.venn.genes, "./data/goSeqData/BA_venn_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BA.30perc.cluster5.genes, "./data/goSeqData/BA_30perc_cluster5_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BA.k.cluster3.genes, "./data/goSeqData/BA_k_cluster3_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BA.mfuzz.cluster5.genes, "./data/goSeqData/BA_mfuzz_cluster5_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BA.SwF.vs.SwP.up.genes, "./data/goSeqData/BA_SwF.vs.SwP_up_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BA.SwF.vs.SwP.down.genes, "./data/goSeqData/BA_SwF.vs.SwP_down_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BA.StF.vs.SwF.up.genes, "./data/goSeqData/BA_StF.vs.SwF_up_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BA.StF.vs.SwF.down.genes, "./data/goSeqData/BA_StF.vs.SwF_down_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(BN.background, "./data/goSeqData/BN_background.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BA.background, "./data/goSeqData/BA_background.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

