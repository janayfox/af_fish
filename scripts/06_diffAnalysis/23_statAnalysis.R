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
library(reshape2)

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

## read in functions ## 
#function to extract expression of DEGs 
extractDEG <- function(DEG.df, DEG.list){
  return(DEG.df[DEG.list,])
}

#function to log transform data
log.trans <- function(df){
  return(log((df + 1), 2))
}

#function to median center data
med.center <- function(df){
  rowMed <- apply(df, 1, median)
  return(df-rowMed)
}

#function to transpose data
transpose.names <- function(data){
  t.data <- t(data)
  rownames(t.data) <- colnames(data)
  colnames(t.data) <- rownames(data)
  return(as.data.frame(t.data))
}

#functions to add on site and collection data for each species
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

#function to plot PCAs
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
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16), legend.text = element_text(size=18), 
          legend.title = element_text(size = 18))
}

# function to plot PCA with pond information 
plot.pond.pca <- function(pca.data, pca.x, pca.y, xlabel, ylabel){
  ggplot(data = pca.data, aes(x = pca.x, y = pca.y, color = pond, shape = group)) + theme_bw() + geom_point() +
    scale_color_manual(values = c("#ad2e24", "#27187E", "#82DCDA")) + 
    scale_shape_manual(values = c(16, 15, 17, 4), labels = c("H-I", "H-A", 
                                                             "L-I", "L-A")) +
    labs(color = "Pond", shape = "Group") + geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") + stat_ellipse(aes(group = pond)) + xlab(xlabel) + ylab(ylabel) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), legend.text = element_text(size=11))
}

#function to plot mFuzz clusters for both species
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

#function to set breaks for heatmaps
set.breaks <- function(expr.data){
  new.breaks <- breaks
  new.breaks[length(breaks)] <- max(max(expr.data),max(breaks))
  new.breaks[1] <- min(min(expr.data),min(breaks))
  return(new.breaks)
}

#function to plot heatmaps
make.heatmap <- function(DEG, sample.anno, gene.anno, sample.hclust, gene.hclust, breaks.samp, anno.colors){
  pheatmap(DEG, cluster_rows = gene.hclust, cluster_cols = sample.hclust, show_colnames = FALSE, show_rownames = FALSE,
           color = cor.color, breaks = breaks.samp, annotation_row = gene.anno,
           annotation_col = sample.anno, annotation_colors = anno.colors,
           annotation_names_col = FALSE, annotation_names_row = FALSE, legend_breaks = c(-4,-2,0,2,4), border_color = NA)
}

# function to add in information on up/down regulation
add.diffExpr <- function(data){
  data$diffExpr <- "NO"
  data$diffExpr[data$logFC >= 2 & data$FDR <= 0.01] <- "UP"
  data$diffExpr[data$logFC <= -2 & data$FDR <= 0.01] <- "DOWN"
  return(data)
}

#function for plotting volcano plots
plot.volcano <- function(data){
  ggplot(data = data, aes(x = logFC, y = -log2(FDR), col = diffExpr)) + geom_point(alpha = 0.7) + 
    theme_bw() + scale_color_manual(values =  c("#e76f51", "black", "#82DCDA")) + xlab("log2(FC)") + ylab("-log2(FDR)") +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), legend.position = "none")
}

## Check number of genes retained in DEG analysis ## 
BA_genes <- c(rownames(BA.StF.vs.StP.result), rownames(BA.StP.vs.SwP.result), rownames(BA.SwF.vs.SwP.result))
BA_genes <- unique(BA_genes)

BN_genes <- c(rownames(BN.StF.vs.StP.result), rownames(BN.StP.vs.SwP.result), rownames(BN.SwF.vs.SwP.result))
BN_genes <- unique(BN_genes)

##  Make UpSet plots and venn diagrams ## 
#get dataframe in proper format 
BN.StF.vs.StP.DEG.list <- data.frame(gene = rownames(BN.StF.vs.StP.DEG), comparison = "H-I vs H-A")
BN.StP.vs.SwP.DEG.list <- data.frame(gene = rownames(BN.StP.vs.SwP.DEG), comparison = "H-A vs L-A")
BN.SwF.vs.SwP.DEG.list <- data.frame(gene = rownames(BN.SwF.vs.SwP.DEG), comparison = "L-I vs L-A")

BN.DEG.upset <- rbind(BN.StF.vs.StP.DEG.list, BN.StP.vs.SwP.DEG.list, BN.SwF.vs.SwP.DEG.list)
BN.DEG.upset <- BN.DEG.upset %>% group_by(gene) %>% summarize(comparisons = list(comparison))

BA.StF.vs.StP.DEG.list <- data.frame(gene = rownames(BA.StF.vs.StP.DEG), comparison = "H-I vs H-A")
BA.StP.vs.SwP.DEG.list <- data.frame(gene = rownames(BA.StP.vs.SwP.DEG), comparison = "H-A vs L-A")
BA.SwF.vs.SwP.DEG.list <- data.frame(gene = rownames(BA.SwF.vs.SwP.DEG), comparison = "L-I vs L-A")

BA.DEG.upset <- rbind(BA.StF.vs.StP.DEG.list, BA.StP.vs.SwP.DEG.list, BA.SwF.vs.SwP.DEG.list)
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

## Run and Plot PCA on DEG##
#make list of DEG from all four important comparisons 
BN.DEG.list <- c(BN.StF.vs.StP.DEG.list$gene,BN.SwF.vs.SwP.DEG.list$gene, BN.StP.vs.SwP.DEG.list$gene)
BN.DEG.list <- unique(BN.DEG.list) #remove duplicates

BA.DEG.list <- c(BA.StF.vs.StP.DEG.list$gene, BA.SwF.vs.SwP.DEG.list$gene, BA.StP.vs.SwP.DEG.list$gene)
BA.DEG.list <- unique(BA.DEG.list) #remove duplicates

#extract DEG
BN.DEG <- extractDEG(allExprBN, BN.DEG.list)
BA.DEG <- extractDEG(allExprBA, BA.DEG.list)

#log2 transform data 
BN.DEG.log <- log.trans(BN.DEG)
BA.DEG.log <- log.trans(BA.DEG)

#median center data 
BN.DEG.log.center <- med.center(BN.DEG.log)
BA.DEG.log.center <- med.center(BA.DEG.log)

#transpose data 
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
t.BN.DEG.log.center <- add.site.collection.BN(t.BN.DEG.log.center)
t.BA.DEG.log.center <- add.site.collection.BA(t.BA.DEG.log.center)

#plot PCA
tiff("BN_pca12_DEG_filtered.tiff", units="in", width = 8, height = 6, res = 600)
BN.DEG.pca12.plot <- plot.pca(t.BN.DEG.log.center, t.BN.DEG.log.center$PC1, 
                                  t.BN.DEG.log.center$PC2, "PC1 (46.2%)", "PC2 (11.6%)")
BN.DEG.pca12.plot
dev.off()

tiff("BN_pca23_DEG.tiff", units="in", width = 8, height = 6, res = 600)
BN.DEG.pca23.plot <- plot.pca(t.BN.DEG.log.center, t.BN.DEG.log.center$PC2, 
                                  t.BN.DEG.log.center$PC3, "PC2 (11.6%)", "PC3 (9.5%)")
BN.DEG.pca23.plot
dev.off()

tiff("BA_pca12_DEG_filtered.tiff", units="in", width = 8, height = 6, res = 600)
BA.DEG.pca12.plot <- plot.pca(t.BA.DEG.log.center, t.BA.DEG.log.center$PC1, 
                                  t.BA.DEG.log.center$PC2, "PC1 (66.1%)", "PC2 (13.1%)")
BA.DEG.pca12.plot
dev.off()

tiff("BA_pca23_DEG.tiff", units="in", width = 8, height = 6, res = 600)
BA.DEG.pca23.plot <- plot.pca(t.BA.DEG.log.center, t.BA.DEG.log.center$PC2, 
                                  t.BA.DEG.log.center$PC3, "PC2 (13.1%)", "PC3 (5.5%)")
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

#plot with pond info
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
tiff("BN_pca12_allExpr_filtered.tiff", units="in", width = 8, height = 6, res = 600)
allExprBN.pca12.plot <- plot.pca(t.allExprBN.log.center, t.allExprBN.log.center$PC1, t.allExprBN.log.center$PC2, "PC1 (9.6%)", "PC2 (5.8%)")
allExprBN.pca12.plot
dev.off()

tiff("BN_pca23_allExpr.tiff", units="in", width = 8, height = 6, res = 600)
allExprBN.pca23.plot <- plot.pca(t.allExprBN.log.center, t.allExprBN.log.center$PC2, t.allExprBN.log.center$PC3, "PC2 (5.8%)", "PC3 (4.8%)")
allExprBN.pca23.plot
dev.off()

tiff("BA_pca12_allExpr_filtered.tiff", units="in", width = 8, height = 6, res = 600)
allExprBA.pca12.plot <- plot.pca(t.allExprBA.log.center, t.allExprBA.log.center$PC1, t.allExprBA.log.center$PC2, "PC1 (18.9%)", "PC2 (7.7%)")
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

## check for effect of pond on overall gene expression ##
#convert dataframe to correct format for analysis 
BA.matrix <- t.allExprBA.log.center[t.allExprBA.log.center$pond != "Field",] #get only pond samples
BN.matrix <- t.allExprBN.log.center[t.allExprBN.log.center$pond != "Field",]

# check for effect of pond on PCs 
summary(aov(PC1 ~ pond, data = BA.matrix))
summary(aov(PC2 ~ pond, data = BA.matrix))
summary(aov(PC3 ~ pond, data = BA.matrix))
summary(aov(PC4 ~ pond, data = BA.matrix))

summary(aov(PC1 ~ pond, data = BN.matrix))
summary(aov(PC2 ~ pond, data = BN.matrix))
summary(aov(PC3 ~ pond, data = BN.matrix))
summary(aov(PC4 ~ pond, data = BN.matrix))

#convert dataframe to correct format for analysis 
BA.matrix <- t.BA.DEG.log.center[t.BA.DEG.log.center$pond != "Field",] #get only pond samples
BN.matrix <- t.BN.DEG.log.center[t.BN.DEG.log.center$pond != "Field",]

col.drop <- c("PC1", "PC2", "PC3", "PC4", "collection", "group") #drop other variables
BA.matrix <- BA.matrix[, !(names(BA.matrix) %in% col.drop)]
BN.matrix <- BN.matrix[, !(names(BN.matrix) %in% col.drop)]

BA.matrix <- rownames_to_column(BA.matrix, var = "individual")
BN.matrix <- rownames_to_column(BN.matrix, var = "individual")

melt_BA <- melt(BA.matrix, id.vars = c("pond", "site", "individual"), variable.name = "gene_ID", value.name = "gene_expression")
melt_BN <- melt(BN.matrix, id.vars = c("pond", "site", "individual"), variable.name = "gene_ID", value.name = "gene_expression")

#fit model 
fit.BA <- lm(gene_expression ~ gene_ID + individual + site + pond, data = melt_BA)
fit.BN <- lm(gene_expression ~ gene_ID + individual + site + pond, data = melt_BN)

an.res.BA <- anova(fit.BA)
an.res.BN <- anova(fit.BN)

an.res.BA
an.res.BN

## Plot Bar Graphs comparing numbers of DEGS ## 
#construct dataframe 
comparison <- c("H-I vs H-A", "L-I vs L-A", "H-A vs L-A",
                    "H-I vs H-A", "L-I vs L-A", "H-A vs L-A")

species <- c("EN", "EN", "EN", "EA", "EA", "EA")

numDEG <- c(nrow(BN.StF.vs.StP.DEG), nrow(BN.SwF.vs.SwP.DEG), nrow(BN.StP.vs.SwP.DEG), 
                nrow(BA.StF.vs.StP.DEG), nrow(BA.SwF.vs.SwP.DEG), nrow(BA.StP.vs.SwP.DEG))

compareDEG <- data.frame(numDEG, species, comparison)

#order factor in order I want in barplot
compareDEG$comparison <- factor(compareDEG$comparison, levels = c("H-I vs H-A", "L-I vs L-A", "H-A vs L-A"))

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

plot.compareDEG.BN<- ggplot(data = compareDEG.BN, aes(x = comparison, y = numDEG)) + geom_bar(stat = "identity", position = position_dodge(), fill = "#27187E") + 
  geom_text(aes(label=numDEG), vjust=1.6, color="white",position = position_dodge(0.9), size=3.5)+
  theme_bw() + 
  xlab("Comparison Type") + ylab ("Number of DEGs") + labs(fill = "Species") +
  theme(legend.position = "bottom", axis.title = element_text(size = 14), axis.text = element_text(size = 11),
        legend.text = element_text(size=11)) + scale_x_discrete(labels = function(x) str_wrap(x, width = 5))

plot.compareDEG.BA<- ggplot(data = compareDEG.BA, aes(x = comparison, y = numDEG)) + geom_bar(stat = "identity", position = position_dodge(), fill = "#27187E") + 
  geom_text(aes(label=numDEG), vjust=1.6, color="white",position = position_dodge(0.9), size=3.5)+
  theme_bw() + 
  xlab("Comparison Type") + ylab ("Number of DEGs") + labs(fill = "Species") +
  theme(legend.position = "bottom", axis.title = element_text(size = 14), axis.text = element_text(size = 11),
        legend.text = element_text(size=11)) + scale_x_discrete(labels = function(x) str_wrap(x, width = 5))

tiff("BA_barplot.tiff", units="in", width = 6, height = 7, res = 600)
plot.compareDEG.BA
dev.off()

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
  Group = c("H-I" = "#2a9d8f", "H-A" = "#e9c46a",
            "L-I" = "#f4a261", "L-A" = "#e76f51"),
  Cluster = c(Cluster1 = "#ea5545", Cluster2 = "#f46a9b", Cluster3 = "#ffd500", Cluster4 = "#87bc45", Cluster5 =  "#006400",
              Cluster6 = "#27aeef", Cluster7 = "#00509d", Cluster8 = "#b33dc6", Cluster9 = "#5a189a"))

BA.mfuzz.anno.colors <- list(
  Group = c("H-I" = "#2a9d8f", "H-A" = "#e9c46a",
            "L-I" = "#f4a261", "L-A" = "#e76f51"),
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

#do cluster analysis on samples 
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

#generate color pallette 
cor.color <- colorRampPalette(c("#FFFF00", "black", "#0000FF"))(60) 

#remove NAs
BN.DEG.log.center.mfuzz <- na.omit(BN.DEG.log.center.mfuzz)
BA.DEG.log.center.mfuzz <- na.omit(BA.DEG.log.center.mfuzz)

#calculate breaks 
breaks = c(seq(-4, 0, length.out = 30),seq(0.01, 4, length.out = 30)) #changed to 4 to fix ramping of colors

#set breaks
BN.DEG.breaks <- set.breaks(BN.DEG.log.center)
BA.DEG.breaks <- set.breaks(BA.DEG.log.center)

#plot
tiff("BN_heatmap_DEG_mfuzz.tiff", units="in", width = 7, height = 8, res = 600)
BN.DEG.heatmap.mfuzz <- make.heatmap(BN.DEG.log.center.mfuzz, sample.info.BN, BN.mfuzz.clusters, BN.DEG.sample.hclust, FALSE, BN.DEG.breaks, BN.mfuzz.anno.colors)
BN.DEG.heatmap.mfuzz
dev.off()

tiff("BA_heatmap_DEG_mfuzz_4.tiff", units="in", width = 7, height = 8, res = 600)
BA.DEG.heatmap.mfuzz <- make.heatmap(BA.DEG.log.center.mfuzz, sample.info.BA, BA.mfuzz.clusters, BA.DEG.sample.hclust, FALSE, BA.DEG.breaks, BA.mfuzz.anno.colors)
BA.DEG.heatmap.mfuzz
dev.off()

## Plot volcano plots ##
#add on up/down info
BN.StF.vs.StP.result <- add.diffExpr(BN.StF.vs.StP.result)
BN.StP.vs.SwP.result <- add.diffExpr(BN.StP.vs.SwP.result)
BN.SwF.vs.SwP.result <- add.diffExpr(BN.SwF.vs.SwP.result)

BA.StF.vs.StP.result <- add.diffExpr(BA.StF.vs.StP.result)
BA.StP.vs.SwP.result <- add.diffExpr(BA.StP.vs.SwP.result)
BA.SwF.vs.SwP.result <- add.diffExpr(BA.SwF.vs.SwP.result)

# plot
tiff("BN_StF.vs.StP_volcano.tiff", units="in", width = 6, height = 6, res = 600)
BN.StF.vs.StP.volcano <- plot.volcano(BN.StF.vs.StP.result)
BN.StF.vs.StP.volcano
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
BN.volcano.plot <- ggarrange(BN.StP.vs.SwP.volcano, BN.StF.vs.StP.volcano, BN.SwF.vs.SwP.volcano,
                             labels = c("A", "B", "C"),
                             ncol = 2, nrow = 2, common.legend = FALSE)

BN.volcano.plot
dev.off()

tiff("BA_volcano_panel.tiff", units="in", width = 10, height = 10, res = 600)
BA.volcano.plot <- ggarrange(BA.StP.vs.SwP.volcano, BA.StF.vs.StP.volcano, BA.SwF.vs.SwP.volcano,
                             labels = c("A", "B", "C"),
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
BN.mfuzz.cluster5.genes <- data.frame(factor = "BN_mfuzz_cluster5", gene_id = BN.cluster5)
BN.mfuzz.cluster6.genes <- data.frame(factor = "BN_mfuzz_cluster6", gene_id = BN.cluster6)

BA.mfuzz.cluster5.genes <- data.frame(factor = "BA_mfuzz_cluster5", gene_id = BA.cluster5)

#make background lists 
BN.StF.vs.StP.background <- rownames(BN.StF.vs.StP.result) 
BN.StP.vs.SwP.background <- rownames(BN.StP.vs.SwP.result) 
BN.SwF.vs.SwP.background <- rownames(BN.SwF.vs.SwP.result)

BN.background <- c(BN.StF.vs.StP.background, BN.StP.vs.SwP.background, BN.SwF.vs.SwP.background)
BN.background <- unique(BN.background)

BA.StF.vs.StP.background <- rownames(BA.StF.vs.StP.result) 
BA.StP.vs.SwP.background <- rownames(BA.StP.vs.SwP.result) 
BA.SwF.vs.SwP.background <- rownames(BA.SwF.vs.SwP.result) 

BA.background <- c(BA.StF.vs.StP.backgroun, BA.StP.vs.SwP.background, BA.SwF.vs.SwP.background)
BA.background <- unique(BA.background)

#save gene lists 
write.table(BN.mfuzz.cluster5.genes, "./data/goSeqData/BN_mfuzz_cluster5_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BN.mfuzz.cluster6.genes, "./data/goSeqData/BN_mfuzz_cluster6_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(BA.k.cluster3.genes, "./data/goSeqData/BA_k_cluster3_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(BN.background, "./data/goSeqData/BN_background.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BA.background, "./data/goSeqData/BA_background.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)


## Calculate variance for each group ##
#format data frame for IQR calculation 
t.BN.DEG <- transform(merge(t(BN.DEG), sample.info.BN, by=0), row.names = Row.names, Row.names=NULL)
t.BA.DEG <- transform(merge(t(BA.DEG), sample.info.BA, by=0), row.names = Row.names, Row.names=NULL)

t.BN.allEXPR <- transform(merge(t(allExprBN), sample.info.BN, by=0), row.names = Row.names, Row.names=NULL)
t.BA.allEXPR <- transform(merge(t(allExprBA), sample.info.BA, by=0), row.names = Row.names, Row.names=NULL)

#calculate IQR for each DEG 
BN.IQR <- sapply(t.BN.DEG[,1:887], function(x) tapply(x, t.BN.DEG$Group, IQR))
BA.IQR <- sapply(t.BA.DEG[,1:776], function(x) tapply(x, t.BA.DEG$Group, IQR))

BN.IQR.all <- sapply(t.BN.allEXPR[,1:611157], function(x) tapply(x, t.BN.allEXPR$Group, IQR))
BA.IQR.all <- sapply(t.BA.allEXPR[,1:546320], function(x) tapply(x, t.BA.allEXPR$Group, IQR))

#edit for plotting 
BN.IQR <- as.data.frame(BN.IQR)
BA.IQR <- as.data.frame(BA.IQR)

BN.IQR$Group <- factor(row.names(BN.IQR))
BA.IQR$Group <- factor(row.names(BA.IQR))

BN.IQR.long <- gather(BN.IQR, gene, IQR, 1:887, factor_key = TRUE)
BA.IQR.long <- gather(BA.IQR, gene, IQR, 1:776, factor_key = TRUE)

BN.IQR.long.sub <- BN.IQR.long[!(BN.IQR.long$gene == "TRINITY_DN427527_c0_g1"),]

#plot boxplots 
ggplot(BN.IQR.long.sub, aes(x=Group, y = IQR)) + geom_boxplot()
ggplot(BA.IQR.long, aes(x=Group, y = IQR)) + geom_boxplot()

#get average IQRs
aggregate(BN.IQR.long$IQR, list(BN.IQR.long$Group), FUN=mean)
aggregate(BA.IQR.long$IQR, list(BA.IQR.long$Group), FUN=mean)

#calculate coefficient of variation 
BN.cv <- sapply(t.BN.DEG[,1:887], function(x) tapply(x, t.BN.DEG$Group, sd(x)/mean(x)))

cv <- function(x) sd(x) / mean(x) 

BN.DEG.long <- gather(t.BN.DEG, gene, expr, 1:887, factor_key = TRUE)
BA.DEG.long <- gather(t.BA.DEG, gene, expr, 1:776, factor_key = TRUE)

CV.BN <- BN.DEG.long %>% group_by(Group, gene) %>%
            summarise(cv = sd(expr)/mean(expr))

CV.BA <- BA.DEG.long %>% group_by(Group, gene) %>%
  summarise(cv = sd(expr)/mean(expr))

#calculate average coefficent of variation for each species 
mean(CV.BA$cv, na.rm = TRUE)
mean(CV.BN$cv, na.rm = TRUE)

#calculate average for populations 
CV.BN$pop <- substr(CV.BN$Group, 1, 1)
CV.BA$pop <- substr(CV.BA$Group, 1, 1)

CV.BA %>% group_by(pop) %>% 
  summarise(mean(cv, na.rm = TRUE))

CV.BN %>% group_by(pop) %>% 
  summarise(mean(cv, na.rm = TRUE))
