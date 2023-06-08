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
library(UpSetR)
library(VennDiagram)
library(RColorBrewer)
library(GGally)
library(Mfuzz)
library(dplyr)

#read in data 
degBN <- read.table("./data/DEG/BN/salmon/diffExpr.P0.01_C2.matrix")
degBA <- read.table("./data/DEG/BA/salmon/diffExpr.P0.01_C2.matrix")

allExprBN <- read.table("./data/totalExpr/BN_bf_new_sal.gene.TMM.EXPR.matrix")
allExprBA <- read.table("./data/totalExpr/BA_bf_new_sal.gene.TMM.EXPR.matrix")

allExprBN.TPM <- read.table("./data/totalExpr/BN.TMM.TPM.matrix")
allExprBA.TPM <- read.table("./data/totalExpr/BA.TMM.TPM.matrix")

allExprBN.RPKM <- read.table("./data/totalExpr/BN.TMM.RPKM.matrix")
allExprBA.RPKM <- read.table("./data/totalExpr/BA.TMM.RPKM.matrix")

BN.StF.vs.StP.DEG <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.stream_field_vs_stream_pond.edgeR.DE_results.P0.01_C2.DE.subset")
BN.StF.vs.SwF.DEG <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.stream_field_vs_swamp_field.edgeR.DE_results.P0.01_C2.DE.subset")
BN.StP.vs.SwP.DEG <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.stream_pond_vs_swamp_pond.edgeR.DE_results.P0.01_C2.DE.subset")
BN.SwF.vs.SwP.DEG <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.swamp_field_vs_swamp_pond.edgeR.DE_results.P0.01_C2.DE.subset")

BA.StF.vs.StP.DEG <- read.table("./data/DEG/BA/salmon/BA_bf_new_sal.gene.counts.matrix.stream_field_vs_stream_pond.edgeR.DE_results.P0.01_C2.DE.subset")
BA.StF.vs.SwF.DEG <- read.table("./data/DEG/BA/salmon/BA_bf_new_sal.gene.counts.matrix.stream_field_vs_swamp_field.edgeR.DE_results.P0.01_C2.DE.subset")
BA.StP.vs.SwP.DEG <- read.table("./data/DEG/BA/salmon/BA_bf_new_sal.gene.counts.matrix.stream_pond_vs_swamp_pond.edgeR.DE_results.P0.01_C2.DE.subset")
BA.SwF.vs.SwP.DEG <- read.table("./data/DEG/BA/salmon/BA_bf_new_sal.gene.counts.matrix.swamp_field_vs_swamp_pond.edgeR.DE_results.P0.01_C2.DE.subset")

BN.pond.data <- read.csv("./data/pond_data_BN.csv")
rownames(BN.pond.data) <- BN.pond.data[,1] #convert first column to row names 
BN.pond.data[,1] <- NULL

BA.pond.data <- read.csv("./data/pond_data_BA.csv")
rownames(BA.pond.data) <- BA.pond.data[,1] 
BA.pond.data[,1] <- NULL

##  Make UpSet plots and venn diagrams ## 
#get list of DEG
BN.StF.vs.StP.DEG.list <- rownames(BN.StF.vs.StP.DEG)
BN.StF.vs.SwF.DEG.list <- rownames(BN.StF.vs.SwF.DEG)
BN.StP.vs.SwP.DEG.list <- rownames(BN.StP.vs.SwP.DEG)
BN.SwF.vs.SwP.DEG.list <- rownames(BN.SwF.vs.SwP.DEG)

BN.DEG.list.upset <- list("H.DO-I vs H.DO-A" = BN.StF.vs.StP.DEG.list, 
                    "L.DO-I vs L.DO-A" = BN.SwF.vs.SwP.DEG.list,
                    "H.DO-I vs L.DO-I" = BN.StF.vs.SwF.DEG.list,
                    "H.DO-A vs L.DO-A" = BN.StP.vs.SwP.DEG.list)

BA.StF.vs.StP.DEG.list <- rownames(BA.StF.vs.StP.DEG)
BA.StF.vs.SwF.DEG.list <- rownames(BA.StF.vs.SwF.DEG)
BA.StP.vs.SwP.DEG.list <- rownames(BA.StP.vs.SwP.DEG)
BA.SwF.vs.SwP.DEG.list <- rownames(BA.SwF.vs.SwP.DEG)

BA.DEG.list.upset <- list("H.DO-I vs H.DO-A" = BA.StF.vs.StP.DEG.list, 
                    "L.DO-I vs L.DO-A" = BA.SwF.vs.SwP.DEG.list,
                    "H.DO-I vs L.DO-I" = BA.StF.vs.SwF.DEG.list,
                    "H.DO-A vs L.DO-A" = BA.StP.vs.SwP.DEG.list)

#Plot upset 
BN.upset <- upset(fromList(BN.DEG.list.upset), nsets = 4, order.by = "freq", main.bar.color = "#27187E",
                      matrix.color = "#0081a7", sets.bar.color = "#2a9d8f", text.scale = c(1.7,1.5,1.7,1.5,1.5,1.2))
BN.upset

BA.upset <- upset(fromList(BA.DEG.list.upset), nsets = 4, order.by = "freq", main.bar.color = "#27187E",
                      matrix.color = "#0081a7", sets.bar.color = "#2a9d8f", text.scale = c(1.7,1.5,1.7,1.5,1.5,1.2))
BA.upset

#Plot venn diagram 
#function to display venn diagram 
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

display_venn(BN.DEG.list, lwd = 2, col = c("#e76f51", "#0081a7", "#FFBC42", "#a7c957"),
             fill = c(alpha("#e76f51", 0.8), alpha("#0081a7", 0.8), alpha("#FFBC42",0.8), alpha("#a7c957",0.8)),
             fontface = "bold", cex = 1.5, cat.cex = 1.5)

display_venn(BA.DEG.list, lwd = 2, col = c("#e76f51", "#0081a7", "#FFBC42", "#a7c957"),
             fill = c(alpha("#e76f51", 0.8), alpha("#0081a7", 0.8), alpha("#FFBC42",0.8), alpha("#a7c957",0.8)),
             fontface = "bold", cex = 1.5, cat.cex = 1.5)

## Run and Plot PCA on DEG##
#make list of DEG from all four important comparisons 
BN.DEG.list <- c(BN.StF.vs.StP.DEG.list, BN.StF.vs.SwF.DEG.list, BN.SwF.vs.SwP.DEG.list, BN.StP.vs.SwP.DEG.list)
BN.DEG.list <- unique(BN.DEG.list) #remove duplicates

BA.DEG.list <- c(BA.StF.vs.StP.DEG.list, BA.StF.vs.SwF.DEG.list, BA.SwF.vs.SwP.DEG.list, BA.StP.vs.SwP.DEG.list)
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
                       labels = c("High-DO Source, Imm.", "High-DO Source, Acc.", 
                                  "Low-DO Source, Imm.", "Low-DO Source, Acc.")) + 
    scale_shape_manual(name = "Source and Collection",
                       labels = c("High-DO Source, Imm.", "High-DO Source, Acc.", 
                                  "Low-DO Source, Imm.", "Low-DO Source, Acc."),
                       values = c(15, 17, 15, 17)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
    stat_ellipse(aes(group = group)) + xlab(xlabel) + ylab(ylabel) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), legend.text = element_text(size=12), 
          legend.title = element_text(size = 13))
}

BN.DEG.pca12.plot <- plot.pca(t.BN.DEG.log.center, t.BN.DEG.log.center$PC1, 
                                  t.BN.DEG.log.center$PC2, "PC1 (43.4%)", "PC2 (13.2%)")
BN.DEG.pca12.plot

BN.DEG.pca23.plot <- plot.pca(t.BN.DEG.log.center, t.BN.DEG.log.center$PC2, 
                                  t.BN.DEG.log.center$PC3, "PC2 (13.2%)", "PC3 (9.3%)")
BN.DEG.pca23.plot


BA.DEG.pca12.plot <- plot.pca(t.BA.DEG.log.center, t.BA.DEG.log.center$PC1, 
                                  t.BA.DEG.log.center$PC2, "PC1 (61.5%)", "PC2 (13.1%)")
BA.DEG.pca12.plot

BA.DEG.pca23.plot <- plot.pca(t.BA.DEG.log.center, t.BA.DEG.log.center$PC2, 
                                  t.BA.DEG.log.center$PC3, "PC2 (13.1%)", "PC3 (7.6%)")
BA.DEG.pca23.plot

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
    scale_shape_manual(values = c(16, 15, 17, 4), labels = c("High-DO Source, Imm.", "High-DO Source, Acc.", 
                                                             "Low-DO Source, Imm.", "Low-DO Source, Acc.")) +
    labs(color = "Pond", shape = "Group") + geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") + stat_ellipse(aes(group = pond)) + xlab(xlabel) + ylab(ylabel) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), legend.text = element_text(size=11))
}

BN.DEG.pca12.pond.plot <- plot.pond.pca(t.BN.DEG.log.center, t.BN.DEG.log.center$PC1, t.BN.DEG.log.center$PC2, "PC1 (43.4%)", "PC2 (13.2%)")
BN.DEG.pca12.pond.plot

BN.DEG.pca23.pond.plot <- plot.pond.pca(t.BN.DEG.log.center, t.BN.DEG.log.center$PC2, t.BN.DEG.log.center$PC3, "PC2 (13.2%)", "PC3 (9.3%)")
BN.DEG.pca23.pond.plot

BA.DEG.pca12.pond.plot <- plot.pond.pca(t.BA.DEG.log.center, t.BA.DEG.log.center$PC1, t.BA.DEG.log.center$PC2, "PC1 (61.5%)", "PC2 (13.1%)")
BA.DEG.pca12.pond.plot

BA.DEG.pca23.pond.plot <- plot.pond.pca(t.BA.DEG.log.center, t.BA.DEG.log.center$PC2, t.BA.DEG.log.center$PC3, "PC2 (13.1%)", "PC3 (7.6%)")
BA.DEG.pca23.pond.plot

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
allExprBN.pca12.plot <- plot.pca(t.allExprBN.log.center, t.allExprBN.log.center$PC1, t.allExprBN.log.center$PC2, "PC1 (9%)", "PC2 (5.8%)")
allExprBN.pca12.plot

allExprBN.pca23.plot <- plot.pca(t.allExprBN.log.center, t.allExprBN.log.center$PC2, t.allExprBN.log.center$PC3, "PC2 (5.8%)", "PC3 (4.8%)")
allExprBN.pca23.plot

allExprBA.pca12.plot <- plot.pca(t.allExprBA.log.center, t.allExprBA.log.center$PC1, t.allExprBA.log.center$PC2, "PC1 (14.4%)", "PC2 (7.3%)")
allExprBA.pca12.plot

allExprBA.pca23.plot <- plot.pca(t.allExprBA.log.center, t.allExprBA.log.center$PC2, t.allExprBA.log.center$PC3, "PC2 (7.3%)", "PC3 (6.2%)")
allExprBA.pca23.plot

#plot pond information 
#add on pond information
t.allExprBN.log.center <- merge(t.allExprBN.log.center, BN.pond.data, by = "row.names", all = TRUE)
rownames(t.allExprBN.log.center) <- t.allExprBN.log.center[,1]
t.allExprBN.log.center[,1] <- NULL

t.allExprBA.log.center <- merge(t.allExprBA.log.center, BA.pond.data, by = "row.names", all = TRUE)
rownames(t.allExprBA.log.center) <- t.allExprBA.log.center[,1]
t.allExprBA.log.center[,1] <- NULL

#plot pca with pond information
allExprBN.pca12.pond.plot <- plot.pond.pca(t.allExprBN.log.center, t.allExprBN.log.center$PC1, t.allExprBN.log.center$PC2, "PC1 (9%)", "PC2 (5.8%)")
allExprBN.pca12.pond.plot

allExprBN.pca23.pond.plot <- plot.pond.pca(t.allExprBN.log.center, t.allExprBN.log.center$PC2, t.allExprBN.log.center$PC3, "PC2 (5.8%)", "PC3 (4.8%)")
allExprBN.pca23.pond.plot

allExprBA.pca12.pond.plot <- plot.pond.pca(t.allExprBA.log.center, t.allExprBA.log.center$PC1, t.allExprBA.log.center$PC2, "PC1 (14.4%)", "PC2 (7.3%)")
allExprBA.pca12.pond.plot

allExprBA.pca23.pond.plot <- plot.pond.pca(t.allExprBA.log.center, t.allExprBA.log.center$PC2, t.allExprBA.log.center$PC3, "PC2 (7.3%)", "PC3 (6.2%)")
allExprBA.pca23.pond.plot

## Run cluster analysis and plot heatmap on DEG##
#do cluster analysis on samplesand genes
BN.DEG.sample.hclust <- hclust(dist(t(BN.DEG.log.center), method = "euclidean"), method = "ward.D2")
BA.DEG.sample.hclust <- hclust(dist(t(BA.DEG.log.center), method = "euclidean"), method = "ward.D2")

BN.DEG.gene.hclust <- hclust(dist(BN.DEG.log.center, method = "euclidean"), method = "ward.D2")
BA.DEG.gene.hclust <- hclust(dist(BA.DEG.log.center, method = "euclidean"), method = "ward.D2")

#make dataframe for sample annotations 
sample.info.BN <- data.frame(Group = t.BN.DEG.log.center$group)
rownames(sample.info.BN) <- row.names(t.BN.DEG.log.center)
sample.info.BN$Group <- gsub("stream_field", "H-DO,Imm.", sample.info.BN$Group)
sample.info.BN$Group <- gsub("stream_pond", "H-DO,Acc.", sample.info.BN$Group)
sample.info.BN$Group <- gsub("swamp_field", "L-DO,Imm.", sample.info.BN$Group)
sample.info.BN$Group <- gsub("swamp_pond", "L-DO,Acc.", sample.info.BN$Group)

sample.info.BA <- data.frame(Group = t.BA.DEG.log.center$group)
rownames(sample.info.BA) <- row.names(t.BA.DEG.log.center)
sample.info.BA$Group <- gsub("stream_field", "H-DO,Imm.", sample.info.BA$Group)
sample.info.BA$Group <- gsub("stream_pond", "H-DO,Acc.", sample.info.BA$Group)
sample.info.BA$Group <- gsub("swamp_field", "L-DO,Imm.", sample.info.BA$Group)
sample.info.BA$Group <- gsub("swamp_pond", "L-DO,Acc.", sample.info.BA$Group)

#set colors for annotations
sample.colors <- list(
  Group = c("H-DO,Imm." = "#2a9d8f", "H-DO,Acc." = "#e9c46a",
            "L-DO,Imm." = "#f4a261", "L-DO,Acc." = "#e76f51"))

#generate color pallette
cor.color <- colorRampPalette(c("blue", "green", "black", "orange", "red"))(30) 

#extract TPM
BN.DEG.TPM <- extractDEG(allExprBN.TPM, BN.DEG.list)
BA.DEG.TPM <- extractDEG(allExprBA.TPM, BA.DEG.list)

#extract RPKM
BN.DEG.RPKM <- extractDEG(allExprBN.RPKM, BN.DEG.list)
BA.DEG.RPKM <- extractDEG(allExprBA.RPKM, BA.DEG.list)

#log and center
BN.DEG.TPM.log <- log.trans(BN.DEG.TPM)
BN.DEG.TPM.log.center <- med.center(BN.DEG.TPM.log)

BA.DEG.TPM.log <- log.trans(BA.DEG.TPM)
BA.DEG.TPM.log.center <- med.center(BA.DEG.TPM.log)

BN.DEG.RPKM.log <- log.trans(BN.DEG.RPKM)
BN.DEG.RPKM.log.center <- med.center(BN.DEG.RPKM.log)

BA.DEG.RPKM.log <- log.trans(BA.DEG.RPKM)
BA.DEG.RPKM.log.center <- med.center(BA.DEG.RPKM.log)

#calculate breaks 
BN.TPM.log.center.breaks <- max(abs(BN.DEG.TPM.log.center))
BA.TPM.log.center.breaks <- max(abs(BA.DEG.TPM.log.center))

BN.RPKM.log.center.breaks <- max(abs(BN.DEG.RPKM.log.center))
BA.RPKM.log.center.breaks <- max(abs(BA.DEG.RPKM.log.center))

BN.log.center.breaks <- max(abs(BN.DEG.log.center))
BA.log.center.breaks <- max(abs(BA.DEG.log.center))

#plot heatmaps
make.heatmap <- function(DEG, sample.anno, sample.hclust, gene.hclust, breaks.samp){
  pheatmap(DEG, cluster_rows = gene.hclust, cluster_cols = sample.hclust, show_colnames = FALSE, show_rownames = FALSE,
           color = cor.color, breaks = seq(-breaks.samp, breaks.samp, length.out = 30),
           annotation_col = sample.anno, annotation_colors = sample.colors,
           annotation_names_col = FALSE, annotation_names_row = FALSE, border_color = "grey")
}

BN.DEG.heatmap.TPM.log.center <- make.heatmap(BN.DEG.TPM.log.center, sample.info.BN, BN.DEG.sample.hclust, BN.DEG.gene.hclust, BN.TPM.log.center.breaks)
BN.DEG.heatmap.TPM.log.center

BN.DEG.heatmap.RPKM.log.center <- make.heatmap(BN.DEG.RPKM.log.center, sample.info.BN, BN.DEG.sample.hclust, BN.DEG.gene.hclust, BN.RPKM.log.center.breaks)
BN.DEG.heatmap.RPKM.log.center

BN.DEG.heatmap.TMM.log.center <- make.heatmap(BN.DEG.log.center, sample.info.BN, BN.DEG.sample.hclust, BN.DEG.gene.hclust, BN.log.center.breaks)
BN.DEG.heatmap.TMM.log.center


BA.DEG.heatmap.TPM.log.center <- make.heatmap(BA.DEG.TPM.log.center, sample.info.BA, BA.DEG.sample.hclust, BA.DEG.gene.hclust, BA.TPM.log.center.breaks)
BA.DEG.heatmap.TPM.log.center

BA.DEG.heatmap.RPKM.log.center <- make.heatmap(BA.DEG.RPKM.log.center, sample.info.BA, BA.DEG.sample.hclust, BA.DEG.gene.hclust, BA.RPKM.log.center.breaks)
BA.DEG.heatmap.RPKM.log.center

BA.DEG.heatmap.TMM.log.center <- make.heatmap(BA.DEG.log.center, sample.info.BA, BA.DEG.sample.hclust, BA.DEG.gene.hclust, BA.log.center.breaks)
BA.DEG.heatmap.TMM.log.center

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
  theme_bw() + scale_fill_manual(values = c("#2a9d8f", "#27187E"), labels = c("E. apluerogramma", "E. neumayeri")) + 
  xlab("Comparison Type") + ylab ("Number of Differentially Expressed Genes") + labs(fill = "Species") +
  theme(legend.position = "bottom", axis.title = element_text(size = 14), axis.text = element_text(size = 11),
        legend.text = element_text(size=11)) + scale_x_discrete(labels = function(x) str_wrap(x, width = 5))

plot.compareDEG

## Partition genes into clusters ## 
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
  pdf(filename)
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
              mfrow = c(r,c), x11 = FALSE)
  dev.off()
  }

mfuzz.plot.BA <- function(filename, mfuzz.cl, r, c){
  pdf(filename)
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
              mfrow = c(r,c), x11 = FALSE)
  dev.off()
}

mfuzz.plot.BN("BN2_mfuzz_plot.pdf", mfuzz.BN.2, 2, 1)
mfuzz.plot.BN("BN3_mfuzz_plot.pdf", mfuzz.BN.3, 2, 2)
mfuzz.plot.BN("BN4_mfuzz_plot.pdf", mfuzz.BN.4, 2, 2)
mfuzz.plot.BN("BN5_mfuzz_plot.pdf", mfuzz.BN.5, 3, 2)
mfuzz.plot.BN("BN6_mfuzz_plot.pdf", mfuzz.BN.6, 3, 2)
mfuzz.plot.BN("BN7_mfuzz_plot.pdf", mfuzz.BN.7, 3, 3)
mfuzz.plot.BN("BN8_mfuzz_plot.pdf", mfuzz.BN.8, 3, 3)
mfuzz.plot.BN("BN9_mfuzz_plot.pdf", mfuzz.BN.9, 3, 3)
mfuzz.plot.BN("BN10_mfuzz_plot.pdf", mfuzz.BN.10, 4, 3)
mfuzz.plot.BN("BN11_mfuzz_plot.pdf", mfuzz.BN.11, 4, 3)
mfuzz.plot.BN("BN12_mfuzz_plot.pdf", mfuzz.BN.12, 4, 3)
mfuzz.plot.BN("BN13_mfuzz_plot.pdf", mfuzz.BN.13, 4, 4)
mfuzz.plot.BN("BN14_mfuzz_plot.pdf", mfuzz.BN.14, 4, 4)
mfuzz.plot.BN("BN15_mfuzz_plot.pdf", mfuzz.BN.15, 4, 4)
mfuzz.plot.BN("BN16_mfuzz_plot.pdf", mfuzz.BN.16, 4, 4)
mfuzz.plot.BN("BN17_mfuzz_plot.pdf", mfuzz.BN.17, 5, 4)
mfuzz.plot.BN("BN18_mfuzz_plot.pdf", mfuzz.BN.18, 5, 4)
mfuzz.plot.BN("BN19_mfuzz_plot.pdf", mfuzz.BN.19, 5, 4)
mfuzz.plot.BN("BN20_mfuzz_plot.pdf", mfuzz.BN.20, 5, 4)

mfuzz.plot.BA("BA2_mfuzz_plot.pdf", mfuzz.BA.2, 2, 1)
mfuzz.plot.BA("BA3_mfuzz_plot.pdf", mfuzz.BA.3, 2, 2)
mfuzz.plot.BA("BA4_mfuzz_plot.pdf", mfuzz.BA.4, 2, 2)
mfuzz.plot.BA("BA5_mfuzz_plot.pdf", mfuzz.BA.5, 3, 2)
mfuzz.plot.BA("BA6_mfuzz_plot.pdf", mfuzz.BA.6, 3, 2)
mfuzz.plot.BA("BA7_mfuzz_plot.pdf", mfuzz.BA.7, 3, 3)
mfuzz.plot.BA("BA8_mfuzz_plot.pdf", mfuzz.BA.8, 3, 3)
mfuzz.plot.BA("BA9_mfuzz_plot.pdf", mfuzz.BA.9, 3, 3)
mfuzz.plot.BA("BA10_mfuzz_plot.pdf", mfuzz.BA.10, 4, 3)
mfuzz.plot.BA("BA11_mfuzz_plot.pdf", mfuzz.BA.11, 4, 3)
mfuzz.plot.BA("BA12_mfuzz_plot.pdf", mfuzz.BA.12, 4, 3)
mfuzz.plot.BA("BA13_mfuzz_plot.pdf", mfuzz.BA.13, 4, 4)
mfuzz.plot.BA("BA14_mfuzz_plot.pdf", mfuzz.BA.14, 4, 4)
mfuzz.plot.BA("BA15_mfuzz_plot.pdf", mfuzz.BA.15, 4, 4)
mfuzz.plot.BA("BA16_mfuzz_plot.pdf", mfuzz.BA.16, 4, 4)
mfuzz.plot.BA("BA17_mfuzz_plot.pdf", mfuzz.BA.17, 5, 4)
mfuzz.plot.BA("BA18_mfuzz_plot.pdf", mfuzz.BA.18, 5, 4)
mfuzz.plot.BA("BA19_mfuzz_plot.pdf", mfuzz.BA.19, 5, 4)
mfuzz.plot.BA("BA20_mfuzz_plot.pdf", mfuzz.BA.20, 5, 4)

#extract cluster info  
BN.gene.clusters <- as.data.frame(mfuzz.BN.9$cluster)
colnames(BN.gene.clusters) <- "cluster"
BN.membership <- as.data.frame(mfuzz.BN.9$membership)

BA.gene.clusters <- as.data.frame(mfuzz.BA.8$cluster)
colnames(BA.gene.clusters) <- "cluster"
BA.membership <- as.data.frame(mfuzz.BA.8$membership)

#merge datasets
BN.gene.cluster.info <- merge(BN.gene.clusters, BN.membership, by = "row.names", all = TRUE)
rownames(BN.gene.cluster.info) <- BN.gene.cluster.info[,1]
BN.gene.cluster.info[,1] <- NULL

BA.gene.cluster.info <- merge(BA.gene.clusters, BA.membership, by = "row.names", all = TRUE)
rownames(BA.gene.cluster.info) <- BA.gene.cluster.info[,1]
BA.gene.cluster.info[,1] <- NULL

#make lists of genes in each interesting cluster 
BN.cluster2 <- rownames(subset(BN.gene.cluster.info, cluster == 2))
BN.cluster9 <- rownames(subset(BN.gene.cluster.info, cluster == 9))
BA.cluster3 <- rownames(subset(BA.gene.cluster.info, cluster == 3))

#extract expression data for clusters
BN.DEG.TPM.log.center.cluster2 <- subset(BN.DEG.TPM.log.center, rownames(BN.DEG.TPM.log.center) %in% BN.cluster2)
BN.DEG.TPM.log.center.cluster9 <- subset(BN.DEG.TPM.log.center, rownames(BN.DEG.TPM.log.center) %in% BN.cluster9)
BA.DEG.TPM.log.center.cluster3 <- subset(BA.DEG.TPM.log.center, rownames(BA.DEG.TPM.log.center) %in% BA.cluster3)

#plot heatmaps of these genes
make.heatmap.geneclust <- function(DEG, sample.anno, sample.hclust, breaks.samp, clust_rows){
                            pheatmap(DEG, cluster_rows = clust_rows, cluster_cols = sample.hclust, show_colnames = FALSE, show_rownames = FALSE,
                                     color = cor.color, breaks = seq(-breaks.samp, breaks.samp, length.out = 30),
                                     annotation_col = sample.anno, annotation_colors = sample.colors,
                                     annotation_names_col = FALSE, annotation_names_row = FALSE, border_color = "NA")
                          }

BN.TPM.log.center.breaks.cluster2 <- max(abs(BN.DEG.TPM.log.center.cluster2))
BN.TPM.log.center.breaks.cluster9 <- max(abs(BN.DEG.TPM.log.center.cluster9))
BA.TPM.log.center.breaks.cluster3 <- max(abs(BA.DEG.TPM.log.center.cluster3))

BN.DEG.heatmap.TPM.log.center.cluster2 <- make.heatmap.geneclust(BN.DEG.TPM.log.center.cluster2, sample.info.BN, BN.DEG.sample.hclust, BN.TPM.log.center.breaks.cluster2, FALSE)
BN.DEG.heatmap.TPM.log.center.cluster2

BN.DEG.heatmap.TPM.log.center.cluster9 <- make.heatmap.geneclust(BN.DEG.TPM.log.center.cluster9, sample.info.BN, BN.DEG.sample.hclust, BN.TPM.log.center.breaks.cluster9, FALSE)
BN.DEG.heatmap.TPM.log.center.cluster9

BA.DEG.heatmap.TPM.log.center.cluster3 <- make.heatmap.geneclust(BA.DEG.TPM.log.center.cluster3, sample.info.BA, BA.DEG.sample.hclust, BA.TPM.log.center.breaks.cluster3, FALSE)
BA.DEG.heatmap.TPM.log.center.cluster3

## plot heatmap of genes identified as interesting in venn diagram ##
#get genes of interest 
BN.int.gene <- Reduce(intersect, list(BN.StF.vs.SwF.DEG.list,BN.SwF.vs.SwP.DEG.list))
BN.int.gene <- setdiff(BN.int.gene, BN.StP.vs.SwP.DEG.list)
BN.int.gene <- setdiff(BN.int.gene, BN.StF.vs.StP.DEG.list)

BA.int.gene <- Reduce(intersect, list(BA.StF.vs.SwF.DEG.list,BA.SwF.vs.SwP.DEG.list))
BA.int.gene <- setdiff(BA.int.gene, BA.StP.vs.SwP.DEG.list)
BA.int.gene <- setdiff(BA.int.gene, BA.StF.vs.StP.DEG.list)

#extract expression data for clusters
BN.DEG.TPM.log.center.int.gene <- subset(BN.DEG.TPM.log.center, rownames(BN.DEG.TPM.log.center) %in% BN.int.gene)
BA.DEG.TPM.log.center.int.gene <- subset(BA.DEG.TPM.log.center, rownames(BA.DEG.TPM.log.center) %in% BA.int.gene)

#plot 
BN.TPM.log.center.breaks.int.gene <- max(abs(BN.DEG.TPM.log.center.int.gene))
BA.TPM.log.center.breaks.int.gene <- max(abs(BA.DEG.TPM.log.center.int.gene))

BN.DEG.heatmap.TPM.log.center.int.gene <- make.heatmap.geneclust(BN.DEG.TPM.log.center.int.gene, sample.info.BN, BN.DEG.sample.hclust, BN.TPM.log.center.breaks.int.gene, TRUE)
BN.DEG.heatmap.TPM.log.center.int.gene

BA.DEG.heatmap.TPM.log.center.int.gene <- make.heatmap.geneclust(BA.DEG.TPM.log.center.int.gene, sample.info.BA, BA.DEG.sample.hclust, BA.TPM.log.center.breaks.int.gene, TRUE)
BA.DEG.heatmap.TPM.log.center.int.gene
