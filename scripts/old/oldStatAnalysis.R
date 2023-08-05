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
library(countToFPKM)

#read in data 
degBN <- read.table("./data/DEG/BN/salmon/diffExpr.P0.01_C2.matrix")
degBA <- read.table("./data/DEG/BA/salmon/diffExpr.P0.01_C2.matrix")

allExprBN <- read.table("./data/totalExpr/BN_bf_new_sal.gene.TMM.EXPR.matrix")
allExprBA <- read.table("./data/totalExpr/BA_bf_new_sal.gene.TMM.EXPR.matrix")

allExprBA.TPM <- read.table("./data/totalExpr/BA_bf_new_sal.gene.TPM.not_cross_norm")

BN.StF.vs.StP.DEG <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.stream_field_vs_stream_pond.edgeR.DE_results.P0.01_C2.DE.subset")
BN.StF.vs.SwF.DEG <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.stream_field_vs_swamp_field.edgeR.DE_results.P0.01_C2.DE.subset")
BN.StF.vs.SwP.DEG <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.stream_field_vs_swamp_pond.edgeR.DE_results.P0.01_C2.DE.subset")
BN.StP.vs.SwF.DEG <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.stream_pond_vs_swamp_field.edgeR.DE_results.P0.01_C2.DE.subset")
BN.StP.vs.SwP.DEG <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.stream_pond_vs_swamp_pond.edgeR.DE_results.P0.01_C2.DE.subset")
BN.SwF.vs.SwP.DEG <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.swamp_field_vs_swamp_pond.edgeR.DE_results.P0.01_C2.DE.subset")

BA.StF.vs.StP.DEG <- read.table("./data/DEG/BA/salmon/BA_bf_new_sal.gene.counts.matrix.stream_field_vs_stream_pond.edgeR.DE_results.P0.01_C2.DE.subset")
BA.StF.vs.SwF.DEG <- read.table("./data/DEG/BA/salmon/BA_bf_new_sal.gene.counts.matrix.stream_field_vs_swamp_field.edgeR.DE_results.P0.01_C2.DE.subset")
BA.StF.vs.SwP.DEG <- read.table("./data/DEG/BA/salmon/BA_bf_new_sal.gene.counts.matrix.stream_field_vs_swamp_pond.edgeR.DE_results.P0.01_C2.DE.subset")
BA.StP.vs.SwF.DEG <- read.table("./data/DEG/BA/salmon/BA_bf_new_sal.gene.counts.matrix.stream_pond_vs_swamp_field.edgeR.DE_results.P0.01_C2.DE.subset")
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
BN.StF.vs.SwP.DEG.list <- rownames(BN.StF.vs.SwP.DEG) 
BN.StP.vs.SwF.DEG.list <- rownames(BN.StP.vs.SwF.DEG)
BN.StP.vs.SwP.DEG.list <- rownames(BN.StP.vs.SwP.DEG)
BN.SwF.vs.SwP.DEG.list <- rownames(BN.SwF.vs.SwP.DEG)

BN.allDEG.list <- list(StF.vs.StCG = BN.StF.vs.StP.DEG.list, StF.vs.SwF = BN.StF.vs.SwF.DEG.list, 
                       StF.vs.SwCG = BN.StF.vs.SwP.DEG.list, StCG.vs.SwF = BN.StP.vs.SwF.DEG.list, 
                       StCG.vs.SwCG = BN.StP.vs.SwP.DEG.list, SwF.vs.SwCG = BN.SwF.vs.SwP.DEG.list)

BN.impDEG.list <- list(StF.vs.StCG = BN.StF.vs.StP.DEG.list, 
                       SwF.vs.SwCG = BN.SwF.vs.SwP.DEG.list,
                       StF.vs.SwF = BN.StF.vs.SwF.DEG.list,
                       StCG.vs.SwCG = BN.StP.vs.SwP.DEG.list)

BA.StF.vs.StP.DEG.list <- rownames(BA.StF.vs.StP.DEG)
BA.StF.vs.SwF.DEG.list <- rownames(BA.StF.vs.SwF.DEG)
BA.StF.vs.SwP.DEG.list <- rownames(BA.StF.vs.SwP.DEG)
BA.StP.vs.SwF.DEG.list <- rownames(BA.StP.vs.SwF.DEG)
BA.StP.vs.SwP.DEG.list <- rownames(BA.StP.vs.SwP.DEG)
BA.SwF.vs.SwP.DEG.list <- rownames(BA.SwF.vs.SwP.DEG)

BA.allDEG.list <- list(StF.vs.StCG = BA.StF.vs.StP.DEG.list, StF.vs.SwF = BA.StF.vs.SwF.DEG.list, 
                       StF.vs.SwCG = BA.StF.vs.SwP.DEG.list, StCG.vs.SwF = BA.StP.vs.SwF.DEG.list, 
                       StCG.vs.SwCG = BA.StP.vs.SwP.DEG.list, SwF.vs.SwCG = BA.SwF.vs.SwP.DEG.list)

BA.impDEG.list <- list(StF.vs.StCG = BA.StF.vs.StP.DEG.list, 
                       SwF.vs.SwCG = BA.SwF.vs.SwP.DEG.list,
                       StF.vs.SwF = BA.StF.vs.SwF.DEG.list,
                       StCG.vs.SwCG = BA.StP.vs.SwP.DEG.list)

#Plot upset 
BN.upset <- upset(fromList(BN.allDEG.list), nsets = 6, order.by = "freq", main.bar.color = "#27187E",
                  matrix.color = "#0081a7", sets.bar.color = "#2a9d8f", text.scale = c(1.7,1.5,1.7,1.5,1.5,1.2))
BN.upset

BN.upset.imp <- upset(fromList(BN.impDEG.list), nsets = 4, order.by = "freq", main.bar.color = "#27187E",
                  matrix.color = "#0081a7", sets.bar.color = "#2a9d8f", text.scale = c(1.7,1.5,1.7,1.5,1.5,1.2))
BN.upset.imp


BA.upset <- upset(fromList(BA.allDEG.list), nsets = 6, order.by = "freq", main.bar.color = "#27187E",
                  matrix.color = "#0081a7", sets.bar.color = "#2a9d8f", text.scale = c(1.7,1.5,1.7,1.5,1.5,1.2))
BA.upset

BA.upset.imp <- upset(fromList(BA.impDEG.list), nsets = 4, order.by = "freq", main.bar.color = "#27187E",
                      matrix.color = "#0081a7", sets.bar.color = "#2a9d8f", text.scale = c(1.7,1.5,1.7,1.5,1.5,1.2))
BA.upset.imp

#Plot venn diagram 
#function to display venn diagram 
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

display_venn(BN.impDEG.list, lwd = 2, col = c("#e76f51", "#0081a7", "#FFBC42", "#a7c957"),
             fill = c(alpha("#e76f51", 0.8), alpha("#0081a7", 0.8), alpha("#FFBC42",0.8), alpha("#a7c957",0.8)),
             fontface = "bold", cex = 1.5, cat.cex = 1.5)

display_venn(BA.impDEG.list, lwd = 2, col = c("#e76f51", "#0081a7", "#FFBC42", "#a7c957"),
             fill = c(alpha("#e76f51", 0.8), alpha("#0081a7", 0.8), alpha("#FFBC42",0.8), alpha("#a7c957",0.8)),
             fontface = "bold", cex = 1.5, cat.cex = 1.5)

## Run and Plot PCA on all DEG##
#log2 transform data 
log.trans <- function(df){
  return(log((df + 1), 2))
}

degBN.log <- log.trans(degBN)
degBA.log <- log.trans(degBA)

#median center data 
med.center <- function(df){
  rowMed <- apply(df, 1, median)
  return(df-rowMed)
}

degBN.log.center <- med.center(degBN.log)
degBA.log.center <- med.center(degBA.log)

#transpose data 
transpose.names <- function(data){
  t.data <- t(data)
  rownames(t.data) <- colnames(data)
  colnames(t.data) <- rownames(data)
  return(as.data.frame(t.data))
}

t.degBN.log.center <- transpose.names(degBN.log.center)
t.degBA.log.center <- transpose.names(degBA.log.center)

#run PCA
pca.BN <- prcomp(t.degBN.log.center)
pca.BA <- prcomp(t.degBA.log.center)

#add on site and collection data 
add.site.collection.BN <- function(data.BN) {
  data.BN <- tibble::rownames_to_column(data.BN, "group")
  data.BN$group <- str_sub(data.BN$group, end = -13)
  
  data.BN$site <- data.BN$group
  data.BN$site <- str_sub(data.BN$site, end = -6)
  data.BN$site <- gsub("_", "", data.BN$site)
  
  data.BN$collection <- data.BN$group
  data.BN$collection <- str_sub(data.BN$collection, start = 7)
  data.BN$collection <- gsub("_", "", data.BN$collection)
  
  return(data.BN)
}

add.site.collection.BA <- function(data.BA) {
  data.BA <- tibble::rownames_to_column(data.BA, "group")
  data.BA$group <- str_sub(data.BA$group, end = -10)
  data.BA$group <- gsub("d_", "d", data.BA$group)
  
  data.BA$site <- data.BA$group
  data.BA$site <- str_sub(data.BA$site, end = -6)
  data.BA$site <- gsub("_", "", data.BA$site)
  
  data.BA$collection <- data.BA$group
  data.BA$collection <- str_sub(data.BA$collection, start = 7)
  data.BA$collection <- gsub("_", "", data.BA$collection)
  
  return(data.BA)
}

t.degBN.log.center <- add.site.collection.BN(t.degBN.log.center)
t.degBA.log.center <- add.site.collection.BA(t.degBA.log.center)

#view results 
pca.BN.eig <- get_eigenvalue(pca.BN) #get eignenvalues
pca.BA.eig <- get_eigenvalue(pca.BA) 
pca.BN.var <- get_pca_var(pca.BN) #get variance
pca.BA.var <- get_pca_var(pca.BA)

summary(pca.BN)
summary(pca.BA)

#extract first 4 PCs
t.degBN.log.center <- cbind(t.degBN.log.center, pca.BN$x[,1:4])
t.degBA.log.center <- cbind(t.degBA.log.center, pca.BA$x[,1:4])

#plot PCA
plot.pca <- function(pca.data, pca.x, pca.y, xlabel, ylabel){
  ggplot(data = pca.data, aes(x = pca.x, y = pca.y, color = group, shape = collection)) + theme_bw() + geom_point() +
    scale_color_manual(values = c("#27187E", "#82DCDA", "#9d0208", "#ff5a5f"), 
                       labels = c("Stream, Field", "Stream, Common Garden", "Swamp, Field", "Swamp, Common Garden")) + 
    scale_shape_manual(values = c(15, 17), labels = c("Field", "Common Garden")) +
    labs(color = "Site/Collection", shape = "Collection") + geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") + stat_ellipse(aes(group = group)) + xlab(xlabel) + ylab(ylabel) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), legend.text = element_text(size=11))
}

BN.pca12.plot <- plot.pca(t.degBN.log.center, t.degBN.log.center$PC1, t.degBN.log.center$PC2, "PC1 (41.2%)", "PC2 (13.4%)")
BN.pca12.plot

BN.pca23.plot <- plot.pca(t.degBN.log.center, t.degBN.log.center$PC2, t.degBN.log.center$PC3, "PC2 (13.4%)", "PC3 (7.7%)")
BN.pca23.plot

BA.pca12.plot <- plot.pca(t.degBA.log.center, t.degBA.log.center$PC1, t.degBA.log.center$PC2, "PC1 (59.6%)", "PC2 (13.2%)")
BA.pca12.plot

BA.pca23.plot <- plot.pca(t.degBA.log.center, t.degBA.log.center$PC2, t.degBA.log.center$PC3, "PC2 (13.2%)", "PC3 (7.7%)")
BA.pca23.plot

#save plots 
save.pca <- function(pca.plot, name){
  ggsave(name, plot = pca.plot, path = "./plots", width = 8, height = 7)
}

save.pca(BN.pca12.plot, "BN_PCA12_allDEG.jpg")
save.pca(BN.pca23.plot, "BN_PCA23_allDEG.jpg")
save.pca(BA.pca12.plot, "BA_PCA12_allDEG.jpg")
save.pca(BA.pca23.plot, "BA_PCA23_allDEG.jpg")

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
pca.BN.eig <- get_eigenvalue(pca.allExpr.BN) #get eignenvalues
pca.BA.eig <- get_eigenvalue(pca.allExpr.BA) 
pca.BN.var <- get_pca_var(pca.allExpr.BN) #get variance
pca.BA.var <- get_pca_var(pca.allExpr.BA)

summary(pca.allExpr.BN)
summary(pca.allExpr.BA)

#extract first 4 PCs
t.allExprBN.log.center <- cbind(t.allExprBN.log.center, pca.allExpr.BN$x[,1:4])
t.allExprBA.log.center <- cbind(t.allExprBA.log.center, pca.allExpr.BA$x[,1:4])

#plot PCA
allExprBN.pca12.plot <- plot.pca(t.allExprBN.log.center, t.allExprBN.log.center$PC1, t.allExprBN.log.center$PC2, "PC1 (8.9%)", "PC2 (5.9%)")
allExprBN.pca12.plot

allExprBN.pca23.plot <- plot.pca(t.allExprBN.log.center, t.allExprBN.log.center$PC2, t.allExprBN.log.center$PC3, "PC2 (5.9%)", "PC3 (4.8%)")
allExprBN.pca23.plot

allExprBA.pca12.plot <- plot.pca(t.allExprBA.log.center, t.allExprBA.log.center$PC1, t.allExprBA.log.center$PC2, "PC1 (14.4%)", "PC2 (7.3%)")
allExprBA.pca12.plot

allExprBA.pca23.plot <- plot.pca(t.allExprBA.log.center, t.allExprBA.log.center$PC2, t.allExprBA.log.center$PC3, "PC2 (7.3%)", "PC3 (6.2%)")
allExprBA.pca23.plot

#saveplots
save.pca(allExprBN.pca12.plot, "BN_PCA12_allExp.jpg")
save.pca(allExprBN.pca23.plot, "BN_PCA23_allExp.jpg")
save.pca(allExprBA.pca12.plot, "BA_PCA12_allExp.jpg")
save.pca(allExprBA.pca23.plot, "BA_PCA23_allExp.jpg")

##Run and plot PCA for four comparisons ##
#make list of DEG from all four important comparisons 
BN.imp.DEG.list <- c(BN.StF.vs.StP.DEG.list, BN.StF.vs.SwF.DEG.list, BN.SwF.vs.SwP.DEG.list, BN.StP.vs.SwP.DEG.list)
BN.imp.DEG.list <- unique(BN.imp.DEG.list) #remove duplicates

BA.imp.DEG.list <- c(BA.StF.vs.StP.DEG.list, BA.StF.vs.SwF.DEG.list, BA.SwF.vs.SwP.DEG.list, BA.StP.vs.SwP.DEG.list)
BA.imp.DEG.list <- unique(BA.imp.DEG.list) #remove duplicates

#extract DEG
extractDEG <- function(DEG.df, DEG.list){
  return(DEG.df[DEG.list,])
}

BN.imp.DEG <- extractDEG(allExprBN, BN.imp.DEG.list)
BA.imp.DEG <- extractDEG(allExprBA, BA.imp.DEG.list)

#transform data 
BN.StF.vs.StP.DEG.log <- log.trans(BN.StF.vs.StP.DEG[,-1:-6])
BN.StF.vs.SwF.DEG.log <- log.trans(BN.StF.vs.SwF.DEG[,-1:-6])
BN.SwF.vs.SwP.DEG.log <- log.trans(BN.SwF.vs.SwP.DEG[,-1:-6])
BN.StP.vs.SwP.DEG.log <- log.trans(BN.StP.vs.SwP.DEG[,-1:-6])
BN.imp.DEG.log <- log.trans(BN.imp.DEG)

BA.StF.vs.StP.DEG.log <- log.trans(BA.StF.vs.StP.DEG[,-1:-6])
BA.StF.vs.SwF.DEG.log <- log.trans(BA.StF.vs.SwF.DEG[,-1:-6])
BA.SwF.vs.SwP.DEG.log <- log.trans(BA.SwF.vs.SwP.DEG[,-1:-6])
BA.StP.vs.SwP.DEG.log <- log.trans(BA.StP.vs.SwP.DEG[,-1:-6])
BA.imp.DEG.log <- log.trans(BA.imp.DEG)

BN.StF.vs.StP.DEG.log.center <- med.center(BN.StF.vs.StP.DEG.log)
BN.StF.vs.SwF.DEG.log.center <- med.center(BN.StF.vs.SwF.DEG.log)
BN.SwF.vs.SwP.DEG.log.center <- med.center(BN.SwF.vs.SwP.DEG.log)
BN.StP.vs.SwP.DEG.log.center <- med.center(BN.StP.vs.SwP.DEG.log)
BN.imp.DEG.log.center <- med.center(BN.imp.DEG.log)

BA.StF.vs.StP.DEG.log.center <- med.center(BA.StF.vs.StP.DEG.log)
BA.StF.vs.SwF.DEG.log.center <- med.center(BA.StF.vs.SwF.DEG.log)
BA.SwF.vs.SwP.DEG.log.center <- med.center(BA.SwF.vs.SwP.DEG.log)
BA.StP.vs.SwP.DEG.log.center <- med.center(BA.StP.vs.SwP.DEG.log)
BA.imp.DEG.log.center <- med.center(BA.imp.DEG.log)

#transpose data 
t.BN.StF.vs.StP.DEG.log.center <- transpose.names(BN.StF.vs.StP.DEG.log.center)
t.BN.StF.vs.SwF.DEG.log.center <- transpose.names(BN.StF.vs.SwF.DEG.log.center)
t.BN.SwF.vs.SwP.DEG.log.center <- transpose.names(BN.SwF.vs.SwP.DEG.log.center)
t.BN.StP.vs.SwP.DEG.log.center <- transpose.names(BN.StP.vs.SwP.DEG.log.center)
t.BN.imp.DEG.log.center <- transpose.names(BN.imp.DEG.log.center)

t.BA.StF.vs.StP.DEG.log.center <- transpose.names(BA.StF.vs.StP.DEG.log.center)
t.BA.StF.vs.SwF.DEG.log.center <- transpose.names(BA.StF.vs.SwF.DEG.log.center)
t.BA.SwF.vs.SwP.DEG.log.center <- transpose.names(BA.SwF.vs.SwP.DEG.log.center)
t.BA.StP.vs.SwP.DEG.log.center <- transpose.names(BA.StP.vs.SwP.DEG.log.center)
t.BA.imp.DEG.log.center <- transpose.names(BA.imp.DEG.log.center)

#run PCA
pca.BN.StF.vs.StP <- prcomp(t.BN.StF.vs.StP.DEG.log.center)
pca.BN.StF.vs.SwF <- prcomp(t.BN.StF.vs.SwF.DEG.log.center)
pca.BN.SwF.vs.SwP <- prcomp(t.BN.SwF.vs.SwP.DEG.log.center)
pca.BN.StP.vs.SwP <- prcomp(t.BN.StP.vs.SwP.DEG.log.center)
pca.BN.impDEG <- prcomp(t.BN.imp.DEG.log.center)

pca.BA.StF.vs.StP <- prcomp(t.BA.StF.vs.StP.DEG.log.center)
pca.BA.StF.vs.SwF <- prcomp(t.BA.StF.vs.SwF.DEG.log.center)
pca.BA.SwF.vs.SwP <- prcomp(t.BA.SwF.vs.SwP.DEG.log.center)
pca.BA.StP.vs.SwP <- prcomp(t.BA.StP.vs.SwP.DEG.log.center)
pca.BA.impDEG <- prcomp(t.BA.imp.DEG.log.center)

#add on site and collection data  
t.BN.StF.vs.StP.DEG.log.center <- add.site.collection.BN(t.BN.StF.vs.StP.DEG.log.center)
t.BN.StF.vs.SwF.DEG.log.center <- add.site.collection.BN(t.BN.StF.vs.SwF.DEG.log.center)
t.BN.SwF.vs.SwP.DEG.log.center <- add.site.collection.BN(t.BN.SwF.vs.SwP.DEG.log.center)
t.BN.StP.vs.SwP.DEG.log.center <- add.site.collection.BN(t.BN.StP.vs.SwP.DEG.log.center)
t.BN.imp.DEG.log.center <- add.site.collection.BN(t.BN.imp.DEG.log.center)

t.BA.StF.vs.StP.DEG.log.center <- add.site.collection.BA(t.BA.StF.vs.StP.DEG.log.center)
t.BA.StF.vs.SwF.DEG.log.center <- add.site.collection.BA(t.BA.StF.vs.SwF.DEG.log.center)
t.BA.SwF.vs.SwP.DEG.log.center <- add.site.collection.BA(t.BA.SwF.vs.SwP.DEG.log.center)
t.BA.StP.vs.SwP.DEG.log.center <- add.site.collection.BA(t.BA.StP.vs.SwP.DEG.log.center)
t.BA.imp.DEG.log.center <- add.site.collection.BA(t.BA.imp.DEG.log.center)

#view results 
summary(pca.BN.StF.vs.StP)
summary(pca.BN.StF.vs.SwF)
summary(pca.BN.SwF.vs.SwP)
summary(pca.BN.StP.vs.SwP)
summary(pca.BN.impDEG)

summary(pca.BA.StF.vs.StP)
summary(pca.BA.StF.vs.SwF)
summary(pca.BA.SwF.vs.SwP)
summary(pca.BA.StP.vs.SwP)
summary(pca.BA.impDEG)

#extract first 4 PCs
t.BN.StF.vs.StP.DEG.log.center <- cbind(t.BN.StF.vs.StP.DEG.log.center, pca.BN.StF.vs.StP$x[,1:4])
t.BN.StF.vs.SwF.DEG.log.center <- cbind(t.BN.StF.vs.SwF.DEG.log.center, pca.BN.StF.vs.SwF$x[,1:4])
t.BN.SwF.vs.SwP.DEG.log.center <- cbind(t.BN.SwF.vs.SwP.DEG.log.center, pca.BN.SwF.vs.SwP$x[,1:4])
t.BN.StP.vs.SwP.DEG.log.center <- cbind(t.BN.StP.vs.SwP.DEG.log.center, pca.BN.StP.vs.SwP$x[,1:4])
t.BN.imp.DEG.log.center <- cbind(t.BN.imp.DEG.log.center, pca.BN.impDEG$x[,1:4])

t.BA.StF.vs.StP.DEG.log.center <- cbind(t.BA.StF.vs.StP.DEG.log.center, pca.BA.StF.vs.StP$x[,1:4])
t.BA.StF.vs.SwF.DEG.log.center <- cbind(t.BA.StF.vs.SwF.DEG.log.center, pca.BA.StF.vs.SwF$x[,1:4])
t.BA.SwF.vs.SwP.DEG.log.center <- cbind(t.BA.SwF.vs.SwP.DEG.log.center, pca.BA.SwF.vs.SwP$x[,1:4])
t.BA.StP.vs.SwP.DEG.log.center <- cbind(t.BA.StP.vs.SwP.DEG.log.center, pca.BA.StP.vs.SwP$x[,1:4])
t.BA.imp.DEG.log.center <- cbind(t.BA.imp.DEG.log.center, pca.BA.impDEG$x[,1:4])

#plot PCA
BN.StF.vs.StP.pca12.plot <- plot.pca(t.BN.StF.vs.StP.DEG.log.center, t.BN.StF.vs.StP.DEG.log.center$PC1, 
                                     t.BN.StF.vs.StP.DEG.log.center$PC2, "PC1 (48%)", "PC2 (10.4%)")
BN.StF.vs.StP.pca12.plot

BN.StF.vs.StP.pca23.plot <- plot.pca(t.BN.StF.vs.StP.DEG.log.center, t.BN.StF.vs.StP.DEG.log.center$PC2, 
                                     t.BN.StF.vs.StP.DEG.log.center$PC3, "PC2 (10.4%)", "PC2 (7.7%)")
BN.StF.vs.StP.pca23.plot

BN.StF.vs.SwF.pca12.plot <- plot.pca(t.BN.StF.vs.SwF.DEG.log.center, t.BN.StF.vs.SwF.DEG.log.center$PC1, 
                                     t.BN.StF.vs.SwF.DEG.log.center$PC2, "PC1 (41.4%)", "PC2 (20.8%)")
BN.StF.vs.SwF.pca12.plot

BN.StF.vs.SwF.pca23.plot <- plot.pca(t.BN.StF.vs.SwF.DEG.log.center, t.BN.StF.vs.SwF.DEG.log.center$PC2, 
                                     t.BN.StF.vs.SwF.DEG.log.center$PC3, "PC2 (20.8%)", "PC3 (8.3%)")
BN.StF.vs.SwF.pca23.plot

BN.SwF.vs.SwP.pca12.plot <- plot.pca(t.BN.SwF.vs.SwP.DEG.log.center, t.BN.SwF.vs.SwP.DEG.log.center$PC1, 
                                     t.BN.SwF.vs.SwP.DEG.log.center$PC2, "PC1 (58.8%)", "PC2 (13.3%)")
BN.SwF.vs.SwP.pca12.plot

BN.SwF.vs.SwP.pca23.plot <- plot.pca(t.BN.SwF.vs.SwP.DEG.log.center, t.BN.SwF.vs.SwP.DEG.log.center$PC2, 
                                     t.BN.SwF.vs.SwP.DEG.log.center$PC3, "PC2 (13.3%)", "PC3 (4.9%)")
BN.SwF.vs.SwP.pca23.plot

BN.StP.vs.SwP.pca12.plot <- plot.pca(t.BN.StP.vs.SwP.DEG.log.center, t.BN.StP.vs.SwP.DEG.log.center$PC1, 
                                     t.BN.StP.vs.SwP.DEG.log.center$PC2, "PC1 (83.3%)", "PC2 (5.9%)")
BN.StP.vs.SwP.pca12.plot

BN.StP.vs.SwP.pca23.plot <- plot.pca(t.BN.StP.vs.SwP.DEG.log.center, t.BN.StP.vs.SwP.DEG.log.center$PC2, 
                                     t.BN.StP.vs.SwP.DEG.log.center$PC3, "PC2 (5.9%)", "PC3 (2.4%)")
BN.StP.vs.SwP.pca23.plot

BN.imp.DEG.pca12.plot <- plot.pca(t.BN.imp.DEG.log.center, t.BN.imp.DEG.log.center$PC1, 
                                  t.BN.imp.DEG.log.center$PC2, "PC1 (43.4%)", "PC2 (13.2%)")
BN.imp.DEG.pca12.plot

BN.imp.DEG.pca23.plot <- plot.pca(t.BN.imp.DEG.log.center, t.BN.imp.DEG.log.center$PC2, 
                                  t.BN.imp.DEG.log.center$PC3, "PC2 (13.2%)", "PC3 (9.3%)")
BN.imp.DEG.pca23.plot


BA.StF.vs.StP.pca12.plot <- plot.pca(t.BA.StF.vs.StP.DEG.log.center, t.BA.StF.vs.StP.DEG.log.center$PC1, 
                                     t.BA.StF.vs.StP.DEG.log.center$PC2, "PC1 (63.8%)", "PC2 (11.3%)")
BA.StF.vs.StP.pca12.plot

BA.StF.vs.StP.pca23.plot <- plot.pca(t.BA.StF.vs.StP.DEG.log.center, t.BA.StF.vs.StP.DEG.log.center$PC2, 
                                     t.BA.StF.vs.StP.DEG.log.center$PC3, "PC2 (11.3%)", "PC2 (8.9%)")
BA.StF.vs.StP.pca23.plot

BA.StF.vs.SwF.pca12.plot <- plot.pca(t.BA.StF.vs.SwF.DEG.log.center, t.BA.StF.vs.SwF.DEG.log.center$PC1, 
                                     t.BA.StF.vs.SwF.DEG.log.center$PC2, "PC1 (48.1%)", "PC2 (24.9%)")
BA.StF.vs.SwF.pca12.plot

BA.StF.vs.SwF.pca23.plot <- plot.pca(t.BA.StF.vs.SwF.DEG.log.center, t.BA.StF.vs.SwF.DEG.log.center$PC2, 
                                     t.BA.StF.vs.SwF.DEG.log.center$PC3, "PC2 (24.9%)", "PC3 (7.9%)")
BA.StF.vs.SwF.pca23.plot

BA.SwF.vs.SwP.pca12.plot <- plot.pca(t.BA.SwF.vs.SwP.DEG.log.center, t.BA.SwF.vs.SwP.DEG.log.center$PC1, 
                                     t.BA.SwF.vs.SwP.DEG.log.center$PC2, "PC1 (73.9%)", "PC2 (12.2%)")
BA.SwF.vs.SwP.pca12.plot

BA.SwF.vs.SwP.pca23.plot <- plot.pca(t.BA.SwF.vs.SwP.DEG.log.center, t.BA.SwF.vs.SwP.DEG.log.center$PC2, 
                                     t.BA.SwF.vs.SwP.DEG.log.center$PC3, "PC2 (12.2%)", "PC3 (2.9%)")
BA.SwF.vs.SwP.pca23.plot

BA.StP.vs.SwP.pca12.plot <- plot.pca(t.BA.StP.vs.SwP.DEG.log.center, t.BA.StP.vs.SwP.DEG.log.center$PC1, 
                                     t.BA.StP.vs.SwP.DEG.log.center$PC2, "PC1 (88.5%)", "PC2 (3.3%)")
BA.StP.vs.SwP.pca12.plot

BA.StP.vs.SwP.pca23.plot <- plot.pca(t.BA.StP.vs.SwP.DEG.log.center, t.BA.StP.vs.SwP.DEG.log.center$PC2, 
                                     t.BA.StP.vs.SwP.DEG.log.center$PC3, "PC2 (3.3%)", "PC3 (2.8%)")
BA.StP.vs.SwP.pca23.plot

BA.imp.DEG.pca12.plot <- plot.pca(t.BA.imp.DEG.log.center, t.BA.imp.DEG.log.center$PC1, 
                                  t.BA.imp.DEG.log.center$PC2, "PC1 (61.5%)", "PC2 (13.1%)")
BA.imp.DEG.pca12.plot

BA.imp.DEG.pca23.plot <- plot.pca(t.BA.imp.DEG.log.center, t.BA.imp.DEG.log.center$PC2, 
                                  t.BA.imp.DEG.log.center$PC3, "PC2 (13.1%)", "PC3 (7.6%)")
BA.imp.DEG.pca23.plot

#save plots
save.pca(BN.StF.vs.StP.pca12.plot, "BN_PCA12_StFvsStP.jpg")
save.pca(BN.StF.vs.StP.pca23.plot, "BN_PCA23_StFvsStP.jpg")
save.pca(BN.StF.vs.SwF.pca12.plot, "BN_PCA12_StFvsSwF.jpg")
save.pca(BN.StF.vs.SwF.pca23.plot, "BN_PCA23_StFvsSwF.jpg")
save.pca(BN.SwF.vs.SwP.pca12.plot, "BN_PCA12_SwFvsSwP.jpg")
save.pca(BN.SwF.vs.SwP.pca23.plot, "BN_PCA23_SwFvsSwP.jpg")
save.pca(BN.StP.vs.SwP.pca12.plot, "BN_PCA12_StPvsSwP.jpg")
save.pca(BN.StP.vs.SwP.pca23.plot, "BN_PCA23_StPvsSwP.jpg")
save.pca(BN.imp.DEG.pca12.plot, "BN_PCA12_impDEG.jpg")
save.pca(BN.imp.DEG.pca23.plot, "BN_PCA23_impDEG.jpg")

save.pca(BA.StF.vs.StP.pca12.plot, "BA_PCA12_StFvsStP.jpg")
save.pca(BA.StF.vs.StP.pca23.plot, "BA_PCA23_StFvsStP.jpg")
save.pca(BA.StF.vs.SwF.pca12.plot, "BA_PCA12_StFvsSwF.jpg")
save.pca(BA.StF.vs.SwF.pca23.plot, "BA_PCA23_StFvsSwF.jpg")
save.pca(BA.SwF.vs.SwP.pca12.plot, "BA_PCA12_SwFvsSwP.jpg")
save.pca(BA.SwF.vs.SwP.pca23.plot, "BA_PCA23_SwFvsSwP.jpg")
save.pca(BA.StP.vs.SwP.pca12.plot, "BA_PCA12_StPvsSwP.jpg")
save.pca(BA.StP.vs.SwP.pca23.plot, "BA_PCA23_StPvsSwP.jpg")
save.pca(BA.imp.DEG.pca12.plot, "BA_PCA12_impDEG.jpg")
save.pca(BA.imp.DEG.pca23.plot, "BA_PCA23_impDEG.jpg")

#plot pond information 
#add on pond information
t.BN.imp.DEG.log.center <- merge(t.BN.imp.DEG.log.center, BN.pond.data, by = "row.names", all = TRUE)
rownames(t.BN.imp.DEG.log.center) <- t.BN.imp.DEG.log.center[,1]
t.BN.imp.DEG.log.center[,1] <- NULL

t.BA.imp.DEG.log.center <- merge(t.BA.imp.DEG.log.center, BA.pond.data, by = "row.names", all = TRUE)
rownames(t.BA.imp.DEG.log.center) <- t.BA.imp.DEG.log.center[,1]
t.BA.imp.DEG.log.center[,1] <- NULL

#plot PCA with pond information 
plot.pond.pca <- function(pca.data, pca.x, pca.y, xlabel, ylabel){
  ggplot(data = pca.data, aes(x = pca.x, y = pca.y, color = pond, shape = group)) + theme_bw() + geom_point() +
    scale_color_manual(values = c("#9d0208", "#27187E", "#82DCDA")) + 
    scale_shape_manual(values = c(16, 15, 17, 4), labels = c("Stream, Field", "Stream, Common Garden","Swamp, Field", "Swamp, Common Garden")) +
    labs(color = "Pond", shape = "Group") + geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") + stat_ellipse(aes(group = pond)) + xlab(xlabel) + ylab(ylabel) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), legend.text = element_text(size=11))
}

allExprBN.pca12.pond.plot <- plot.pond.pca(t.allExprBN.log.center, t.allExprBN.log.center$PC1, t.allExprBN.log.center$PC2, "PC1 (8.9%)", "PC2 (5.9%)")

BN.imp.DEG.pca12.pond.plot <- plot.pond.pca(t.BN.imp.DEG.log.center, t.BN.imp.DEG.log.center$PC1, 
                                  t.BN.imp.DEG.log.center$PC2, "PC1 (43.4%)", "PC2 (13.2%)")
BN.imp.DEG.pca12.pond.plot

BN.imp.DEG.pca23.pond.plot <- plot.pond.pca(t.BN.imp.DEG.log.center, t.BN.imp.DEG.log.center$PC2, 
                                  t.BN.imp.DEG.log.center$PC3, "PC2 (13.2%)", "PC3 (9.3%)")
BN.imp.DEG.pca23.pond.plot

BA.imp.DEG.pca12.pond.plot <- plot.pond.pca(t.BA.imp.DEG.log.center, t.BA.imp.DEG.log.center$PC1, 
                                  t.BA.imp.DEG.log.center$PC2, "PC1 (61.5%)", "PC2 (13.1%)")
BA.imp.DEG.pca12.pond.plot

BA.imp.DEG.pca23.pond.plot <- plot.pond.pca(t.BA.imp.DEG.log.center, t.BA.imp.DEG.log.center$PC2, 
                                  t.BA.imp.DEG.log.center$PC3, "PC2 (13.1%)", "PC3 (7.6%)")
BA.imp.DEG.pca23.pond.plot

save.pca(BN.imp.DEG.pca12.pond.plot, "BN_PCA12_impDEG_pond.jpg")
save.pca(BN.imp.DEG.pca23.pond.plot, "BN_PCA23_impDEG_pond.jpg")
save.pca(BA.imp.DEG.pca12.pond.plot, "BA_PCA12_impDEG_pond.jpg")
save.pca(BA.imp.DEG.pca23.pond.plot, "BA_PCA23_impDEG_pond.jpg")

## Run cluster analysis and plot heatmap on all DEG##
#do cluster analysis 
BN.gene.hclust <- hclust(dist(degBN.log.center), method = "ward.D2")
BN.sample.hclust <- hclust(dist(t(degBN.log.center)), method = "ward.D2")

BA.gene.hclust <- hclust(dist(degBA.log.center), method = "ward.D2")
BA.sample.hclust <- hclust(dist(t(degBA.log.center)), method = "ward.D2")

#find gene clusters by cutting tree at 50% height 
find.gene.clust <- function(DEG.matrix, hclust.dat){
  gene.clust <- cutree(tree = as.dendrogram(hclust.dat), h = (0.4*max(hclust.dat$height))) 
  print(max(gene.clust))
  return(gene.clust)
}

BN.gene.clust.dat <- find.gene.clust(degBN.log.center, BN.gene.hclust)
BA.gene.clust.dat <- find.gene.clust(degBA.log.center, BA.gene.hclust)

#make cluster info into dataframe for annotations
gene.clust.BN <- data.frame(Cluster = BN.gene.clust.dat)
gene.clust.BN$Cluster <- gsub("1", "Cluster1", gene.clust.BN$Cluster)
gene.clust.BN$Cluster <- gsub("2", "Cluster2", gene.clust.BN$Cluster)
gene.clust.BN$Cluster <- gsub("3", "Cluster3", gene.clust.BN$Cluster)
gene.clust.BN$Cluster <- gsub("4", "Cluster4", gene.clust.BN$Cluster)
gene.clust.BN$Cluster <- as.factor(gene.clust.BN$Cluster)

gene.clust.BA <- data.frame(Cluster = BA.gene.clust.dat)
gene.clust.BA$Cluster <- gsub("1", "Cluster1", gene.clust.BA$Cluster)
gene.clust.BA$Cluster <- gsub("2", "Cluster2", gene.clust.BA$Cluster)
gene.clust.BA$Cluster <- gsub("3", "Cluster3", gene.clust.BA$Cluster)
gene.clust.BA$Cluster <- gsub("4", "Cluster4", gene.clust.BA$Cluster)
gene.clust.BA$Cluster <- gsub("5", "Cluster5", gene.clust.BA$Cluster)
gene.clust.BA$Cluster <- as.factor(gene.clust.BA$Cluster)

#make dataframe for sample annotations 
sample.info.BN <- data.frame(Group = t.degBN.log.center$group)
rownames(sample.info.BN) <- row.names(t.degBN.log.center)
sample.info.BN$Group <- gsub("stream_field", "Stream,Field", sample.info.BN$Group)
sample.info.BN$Group <- gsub("stream_pond", "Stream,CG", sample.info.BN$Group)
sample.info.BN$Group <- gsub("swamp_field", "Swamp,Field", sample.info.BN$Group)
sample.info.BN$Group <- gsub("swamp_pond", "Swamp,CG", sample.info.BN$Group)

sample.info.BA <- data.frame(Group = t.degBA.log.center$group)
rownames(sample.info.BA) <- row.names(t.degBA.log.center)
sample.info.BA$Group <- gsub("stream_field", "Stream,Field", sample.info.BA$Group)
sample.info.BA$Group <- gsub("stream_pond", "Stream,CG", sample.info.BA$Group)
sample.info.BA$Group <- gsub("swamp_field", "Swamp,Field", sample.info.BA$Group)
sample.info.BA$Group <- gsub("swamp_pond", "Swamp,CG", sample.info.BA$Group)

#set colors for annotations
four.clust.colors <- list(
  Group = c("Stream,Field" = "#2a9d8f", "Stream,CG" = "#e9c46a",
            "Swamp,Field" = "#f4a261", "Swamp,CG" = "#e76f51"),
  Cluster = c(Cluster1 = "#fb6f92", Cluster2 = "#a7c957", Cluster3 = "#FFBC42", Cluster4 = "#0077b6"))

five.clust.colors <- list(
  Group = c("Stream,Field" = "#2a9d8f", "Stream,CG" = "#e9c46a",
            "Swamp,Field" = "#f4a261", "Swamp,CG" = "#e76f51"),
  Cluster = c(Cluster1 = "#fb6f92", Cluster2 = "#a7c957", Cluster3 = "#FFBC42", Cluster4 = "#0077b6", Cluster5 = "#ff7d00"))


# convert log expression data to z scores for better plotting
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

degBN.log.zscore <- t(apply(degBN.log, 1, cal_z_score))
degBA.log.zscore <- t(apply(degBA.log, 1, cal_z_score))

#plot heatmaps
make.heatmap <- function(DEG, gene.anno, sample.anno, color.anno, gene.hclust, sample.hclust){
  pheatmap(DEG, cluster_rows = gene.hclust, cluster_cols = sample.hclust, show_colnames = FALSE, show_rownames = FALSE,
           annotation_row = gene.anno, annotation_col = sample.anno, annotation_colors = color.anno,
           annotation_names_col = FALSE, annotation_names_row = FALSE, border_color = "NA")
}

BN.heatmap <- make.heatmap(degBN.log.zscore, gene.clust.BN, sample.info.BN, four.clust.colors, BN.gene.hclust, BN.sample.hclust)
BN.heatmap

BA.heatmap <- make.heatmap(degBA.log.zscore, gene.clust.BA, sample.info.BA, five.clust.colors, BA.gene.hclust, BA.sample.hclust)
BA.heatmap

## Make heat maps for each specific comparison and one for all three important comparisons combined ## 
#do cluster analysis 
BN.StF.vs.StP.gene.hclust <- hclust(dist(BN.StF.vs.StP.DEG.log.center), method = "ward.D2")
BN.StF.vs.SwF.gene.hclust <- hclust(dist(BN.StF.vs.SwF.DEG.log.center), method = "ward.D2")
BN.SwF.vs.SwP.gene.hclust <- hclust(dist(BN.SwF.vs.SwP.DEG.log.center), method = "ward.D2")
BN.StP.vs.SwP.gene.hclust <- hclust(dist(BN.StP.vs.SwP.DEG.log.center), method = "ward.D2")
BN.imp.DEG.gene.hclust <- hclust(dist(BN.imp.DEG.log.center), method = "ward.D2")

BN.StF.vs.StP.sample.hclust <- hclust(dist(t(BN.StF.vs.StP.DEG.log.center)), method = "ward.D2")
BN.StF.vs.SwF.sample.hclust <- hclust(dist(t(BN.StF.vs.SwF.DEG.log.center)), method = "ward.D2")
BN.SwF.vs.SwP.sample.hclust <- hclust(dist(t(BN.SwF.vs.SwP.DEG.log.center)), method = "ward.D2")
BN.StP.vs.SwP.sample.hclust <- hclust(dist(t(BN.StP.vs.SwP.DEG.log.center)), method = "ward.D2")
BN.imp.DEG.sample.hclust <- hclust(dist(t(BN.imp.DEG.log.center)), method = "ward.D2")

BA.StF.vs.StP.gene.hclust <- hclust(dist(BA.StF.vs.StP.DEG.log.center), method = "ward.D2")
BA.StF.vs.SwF.gene.hclust <- hclust(dist(BA.StF.vs.SwF.DEG.log.center), method = "ward.D2")
BA.SwF.vs.SwP.gene.hclust <- hclust(dist(BA.SwF.vs.SwP.DEG.log.center), method = "ward.D2")
BA.StP.vs.SwP.gene.hclust <- hclust(dist(BA.StP.vs.SwP.DEG.log.center), method = "ward.D2")
BA.imp.DEG.gene.hclust <- hclust(dist(BA.imp.DEG.log.center), method = "ward.D2")

BA.StF.vs.StP.sample.hclust <- hclust(dist(t(BA.StF.vs.StP.DEG.log.center)), method = "ward.D2")
BA.StF.vs.SwF.sample.hclust <- hclust(dist(t(BA.StF.vs.SwF.DEG.log.center)), method = "ward.D2")
BA.SwF.vs.SwP.sample.hclust <- hclust(dist(t(BA.SwF.vs.SwP.DEG.log.center)), method = "ward.D2")
BA.StP.vs.SwP.sample.hclust <- hclust(dist(t(BA.StP.vs.SwP.DEG.log.center)), method = "ward.D2")
BA.imp.DEG.sample.hclust <- hclust(dist(t(BA.imp.DEG.log.center)), method = "ward.D2")

#cut tree into clusters based on height of tree
BN.StF.vs.StP.gene.clust <- find.gene.clust(BN.StF.vs.StP.DEG.log.center, BN.StF.vs.StP.gene.hclust)
BN.StF.vs.SwF.gene.clust <- find.gene.clust(BN.StF.vs.SwF.DEG.log.center, BN.StF.vs.SwF.gene.hclust)
BN.SwF.vs.SwP.gene.clust <- find.gene.clust(BN.SwF.vs.SwP.DEG.log.center, BN.SwF.vs.SwP.gene.hclust)
BN.StP.vs.SwP.gene.clust <- find.gene.clust(BN.StP.vs.SwP.DEG.log.center, BN.StP.vs.SwP.gene.hclust)
BN.imp.DEG.gene.clust <- find.gene.clust(BN.imp.DEG.log.center, BN.imp.DEG.gene.hclust)

BA.StF.vs.StP.gene.clust <- find.gene.clust(BA.StF.vs.StP.DEG.log.center, BA.StF.vs.StP.gene.hclust)
BA.StF.vs.SwF.gene.clust <- find.gene.clust(BA.StF.vs.SwF.DEG.log.center, BA.StF.vs.SwF.gene.hclust)
BA.SwF.vs.SwP.gene.clust <- find.gene.clust(BA.SwF.vs.SwP.DEG.log.center, BA.SwF.vs.SwP.gene.hclust)
BA.StP.vs.SwP.gene.clust <- find.gene.clust(BA.StP.vs.SwP.DEG.log.center, BA.StP.vs.SwP.gene.hclust)
BA.imp.DEG.gene.clust <- find.gene.clust(BA.imp.DEG.log.center, BA.imp.DEG.gene.hclust)

#make cluster info into dataframe for annotations
anno.gene.clust.BN.StF.vs.StP <- data.frame(Cluster = BN.StF.vs.StP.gene.clust)
anno.gene.clust.BN.StF.vs.StP$Cluster <- gsub("1", "Cluster1", anno.gene.clust.BN.StF.vs.StP$Cluster)
anno.gene.clust.BN.StF.vs.StP$Cluster <- gsub("2", "Cluster2", anno.gene.clust.BN.StF.vs.StP$Cluster)
anno.gene.clust.BN.StF.vs.StP$Cluster <- gsub("3", "Cluster3", anno.gene.clust.BN.StF.vs.StP$Cluster)
anno.gene.clust.BN.StF.vs.StP$Cluster <- gsub("4", "Cluster4", anno.gene.clust.BN.StF.vs.StP$Cluster)
anno.gene.clust.BN.StF.vs.StP$Cluster <- gsub("5", "Cluster5", anno.gene.clust.BN.StF.vs.StP$Cluster)

anno.gene.clust.BN.StF.vs.SwF <- data.frame(Cluster = BN.StF.vs.SwF.gene.clust)
anno.gene.clust.BN.StF.vs.SwF$Cluster <- gsub("1", "Cluster1", anno.gene.clust.BN.StF.vs.SwF$Cluster)
anno.gene.clust.BN.StF.vs.SwF$Cluster <- gsub("2", "Cluster2", anno.gene.clust.BN.StF.vs.SwF$Cluster)
anno.gene.clust.BN.StF.vs.SwF$Cluster <- gsub("3", "Cluster3", anno.gene.clust.BN.StF.vs.SwF$Cluster)
anno.gene.clust.BN.StF.vs.SwF$Cluster <- gsub("4", "Cluster4", anno.gene.clust.BN.StF.vs.SwF$Cluster)
anno.gene.clust.BN.StF.vs.SwF$Cluster <- gsub("5", "Cluster5", anno.gene.clust.BN.StF.vs.SwF$Cluster)
anno.gene.clust.BN.StF.vs.SwF$Cluster <- gsub("6", "Cluster6", anno.gene.clust.BN.StF.vs.SwF$Cluster)

anno.gene.clust.BN.SwF.vs.SwP <- data.frame(Cluster = BN.SwF.vs.SwP.gene.clust)
anno.gene.clust.BN.SwF.vs.SwP$Cluster <- gsub("1", "Cluster1", anno.gene.clust.BN.SwF.vs.SwP$Cluster)
anno.gene.clust.BN.SwF.vs.SwP$Cluster <- gsub("2", "Cluster2", anno.gene.clust.BN.SwF.vs.SwP$Cluster)
anno.gene.clust.BN.SwF.vs.SwP$Cluster <- gsub("3", "Cluster3", anno.gene.clust.BN.SwF.vs.SwP$Cluster)
anno.gene.clust.BN.SwF.vs.SwP$Cluster <- gsub("4", "Cluster4", anno.gene.clust.BN.SwF.vs.SwP$Cluster)
anno.gene.clust.BN.SwF.vs.SwP$Cluster <- gsub("5", "Cluster5", anno.gene.clust.BN.SwF.vs.SwP$Cluster)

anno.gene.clust.BN.StP.vs.SwP <- data.frame(Cluster = BN.StP.vs.SwP.gene.clust)
anno.gene.clust.BN.StP.vs.SwP$Cluster <- gsub("1", "Cluster1", anno.gene.clust.BN.StP.vs.SwP$Cluster)
anno.gene.clust.BN.StP.vs.SwP$Cluster <- gsub("2", "Cluster2", anno.gene.clust.BN.StP.vs.SwP$Cluster)
anno.gene.clust.BN.StP.vs.SwP$Cluster <- gsub("3", "Cluster3", anno.gene.clust.BN.StP.vs.SwP$Cluster)

anno.gene.clust.BN.imp <- data.frame(Cluster = BN.imp.DEG.gene.clust)
anno.gene.clust.BN.imp$Cluster <- gsub("1", "Cluster1", anno.gene.clust.BN.imp$Cluster)
anno.gene.clust.BN.imp$Cluster <- gsub("2", "Cluster2", anno.gene.clust.BN.imp$Cluster)
anno.gene.clust.BN.imp$Cluster <- gsub("3", "Cluster3", anno.gene.clust.BN.imp$Cluster)
anno.gene.clust.BN.imp$Cluster <- gsub("4", "Cluster4", anno.gene.clust.BN.imp$Cluster)

anno.gene.clust.BA.StF.vs.StP <- data.frame(Cluster = BA.StF.vs.StP.gene.clust)
anno.gene.clust.BA.StF.vs.StP$Cluster <- gsub("1", "Cluster1", anno.gene.clust.BA.StF.vs.StP$Cluster)
anno.gene.clust.BA.StF.vs.StP$Cluster <- gsub("2", "Cluster2", anno.gene.clust.BA.StF.vs.StP$Cluster)
anno.gene.clust.BA.StF.vs.StP$Cluster <- gsub("3", "Cluster3", anno.gene.clust.BA.StF.vs.StP$Cluster)
anno.gene.clust.BA.StF.vs.StP$Cluster <- gsub("4", "Cluster4", anno.gene.clust.BA.StF.vs.StP$Cluster)
anno.gene.clust.BA.StF.vs.StP$Cluster <- gsub("5", "Cluster5", anno.gene.clust.BA.StF.vs.StP$Cluster)
anno.gene.clust.BA.StF.vs.StP$Cluster <- gsub("6", "Cluster6", anno.gene.clust.BA.StF.vs.StP$Cluster)

anno.gene.clust.BA.StF.vs.SwF <- data.frame(Cluster = BA.StF.vs.SwF.gene.clust)
anno.gene.clust.BA.StF.vs.SwF$Cluster <- gsub("1", "Cluster1", anno.gene.clust.BA.StF.vs.SwF$Cluster)
anno.gene.clust.BA.StF.vs.SwF$Cluster <- gsub("2", "Cluster2", anno.gene.clust.BA.StF.vs.SwF$Cluster)
anno.gene.clust.BA.StF.vs.SwF$Cluster <- gsub("3", "Cluster3", anno.gene.clust.BA.StF.vs.SwF$Cluster)
anno.gene.clust.BA.StF.vs.SwF$Cluster <- gsub("4", "Cluster4", anno.gene.clust.BA.StF.vs.SwF$Cluster)
anno.gene.clust.BA.StF.vs.SwF$Cluster <- gsub("5", "Cluster5", anno.gene.clust.BA.StF.vs.SwF$Cluster)
anno.gene.clust.BA.StF.vs.SwF$Cluster <- gsub("6", "Cluster6", anno.gene.clust.BA.StF.vs.SwF$Cluster)

anno.gene.clust.BA.SwF.vs.SwP <- data.frame(Cluster = BA.SwF.vs.SwP.gene.clust)
anno.gene.clust.BA.SwF.vs.SwP$Cluster <- gsub("1", "Cluster1", anno.gene.clust.BA.SwF.vs.SwP$Cluster)
anno.gene.clust.BA.SwF.vs.SwP$Cluster <- gsub("2", "Cluster2", anno.gene.clust.BA.SwF.vs.SwP$Cluster)
anno.gene.clust.BA.SwF.vs.SwP$Cluster <- gsub("3", "Cluster3", anno.gene.clust.BA.SwF.vs.SwP$Cluster)
anno.gene.clust.BA.SwF.vs.SwP$Cluster <- gsub("4", "Cluster4", anno.gene.clust.BA.SwF.vs.SwP$Cluster)

anno.gene.clust.BA.StP.vs.SwP <- data.frame(Cluster = BA.StP.vs.SwP.gene.clust)
anno.gene.clust.BA.StP.vs.SwP$Cluster <- gsub("1", "Cluster1", anno.gene.clust.BA.StP.vs.SwP$Cluster)
anno.gene.clust.BA.StP.vs.SwP$Cluster <- gsub("2", "Cluster2", anno.gene.clust.BA.StP.vs.SwP$Cluster)
anno.gene.clust.BA.StP.vs.SwP$Cluster <- gsub("3", "Cluster3", anno.gene.clust.BA.StP.vs.SwP$Cluster)
anno.gene.clust.BA.StP.vs.SwP$Cluster <- gsub("4", "Cluster4", anno.gene.clust.BA.StP.vs.SwP$Cluster)

anno.gene.clust.BA.imp <- data.frame(Cluster = BA.imp.DEG.gene.clust)
anno.gene.clust.BA.imp$Cluster <- gsub("1", "Cluster1", anno.gene.clust.BA.imp$Cluster)
anno.gene.clust.BA.imp$Cluster <- gsub("2", "Cluster2", anno.gene.clust.BA.imp$Cluster)
anno.gene.clust.BA.imp$Cluster <- gsub("3", "Cluster3", anno.gene.clust.BA.imp$Cluster)
anno.gene.clust.BA.imp$Cluster <- gsub("4", "Cluster4", anno.gene.clust.BA.imp$Cluster)
anno.gene.clust.BA.imp$Cluster <- gsub("5", "Cluster5", anno.gene.clust.BA.imp$Cluster)
anno.gene.clust.BA.imp$Cluster <- gsub("6", "Cluster6", anno.gene.clust.BA.imp$Cluster)

#set colors for annotations
thr.clust.colors <- list(
   Group = c("Stream,Field" = "#2a9d8f", "Stream,CG" = "#e9c46a",
             "Swamp,Field" = "#f4a261", "Swamp,CG" = "#e76f51"),
   Cluster = c(Cluster1 = "#fb6f92", Cluster2 = "#a7c957", Cluster3 = "#FFBC42"))

six.clust.colors <- list(
  Group = c("Stream,Field" = "#2a9d8f", "Stream,CG" = "#e9c46a",
            "Swamp,Field" = "#f4a261", "Swamp,CG" = "#e76f51"),
  Cluster = c(Cluster1 = "#fb6f92", Cluster2 = "#a7c957", Cluster3 = "#FFBC42", Cluster4 = "#0077b6", Cluster5 = "#ff7d00", Cluster6 = "#6f2dbd"))

# convert log expression data to z scores for better plotting
BN.StF.vs.StP.DEG.log.zscore <- t(apply(BN.StF.vs.StP.DEG.log, 1, cal_z_score))
BN.StF.vs.SwF.DEG.log.zscore <- t(apply(BN.StF.vs.SwF.DEG.log, 1, cal_z_score))
BN.SwF.vs.SwP.DEG.log.zscore <- t(apply(BN.SwF.vs.SwP.DEG.log, 1, cal_z_score))
BN.StP.vs.SwP.DEG.log.zscore <- t(apply(BN.StP.vs.SwP.DEG.log, 1, cal_z_score))
BN.imp.DEG.log.zscore <- t(apply(BN.imp.DEG.log, 1, cal_z_score))

BA.StF.vs.StP.DEG.log.zscore <- t(apply(BA.StF.vs.StP.DEG.log, 1, cal_z_score))
BA.StF.vs.SwF.DEG.log.zscore <- t(apply(BA.StF.vs.SwF.DEG.log, 1, cal_z_score))
BA.SwF.vs.SwP.DEG.log.zscore <- t(apply(BA.SwF.vs.SwP.DEG.log, 1, cal_z_score))
BA.StP.vs.SwP.DEG.log.zscore <- t(apply(BA.StP.vs.SwP.DEG.log, 1, cal_z_score))
BA.imp.DEG.log.zscore <- t(apply(BA.imp.DEG.log, 1, cal_z_score))

#plot heatmaps 
BN.StF.vs.StP.heatmap <- make.heatmap(BN.StF.vs.StP.DEG.log.zscore, anno.gene.clust.BN.StF.vs.StP, sample.info.BN, five.clust.colors, 
                                      BN.StF.vs.StP.gene.hclust, BN.StF.vs.StP.sample.hclust)
BN.StF.vs.StP.heatmap 

BN.StF.vs.SwF.heatmap <- make.heatmap(BN.StF.vs.SwF.DEG.log.zscore, anno.gene.clust.BN.StF.vs.SwF, sample.info.BN, six.clust.colors,
                                      BN.StF.vs.SwF.gene.hclust, BN.StF.vs.SwF.sample.hclust)
BN.StF.vs.SwF.heatmap

BN.SwF.vs.SwP.heatmap <- make.heatmap(BN.SwF.vs.SwP.DEG.log.zscore, anno.gene.clust.BN.SwF.vs.SwP, sample.info.BN, five.clust.colors,
                                      BN.SwF.vs.SwP.gene.hclust, BN.SwF.vs.SwP.sample.hclust)
BN.SwF.vs.SwP.heatmap

BN.StP.vs.SwP.heatmap <- make.heatmap(BN.StP.vs.SwP.DEG.log.zscore, anno.gene.clust.BN.StP.vs.SwP, sample.info.BN, thr.clust.colors,
                                      BN.StP.vs.SwP.gene.hclust, BN.StP.vs.SwP.sample.hclust)
BN.StP.vs.SwP.heatmap

BN.imp.DEG.heatmap <- make.heatmap(BN.imp.DEG.log.zscore, anno.gene.clust.BN.imp, sample.info.BN, four.clust.colors,
                               BN.imp.DEG.gene.hclust, BN.imp.DEG.sample.hclust)
BN.imp.DEG.heatmap


BA.StF.vs.StP.heatmap <- make.heatmap(BA.StF.vs.StP.DEG.log.zscore, anno.gene.clust.BA.StF.vs.StP, sample.info.BA, six.clust.colors,
                                      BA.StF.vs.StP.gene.hclust, BA.StF.vs.StP.sample.hclust)
BA.StF.vs.StP.heatmap 

BA.StF.vs.SwF.heatmap <- make.heatmap(BA.StF.vs.SwF.DEG.log.zscore, anno.gene.clust.BA.StF.vs.SwF, sample.info.BA, six.clust.colors,
                                      BA.StF.vs.SwF.gene.hclust, BA.StF.vs.SwF.sample.hclust)
BA.StF.vs.SwF.heatmap

BA.SwF.vs.SwP.heatmap <- make.heatmap(BA.SwF.vs.SwP.DEG.log.zscore, anno.gene.clust.BA.SwF.vs.SwP, sample.info.BA, four.clust.colors,
                                      BA.SwF.vs.SwP.gene.hclust, BA.SwF.vs.SwP.sample.hclust)
BA.SwF.vs.SwP.heatmap

BA.StP.vs.SwP.heatmap <- make.heatmap(BA.StP.vs.SwP.DEG.log.zscore, anno.gene.clust.BA.StP.vs.SwP, sample.info.BA, four.clust.colors,
                                      BA.StP.vs.SwP.gene.hclust, BA.StP.vs.SwP.sample.hclust)
BA.StP.vs.SwP.heatmap

BA.imp.DEG.heatmap <- make.heatmap(BA.imp.DEG.log.zscore, anno.gene.clust.BA.imp, sample.info.BA, six.clust.colors,
                               BA.imp.DEG.gene.hclust, BA.imp.DEG.sample.hclust)
BA.imp.DEG.heatmap


#experiment with TPM
BA.imp.DEG.TPM <- extractDEG(allExprBA.TPM, BA.imp.DEG.list)
BA.imp.DEG.TPM.log <- log.trans(BA.imp.DEG.TPM)
BA.imp.DEG.TPM.log.center <- med.center(BA.imp.DEG.TPM.log)

BA.imp.DEG.TPM.z <- t(apply(BA.imp.DEG.TPM, 1, cal_z_score))

BA.imp.DEG.heatmap.TPM <- make.heatmap(BA.imp.DEG.TPM.z, anno.gene.clust.BA.imp, sample.info.BA, six.clust.colors,
                                   BA.imp.DEG.gene.hclust, BA.imp.DEG.sample.hclust)
BA.imp.DEG.heatmap.TPM

BA.imp.DEG.log.center <- med.center(BA.imp.DEG.log)

BA.imp.DEG.heatmap.log.center <- make.heatmap(BA.imp.DEG.log.center, anno.gene.clust.BA.imp, sample.info.BA, six.clust.colors,
                                   BA.imp.DEG.gene.hclust, BA.imp.DEG.sample.hclust)
BA.imp.DEG.heatmap.log.center

BA.breaks <- max(abs(BA.imp.DEG.log.center))

cor.color <- colorRampPalette(c("green", "black", "red"))(40) #generate color pallette

BA.imp.DEG.heatmap.log.center <- pheatmap(BA.imp.DEG.TPM.log.center, cluster_cols = TRUE, cluster_rows = TRUE, show_colnames = FALSE, show_rownames = FALSE, 
                                  color = cor.color, breaks = seq(-BA.breaks, BA.breaks, length.out = 40), border_color = "grey", 
                                  clustering_method = "ward.D2")

BA.imp.DEG.heatmap.log.center

BA.StP.vs.SwP.DEG.log.heatmap <- pheatmap(BA.StP.vs.SwP.DEG.log.center, cluster_cols = TRUE, cluster_rows = TRUE, show_colnames = FALSE, show_rownames = FALSE, 
                                          color = cor.color, breaks = seq(-BA.breaks, BA.breaks, length.out = 40), border_color = "grey", 
                                          clustering_method = "ward.D2")

BA.StP.vs.SwP.DEG.log.heatmap

#try heatmap without gene clusters
six.clust <- list(Cluster = c(Cluster1 = "#fb6f92", Cluster2 = "#a7c957", Cluster3 = "#FFBC42", Cluster4 = "#0077b6", Cluster5 = "#ff7d00", Cluster6 = "#6f2dbd"))

BA.breaks <- max(abs(BA.imp.DEG.log.center))

cor.color <- colorRampPalette(c("green", "black", "red"))(40) #generate color pallette

BA.imp.nogeneclust <-  pheatmap(BA.imp.DEG.log.center, cluster_rows = FALSE, cluster_cols = BA.imp.DEG.sample.hclust, 
                                show_colnames = FALSE, show_rownames = FALSE, annotation_col = sample.info.BA, annotation_colors = six.clust,
                                 annotation_names_col = FALSE, annotation_names_row = FALSE, border_color = "NA", 
                                color = cor.color, breaks = seq(-BA.breaks, BA.breaks, length.out = 40))

BA.imp.withgeneclust <-  pheatmap(BA.imp.DEG.log.center, cluster_rows = BA.imp.DEG.gene.hclust, cluster_cols = BA.imp.DEG.sample.hclust, 
                                  annotation_row = anno.gene.clust.BA.imp,
                                show_colnames = FALSE, show_rownames = FALSE, annotation_col = sample.info.BA, annotation_colors = six.clust.colors,
                                annotation_names_col = FALSE, annotation_names_row = FALSE, border_color = "NA", 
                                color = cor.color, breaks = seq(-BA.breaks, BA.breaks, length.out = 40))

BA.imp.cluster <-  pheatmap(BA.imp.DEG.log.center.Cluster2, cluster_rows = FALSE, cluster_cols = BA.imp.DEG.sample.hclust, 
                                show_colnames = FALSE, show_rownames = FALSE, annotation_col = sample.info.BA, annotation_colors = six.clust,
                                annotation_names_col = FALSE, annotation_names_row = FALSE, border_color = "NA", 
                                color = cor.color, breaks = seq(-BA.breaks, BA.breaks, length.out = 40))

anno.gene.clust.BA.imp.Cluster2 <- subset(anno.gene.clust.BA.imp, Cluster == "Cluster2")
Cluster2.rownames <- rownames(anno.gene.clust.BA.imp.Cluster2)
BA.imp.DEG.log.center.Cluster2 <- subset(BA.imp.DEG.log.center, rownames(BA.imp.DEG.log.center) %in% Cluster2.rownames)

#try cutting tree lower 
find.gene.clust <- function(DEG.matrix, hclust.dat){
  gene.clust <- cutree(tree = as.dendrogram(hclust.dat), h = (0.3*max(hclust.dat$height))) 
  print(max(gene.clust))
  return(gene.clust)
}

BA.imp.DEG.gene.clust <- find.gene.clust(BA.imp.DEG.log.center, BA.imp.DEG.gene.hclust)

## Plot Bar Graphs comparing numbers of DEGS ## 
#construct dataframe 
comparison <- c("Str,F vs Str,CG", "Str,F vs Swp,F", "Str,F vs Swp,CG", 
                "Str,CG vs Swp,F", "Str,CG vs Swp,CG", "Swp,F vs Swp,CG",
                "Str,F vs Str,CG", "Str,F vs Swp,F", "Str,F vs Swp,CG", 
                "Str,CG vs Swp,F", "Str,CG vs Swp,CG", "Swp,F vs Swp,CG")

comparison.imp <- c("Str,F vs Str,CG", "Str,F vs Swp,F", "Swp,F vs Swp,CG", "Str,CG vs Swp,CG",
                "Str,F vs Str,CG", "Str,F vs Swp,F", "Swp,F vs Swp,CG", "Str,CG vs Swp,CG")

species <- c("EN", "EN", "EN", "EN", "EN", "EN", "EA", "EA", "EA", "EA", "EA", "EA")

species.imp <- c("EN", "EN", "EN", "EN", "EA", "EA", "EA", "EA")

numDEG <- c(nrow(BN.StF.vs.StP.DEG), nrow(BN.StF.vs.SwF.DEG), nrow(BN.StF.vs.SwP.DEG), 
            nrow(BN.StP.vs.SwF.DEG), nrow(BN.StP.vs.SwP.DEG), nrow(BN.SwF.vs.SwP.DEG), 
            nrow(BA.StF.vs.StP.DEG), nrow(BA.StF.vs.SwF.DEG), nrow(BA.StF.vs.SwP.DEG), 
            nrow(BA.StP.vs.SwF.DEG), nrow(BA.StP.vs.SwP.DEG), nrow(BA.SwF.vs.SwP.DEG))

numDEG.imp <- c(nrow(BN.StF.vs.StP.DEG), nrow(BN.StF.vs.SwF.DEG), nrow(BN.SwF.vs.SwP.DEG), nrow(BN.StP.vs.SwP.DEG), 
            nrow(BA.StF.vs.StP.DEG), nrow(BA.StF.vs.SwF.DEG), nrow(BA.SwF.vs.SwP.DEG), nrow(BA.StP.vs.SwP.DEG))

compareDEG <- data.frame(numDEG, species, comparison)

compareDEG.imp <- data.frame(numDEG.imp, species.imp, comparison.imp)

#order factor in order I want in barplot
compareDEG$comparison <- factor(compareDEG$comparison, levels = c("Str,F vs Str,CG", "Swp,F vs Swp,CG", 
                                                                  "Str,F vs Swp,F", 
                                                                  "Str,CG vs Swp,CG", "Str,F vs Swp,CG", 
                                                                  "Str,CG vs Swp,F"))

compareDEG.imp$comparison.imp <- factor(compareDEG.imp$comparison.imp, levels = c("Str,F vs Str,CG", "Swp,F vs Swp,CG", 
                                                                  "Str,F vs Swp,F", "Str,CG vs Swp,CG"))

plot.compare.allDEG <- ggplot(data = compareDEG, aes(x = comparison, y = numDEG, fill = species)) + geom_bar(stat = "identity", position = position_dodge()) + 
                    geom_text(aes(label=numDEG), vjust=1.6, color="white",position = position_dodge(0.9), size=3.5) +
                    theme_bw() + scale_fill_manual(values = c("#2a9d8f", "#27187E"), labels = c("E. apluerogramma", "E. neumayeri")) + 
                    xlab("Comparison Type") + ylab ("Number of Differentially Expressed Genes") + labs(fill = "Species") +
                    theme(legend.position = "bottom", axis.title = element_text(size = 14), axis.text = element_text(size = 11),
                          legend.text = element_text(size=11)) + scale_x_discrete(labels = function(x) str_wrap(x, width = 5))
                    
plot.compare.allDEG  

plot.compare.impDEG <- ggplot(data = compareDEG.imp, aes(x = comparison.imp, y = numDEG.imp, fill = species.imp)) + geom_bar(stat = "identity", position = position_dodge()) + 
  geom_text(aes(label=numDEG.imp), vjust=1.6, color="white",position = position_dodge(0.9), size=3.5)+
  theme_bw() + scale_fill_manual(values = c("#2a9d8f", "#27187E"), labels = c("E. apluerogramma", "E. neumayeri")) + 
  xlab("Comparison Type") + ylab ("Number of Differentially Expressed Genes") + labs(fill = "Species") +
  theme(legend.position = "bottom", axis.title = element_text(size = 14), axis.text = element_text(size = 11),
        legend.text = element_text(size=11)) + scale_x_discrete(labels = function(x) str_wrap(x, width = 5))

plot.compare.impDEG 

## Make correlation heatmap for allDEG ##
#calculate Pearson correlation matrix
cor.impDEG.BN <- cor(BN.imp.DEG.log.center)
cor.impDEG.BA <- cor(BA.imp.DEG.log.center)

#make heatmap
BN.breaks <- max(abs(cor.impDEG.BN)) 
BA.breaks <- max(abs(cor.impDEG.BA))

cor.color <- colorRampPalette(c("#003f88", "white", "#bf0603"))(100) #generate color pallette

group.color = list(Group = c("Stream,Field" = "#2a9d8f", "Stream,CG" = "#e9c46a", "Swamp,Field" = "#f4a261", 
                             "Swamp,CG" = "#e76f51"))

BN.cor.allDEG.heatmap <- pheatmap(cor.impDEG.BN, cluster_cols = TRUE, cluster_rows = TRUE, show_colnames = FALSE, show_rownames = FALSE, 
                           color = cor.color, breaks = seq(-BN.breaks, BN.breaks, length.out = 100), border_color = "NA", 
                           clustering_method = "ward.D2", annotation_col = sample.info.BN, annotation_row = sample.info.BN, 
                           annotation_colors = group.color, annotation_names_col = FALSE, annotation_names_row = FALSE)

BN.cor.allDEG.heatmap

BA.cor.allDEG.heatmap <- pheatmap(cor.impDEG.BA, cluster_cols = TRUE, cluster_rows = TRUE, show_colnames = FALSE, show_rownames = FALSE, 
                                  color = cor.color, breaks = seq(-BA.breaks, BA.breaks, length.out = 100), border_color = "NA", 
                                  clustering_method = "ward.D2", annotation_col = sample.info.BA, annotation_row = sample.info.BA, 
                                  annotation_colors = group.color, annotation_names_col = FALSE, annotation_names_row = FALSE)

BA.cor.allDEG.heatmap

##Make parallel coordinate plots ##
#making plots for the three clusters found in the 4 main DEG comparisons and for each comparison
#add cluster information on to dataframe
BN.imp.DEG.log.center.cluster <- merge(BN.imp.DEG.log.center, anno.gene.clust.BN.imp, by = "row.names", all = TRUE)
BN.StF.vs.StP.DEG.log.center.cluster <- merge(BN.StF.vs.StP.DEG.log.center, anno.gene.clust.BN.StF.vs.StP, by = "row.names", all = TRUE)
BN.StF.vs.SwF.DEG.log.center.cluster <- merge(BN.StF.vs.SwF.DEG.log.center, anno.gene.clust.BN.StF.vs.SwF, by = "row.names", all = TRUE)
BN.SwF.vs.SwP.DEG.log.center.cluster <- merge(BN.SwF.vs.SwP.DEG.log.center, anno.gene.clust.BN.SwF.vs.SwP, by = "row.names", all = TRUE)
BN.StP.vs.SwP.DEG.log.center.cluster <- merge(BN.StP.vs.SwP.DEG.log.center, anno.gene.clust.BN.StP.vs.SwP, by = "row.names", all = TRUE)

BA.imp.DEG.log.center.cluster <- merge(BA.imp.DEG.log.center, anno.gene.clust.BA.imp, by = "row.names", all = TRUE)
BA.StF.vs.StP.DEG.log.centercluster <- merge(BA.StF.vs.StP.DEG.log.center, anno.gene.clust.BA.StF.vs.StP, by = "row.names", all = TRUE)
BA.StF.vs.SwF.DEG.log.center.cluster <- merge(BA.StF.vs.SwF.DEG.log.center, anno.gene.clust.BA.StF.vs.SwF, by = "row.names", all = TRUE)
BA.SwF.vs.SwP.DEG.log.center.cluster <- merge(BA.SwF.vs.SwP.DEG.log.center, anno.gene.clust.BA.SwF.vs.SwP, by = "row.names", all = TRUE)
BA.StP.vs.SwP.DEG.log.center.cluster <- merge(BA.StP.vs.SwP.DEG.log.center, anno.gene.clust.BA.StP.vs.SwP, by = "row.names", all = TRUE)

#rename columns 
change.colnames.BN <- function(df){
  colnames(df) <- c("row", "StCG4", "StCG3", "SwF3", "SwCG7", "SwF4", "StF4", "StF3", "SwCG6",
                   "SwF5", "SwF2", "SwCG1", "StCG2", "StCG5", "StF2", "StF5", "StF8", "StF6",
                   "StF1", "SwF1", "SwCG2", "SwCG5", "StCG8", "SwF6", "StCG6", "SwF8", "StCG1",
                   "StF7", "StF9", "SwF9", "StCG7", "SwCG4", "SwCG3", "Cluster")
  return(df)

  }

change.colnames.BA <- function(df){
  colnames(df) <- c("row", "SwCG1", "SwCG6", "SwCG8", "SwCG11", "StCG2", "StCG5", "StF9", "SwF8", 
                    "SwF10", "SwF1", "StF7", "SwF6", "SwCG9", "SwCG7", "StF6", "SwF7", "StF1",
                    "StF8", "SwF9", "StCG4", "StCG3", "SwCG10", "StF2", "SwF3", "StF5", "SwF4",
                    "StCG7", "SwCG3", "SwCG4", "StCG6", "StCG1", "StF4", "SwF5", "StF3", "SwF2",
                    "SwCG5", "SwCG2", "Cluster")
  return(df)
  
}

BN.imp.DEG.log.center.cluster <- change.colnames.BN(BN.imp.DEG.log.center.cluster)
BN.StF.vs.StP.DEG.log.center.cluster <- change.colnames.BN(BN.StF.vs.StP.DEG.log.center.cluster)
BN.StF.vs.SwF.DEG.log.center.cluster <- change.colnames.BN(BN.StF.vs.SwF.DEG.log.center.cluster)
BN.SwF.vs.SwP.DEG.log.center.cluster <- change.colnames.BN(BN.SwF.vs.SwP.DEG.log.center.cluster)
BN.StP.vs.SwP.DEG.log.center.cluster <- change.colnames.BN(BN.StP.vs.SwP.DEG.log.center.cluster)

BA.imp.DEG.log.center.cluster <- change.colnames.BA(BA.imp.DEG.log.center.cluster)
BA.StF.vs.StP.DEG.log.centercluster <- change.colnames.BA(BA.StF.vs.StP.DEG.log.centercluster)
BA.StF.vs.SwF.DEG.log.center.cluster <- change.colnames.BA(BA.StF.vs.SwF.DEG.log.center.cluster)
BA.SwF.vs.SwP.DEG.log.center.cluster <- change.colnames.BA(BA.SwF.vs.SwP.DEG.log.center.cluster)
BA.StP.vs.SwP.DEG.log.center.cluster <- change.colnames.BA(BA.StP.vs.SwP.DEG.log.center.cluster)

#reorder columns into groupings 
reorder.col.BN <- function(df){
  df <- df[,c("row", "Cluster", "SwCG1", "SwCG2", "SwCG3", "SwCG4", "SwCG5", "SwCG6", "SwCG7", "SwF1", 
              "SwF2", "SwF3", "SwF4", "SwF5", "SwF6", "SwF8", "SwF9", "StCG1", "StCG2",
              "StCG3", "StCG4", "StCG5", "StCG6", "StCG7", "StCG8", "StF1", "StF2",
              "StF3", "StF4", "StF5", "StF6", "StF7", "StF8", "StF9")]
  return(df)
}

reorder.col.BA <- function(df){
  df <- df[,c("row", "Cluster", "SwCG1", "SwCG2", "SwCG3", "SwCG4", "SwCG5", "SwCG6", "SwCG7", "SwCG8",
              "SwCG9", "SwCG10", "SwCG11", "SwF1", "SwF2", "SwF3", "SwF4", "SwF5", 
              "SwF6", "SwF7", "SwF8", "SwF9", "SwF10", "StCG1", "StCG2",
              "StCG3", "StCG4", "StCG5", "StCG6", "StCG7", "StF1", "StF2",
              "StF3", "StF4", "StF5", "StF6", "StF7", "StF8", "StF9")]
  return(df)
}

BN.imp.DEG.log.center.cluster <- reorder.col.BN(BN.imp.DEG.log.center.cluster)
BN.StF.vs.StP.DEG.log.center.cluster <- reorder.col.BN(BN.StF.vs.StP.DEG.log.center.cluster)
BN.StF.vs.SwF.DEG.log.center.cluster <- reorder.col.BN(BN.StF.vs.SwF.DEG.log.center.cluster)
BN.SwF.vs.SwP.DEG.log.center.cluster <- reorder.col.BN(BN.SwF.vs.SwP.DEG.log.center.cluster)
BN.StP.vs.SwP.DEG.log.center.cluster <- reorder.col.BN(BN.StP.vs.SwP.DEG.log.center.cluster)

BA.imp.DEG.log.center.cluster <- reorder.col.BA(BA.imp.DEG.log.center.cluster)
BA.StF.vs.StP.DEG.log.center.cluster <- reorder.col.BA(BA.StF.vs.StP.DEG.log.centercluster)
BA.StF.vs.SwF.DEG.log.center.cluster <- reorder.col.BA(BA.StF.vs.SwF.DEG.log.center.cluster)
BA.SwF.vs.SwP.DEG.log.center.cluster <- reorder.col.BA(BA.SwF.vs.SwP.DEG.log.center.cluster)
BA.StP.vs.SwP.DEG.log.center.cluster <- reorder.col.BA(BA.StP.vs.SwP.DEG.log.center.cluster)

#plot 
plot.parcoord.BN <- function(clust.dat, col.values){
  ggparcoord(clust.dat, columns = 3:34, showPoints = TRUE, groupColumn = 2, alphaLines = 0.3) + 
    theme_bw() + theme(legend.position = "bottom", axis.title = element_text(size = 14), axis.text = element_text(size = 9), legend.text = element_text(size = 11)) + 
    xlab("Sample") + ylab("Log2 Median-Centered Gene Expression") + scale_color_manual(values = col.values)
}

plot.parcoord.BA <- function(clust.dat, col.values){
  ggparcoord(clust.dat, columns = 3:39, showPoints = TRUE, groupColumn = 2, alphaLines = 0.3) + 
    theme_bw() + theme(legend.position = "bottom", axis.title = element_text(size = 14), axis.text = element_text(size = 9), legend.text = element_text(size = 11)) + 
    xlab("Sample") + ylab("Log2 Median-Centered Gene Expression") + scale_color_manual(values = col.values)
}

#make color values for different numbers of clusters 
thr.clust <- c("#fb6f92", "#a7c957", "#fca311")
four.clust <- c("#fb6f92", "#a7c957", "#ff9f1c", "#90e0ef")
five.clust <- c("#fb6f92", "#a7c957", "#fca311", "#90e0ef", "#9b5de5")
six.clust <- c("#fb6f92", "#a7c957", "#fca311", "#90e0ef", "#9b5de5", "#d62828")

BN.imp.parcord.plot <- plot.parcoord.BN(BN.imp.DEG.log.center.cluster, four.clust)
BN.imp.parcord.plot

BN.StF.vs.StP.parcord.plot <- plot.parcoord.BN(BN.StF.vs.StP.DEG.log.center.cluster, five.clust)
BN.StF.vs.StP.parcord.plot

BN.StF.vs.SwF.parcord.plot <- plot.parcoord.BN(BN.StF.vs.SwF.DEG.log.center.cluster, six.clust)
BN.StF.vs.SwF.parcord.plot

BN.SwF.vs.SwP.parcord.plot <- plot.parcoord.BN(BN.SwF.vs.SwP.DEG.log.center.cluster, five.clust)
BN.SwF.vs.SwP.parcord.plot

BN.StP.vs.SwP.parcord.plot <- plot.parcoord.BN(BN.StP.vs.SwP.DEG.log.center.cluster, thr.clust)
BN.StP.vs.SwP.parcord.plot

BA.imp.parcord.plot <- plot.parcoord.BA(BA.imp.DEG.log.center.cluster, five.clust)
BA.imp.parcord.plot

BA.StF.vs.StP.parcord.plot <- plot.parcoord.BA(BA.StF.vs.StP.DEG.log.center.cluster, six.clust)
BA.StF.vs.StP.parcord.plot

BA.StF.vs.SwF.parcord.plot <- plot.parcoord.BA(BA.StF.vs.SwF.DEG.log.center.cluster, six.clust)
BA.StF.vs.SwF.parcord.plot

BA.SwF.vs.SwP.parcord.plot <- plot.parcoord.BA(BA.SwF.vs.SwP.DEG.log.center.cluster, four.clust)
BA.SwF.vs.SwP.parcord.plot

BA.StP.vs.SwP.parcord.plot <- plot.parcoord.BA(BA.StP.vs.SwP.DEG.log.center.cluster, four.clust)
BA.StP.vs.SwP.parcord.plot
