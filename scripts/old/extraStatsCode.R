#extra code for trying to plot the mfuzz clusters differently 
#extract cluster info  
BN.gene.clusters <- as.data.frame(mfuzz.BN.7$cluster)
BN.membership <- as.data.frame(mfuzz.BN.7$membership)

#merge datasets
BN.gene.cluster.info <- merge(BN.gene.clusters, BN.membership, by = "row.names", all = TRUE)
rownames(BN.gene.cluster.info) <- BN.gene.cluster.info[,1]
BN.gene.cluster.info[,1] <- NULL

BA.gene.clusters <- as.data.frame(mfuzz.BA.6$cluster)
BA.membership <- as.data.frame(mfuzz.BA.6$membership)

BA.gene.cluster.info <- merge(BA.gene.clusters, BA.membership, by = "row.names", all = TRUE)
rownames(BA.gene.cluster.info) <- BA.gene.cluster.info[,1]
BA.gene.cluster.info[,1] <- NULL

#extract expression data 
BN.std.expr <- as.data.frame(exprs(BN.DEG.log.eset.std))
BA.std.expr <- as.data.frame(exprs(BA.DEG.log.eset.std))

#merge with expression data 
BN.DEG.std.geneclust <- merge(BN.gene.cluster.info, BN.std.expr, by = "row.names", all = TRUE)
rownames(BN.DEG.std.geneclust) <- BN.DEG.std.geneclust[,1]
BN.DEG.std.geneclust[,1] <- NULL

BA.DEG.std.geneclust <- merge(BA.gene.cluster.info, BA.std.expr, by = "row.names", all = TRUE)
rownames(BA.DEG.std.geneclust) <- BA.DEG.std.geneclust[,1]
BA.DEG.std.geneclust[,1] <- NULL

#change column names
# colnames(BN.DEG.std.geneclust) <- c("topCluster", "Cluster1", "Cluster2",
#                                            "Cluster3", "Cluster4", "Cluster5",
#                                            "Cluster6", "Cluster7",
#                                            "HDO.I", "HDO.I", "HDO.I", "HDO.I",
#                                            "HDO.I", "HDO.I", "HDO.I", "HDO.I",
#                                            "HDO.I","HDO.A", "HDO.A", "HDO.A",
#                                            "HDO.A", "HDO.A", "HDO.A", "HDO.A",
#                                            "HDO.A", "LDO.A", "LDO.A", "LDO.A",
#                                            "LDO.A", "LDO.A", "LDO.A", "LDO.A",
#                                            "LDO.I", "LDO.I", "LDO.I", "LDO.I",
#                                            "LDO.I", "LDO.I", "LDO.I", "LDO.I")
# 
# colnames(BA.DEG.std.geneclust) <- c("topCluster", "Cluster1", "Cluster2",
#                                            "Cluster3", "Cluster4", "Cluster5",
#                                            "Cluster6", 
#                                            "HDO.I", "HDO.I", "HDO.I", "HDO.I",
#                                            "HDO.I", "HDO.I", "HDO.I", "HDO.I",
#                                            "HDO.I","HDO.A", "HDO.A", "HDO.A",
#                                            "HDO.A", "HDO.A", "HDO.A", "HDO.A",
#                                            "LDO.A", "LDO.A", "LDO.A", "LDO.A",
#                                            "LDO.A", "LDO.A", "LDO.A", "LDO.A",
#                                            "LDO.A", "LDO.A", "LDO.A", "LDO.I",
#                                            "LDO.I", "LDO.I", "LDO.I", "LDO.I",
#                                            "LDO.I", "LDO.I", "LDO.I", "LDO.I",
#                                            "LDO.I")

#separate into clusters with minimum membership value of 0.7
BN.cluster1 <- subset(BN.DEG.std.geneclust, Cluster1 > 0.7)
BN.cluster2 <- subset(BN.DEG.std.geneclust, Cluster2 > 0.7)
BN.cluster3 <- subset(BN.DEG.std.geneclust, Cluster3 > 0.7)
BN.cluster4 <- subset(BN.DEG.std.geneclust, Cluster4 > 0.7)
BN.cluster5 <- subset(BN.DEG.std.geneclust, Cluster5 > 0.7)
BN.cluster6 <- subset(BN.DEG.std.geneclust, Cluster6 > 0.7)
BN.cluster7 <- subset(BN.DEG.std.geneclust, Cluster7 > 0.7)

BA.cluster1 <- subset(BA.DEG.std.geneclust, Cluster1 > 0.7)
BA.cluster2 <- subset(BA.DEG.std.geneclust, Cluster2 > 0.7)
BA.cluster3 <- subset(BA.DEG.std.geneclust, Cluster3 > 0.7)
BA.cluster4 <- subset(BA.DEG.std.geneclust, Cluster4 > 0.7)
BA.cluster5 <- subset(BA.DEG.std.geneclust, Cluster5 > 0.7)
BA.cluster6 <- subset(BA.DEG.std.geneclust, Cluster6 > 0.7)

#create column for coloring by based on membership value 
color.col <- function(df, cluster){
  df <- mutate(df, color = case_when(
    cluster <= 0.8 ~ "yellow",
    cluster > 0.8 & cluster <= 0.9 ~ "orange",
    cluster > 0.9 ~ "red"
  ))
  return(df)
}

BN.cluster1 <- color.col(BN.cluster1, BN.cluster1$Cluster1)
BN.cluster2 <- color.col(BN.cluster2, BN.cluster2$Cluster2)
BN.cluster3 <- color.col(BN.cluster3, BN.cluster3$Cluster3)
BN.cluster4 <- color.col(BN.cluster4, BN.cluster4$Cluster4)
BN.cluster5 <- color.col(BN.cluster5, BN.cluster5$Cluster5)
BN.cluster6 <- color.col(BN.cluster6, BN.cluster6$Cluster6)
BN.cluster7 <- color.col(BN.cluster7, BN.cluster7$Cluster7)

BA.cluster1 <- color.col(BA.cluster1, BA.cluster1$Cluster1)
BA.cluster2 <- color.col(BA.cluster2, BA.cluster2$Cluster2)
BA.cluster3 <- color.col(BA.cluster3, BA.cluster3$Cluster3)
BA.cluster4 <- color.col(BA.cluster4, BA.cluster4$Cluster4)
BA.cluster5 <- color.col(BA.cluster5, BA.cluster5$Cluster5)
BA.cluster6 <- color.col(BA.cluster6, BA.cluster6$Cluster6)

#plot
plot.parcoord.BN <- function(df){
  ggparcoord(df, columns = 9:40, showPoints = FALSE, groupColumn = 41, alphaLines = 0.3) + 
    theme_bw() + theme(legend.position = "bottom", axis.title = element_text(size = 14), axis.text = element_text(size = 9), legend.text = element_text(size = 11)) + 
    xlab("Sample/Collection") + ylab("Expression Changes") + scale_color_manual(values = c("#f77f00", "#d62828", "#fcbf49"))
}

plot.parcoord.BN(BN.cluster1)
plot.parcoord.BN(BN.cluster2)

#try method of cutting tree from hierarchical clustering 
#cut trees at 30% height 
find.gene.clust <- function(hclust.dat){
  gene.clust <- cutree(tree = as.dendrogram(hclust.dat), h = (0.3*max(hclust.dat$height))) 
  print(max(gene.clust))
  return(gene.clust)
}

BN.geneclusters <- find.gene.clust(BN.DEG.gene.hclust)
BA.geneclusters <- find.gene.clust(BA.DEG.gene.hclust)

#add cluster information on to dataframe
BN.DEG.TPM.log.center.cluster <- merge(BN.DEG.TPM.log.center, BN.geneclusters, by = "row.names", all = TRUE)
BA.DEG.TPM.log.center.cluster <- merge(BA.DEG.TPM.log.center, BA.geneclusters, by = "row.names", all = TRUE)

#rename columns 
colnames(BN.DEG.TPM.log.center.cluster) <- c("row", "HDO.A4", "HDO.A3", "LDO.I3", "LDO.A7", "LDO.I4", "HDO.I4", "HDO.I3", "LDO.A6",
                                             "LDO.I5", "LDO.I2", "LDO.A1", "HDO.A2", "HDO.A5", "HDO.I2", "HDO.I5", "HDO.I8", "HDO.I6",
                                             "HDO.I1", "LDO.I1", "LDO.A2", "LDO.A5", "HDO.A8", "LDO.I6", "HDO.A6", "LDO.I8", "HDO.A1",
                                             "HDO.I7", "HDO.I9", "LDO.I9", "HDO.A7", "LDO.A4", "LDO.A3", "Cluster")



colnames(BA.DEG.TPM.log.center.cluster) <- c("row", "LDO.A1", "LDO.A6", "LDO.A8", "LDO.A11", "HDO.A2", "HDO.A5", "HDO.I9", "LDO.I8", 
                                             "LDO.I10", "LDO.I1", "HDO.I7", "LDO.I6", "LDO.A9", "LDO.A7", "HDO.I6", "LDO.I7", "HDO.I1",
                                             "HDO.I8", "LDO.I9", "HDO.A4", "HDO.A3", "LDO.A10", "HDO.I2", "LDO.I3", "HDO.I5", "LDO.I4",
                                             "HDO.A7", "LDO.A3", "LDO.A4", "HDO.A6", "HDO.A1", "HDO.I4", "LDO.I5", "HDO.I3", "LDO.I2",
                                             "LDO.A5", "LDO.A2", "Cluster")

BN.DEG.TPM.log.center.cluster <- BN.DEG.TPM.log.center.cluster[,c("row", "Cluster", "HDO.I1", "HDO.I2", "HDO.I3","HDO.I4","HDO.I5",  
                                                                  "HDO.I6", "HDO.I7", "HDO.I8", "HDO.I9", "HDO.A1", "HDO.A2", 
                                                                  "HDO.A3", "HDO.A4", "HDO.A5", "HDO.A6", "HDO.A7", "HDO.A8",
                                                                  "LDO.A1", "LDO.A2", "LDO.A3", "LDO.A4", "LDO.A5", "LDO.A6", "LDO.A7", 
                                                                  "LDO.I1", "LDO.I2", "LDO.I3", "LDO.I4", "LDO.I5", "LDO.I6",  "LDO.I8", 
                                                                  "LDO.I9")]

BA.DEG.TPM.log.center.cluster <- BA.DEG.TPM.log.center.cluster[,c("row", "Cluster", "HDO.I1", "HDO.I2", "HDO.I3", "HDO.I4", "HDO.I5", 
                                                                  "HDO.I6", "HDO.I7", "HDO.I8", "HDO.I9", "HDO.A1", "HDO.A2", "HDO.A3", 
                                                                  "HDO.A4", "HDO.A5", "HDO.A6","HDO.A7", "LDO.A1", "LDO.A2", "LDO.A3", 
                                                                  "LDO.A4", "LDO.A5", "LDO.A6", "LDO.A7", "LDO.A8", "LDO.A9", "LDO.A10", 
                                                                  "LDO.A11", "LDO.I1", "LDO.I2", "LDO.I3", "LDO.I4", "LDO.I5", "LDO.I6",
                                                                  "LDO.I7", "LDO.I8", "LDO.I9",  "LDO.I10")]
#plot
plot.parcoord.BN <- function(clust.dat){
  ggparcoord(clust.dat, columns = 3:34, showPoints = FALSE, alphaLines = 0.3) + 
    theme_bw() + theme(legend.position = "bottom", axis.title = element_text(size = 14), axis.text = element_text(size = 9), legend.text = element_text(size = 11)) + 
    xlab("Sample") + ylab("Log2 Median-Centered Gene Expression") + geom_line(color = "#2a9d8f", alpha = 0.3)
}

plot.parcoord.BA <- function(clust.dat){
  ggparcoord(clust.dat, columns = 3:39, showPoints = FALSE, alphaLines = 0.3) + 
    theme_bw() + theme(legend.position = "bottom", axis.title = element_text(size = 14), axis.text = element_text(size = 9), legend.text = element_text(size = 11)) + 
    xlab("Sample") + ylab("Log2 Median-Centered Gene Expression") + geom_line(color = "#2a9d8f", alpha = 0.3)
}

plot.parcoord.BN(subset(BN.DEG.TPM.log.center.cluster, Cluster == 1))
plot.parcoord.BN(subset(BN.DEG.TPM.log.center.cluster, Cluster == 2))
plot.parcoord.BN(subset(BN.DEG.TPM.log.center.cluster, Cluster == 3))
plot.parcoord.BN(subset(BN.DEG.TPM.log.center.cluster, Cluster == 4))
plot.parcoord.BN(subset(BN.DEG.TPM.log.center.cluster, Cluster == 5))
plot.parcoord.BN(subset(BN.DEG.TPM.log.center.cluster, Cluster == 6))
plot.parcoord.BN(subset(BN.DEG.TPM.log.center.cluster, Cluster == 7))
plot.parcoord.BN(subset(BN.DEG.TPM.log.center.cluster, Cluster == 8))
plot.parcoord.BN(subset(BN.DEG.TPM.log.center.cluster, Cluster == 9))

plot.parcoord.BA(subset(BA.DEG.TPM.log.center.cluster, Cluster == 1))
plot.parcoord.BA(subset(BA.DEG.TPM.log.center.cluster, Cluster == 2))
plot.parcoord.BA(subset(BA.DEG.TPM.log.center.cluster, Cluster == 3))
plot.parcoord.BA(subset(BA.DEG.TPM.log.center.cluster, Cluster == 4))
plot.parcoord.BA(subset(BA.DEG.TPM.log.center.cluster, Cluster == 5))
plot.parcoord.BA(subset(BA.DEG.TPM.log.center.cluster, Cluster == 6))
plot.parcoord.BA(subset(BA.DEG.TPM.log.center.cluster, Cluster == 7))
plot.parcoord.BA(subset(BA.DEG.TPM.log.center.cluster, Cluster == 8))
plot.parcoord.BA(subset(BA.DEG.TPM.log.center.cluster, Cluster == 9))

##extra code for other typse of heatmaps 

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

##old way of doing upset plots 
#get list of DEG
BN.StF.vs.StP.DEG.list <- rownames(BN.StF.vs.StP.DEG)
BN.StF.vs.SwF.DEG.list <- rownames(BN.StF.vs.SwF.DEG)
BN.StP.vs.SwP.DEG.list <- rownames(BN.StP.vs.SwP.DEG)
BN.SwF.vs.SwP.DEG.list <- rownames(BN.SwF.vs.SwP.DEG)

BN.DEG.list.upset <- data.frame("H.DO-Imm. vs H.DO-Acc." = BN.StF.vs.StP.DEG.list, 
                                "L.DO-Imm. vs L.DO-Acc." = BN.SwF.vs.SwP.DEG.list,
                                "H.DO-Imm. vs L.DO-Imm." = BN.StF.vs.SwF.DEG.list,
                                "H.DO-Acc. vs L.DO-Acc." = BN.StP.vs.SwP.DEG.list)

BA.StF.vs.StP.DEG.list <- rownames(BA.StF.vs.StP.DEG)
BA.StF.vs.SwF.DEG.list <- rownames(BA.StF.vs.SwF.DEG)
BA.StP.vs.SwP.DEG.list <- rownames(BA.StP.vs.SwP.DEG)
BA.SwF.vs.SwP.DEG.list <- rownames(BA.SwF.vs.SwP.DEG)

BA.DEG.list.upset <- list("H.DO-Imm. vs H.DO-Acc." = BA.StF.vs.StP.DEG.list, 
                          "L.DO-Imm. vs L.DO-Acc." = BA.SwF.vs.SwP.DEG.list,
                          "H.DO-Imm. vs L.DO-Imm." = BA.StF.vs.SwF.DEG.list,
                          "H.DO-Acc. vs L.DO-Acc." = BA.StP.vs.SwP.DEG.list)

BN.upset <- upset(fromList(BN.DEG.list.upset), nsets = 4, order.by = "freq", main.bar.color = "#27187E",
                  matrix.color = "#0081a7", sets.bar.color = "#2a9d8f", text.scale = c(1.7,1.5,1.7,1.5,1.5,1.2))
BN.upset

BA.upset <- upset(fromList(BA.DEG.list.upset), nsets = 4, order.by = "freq", main.bar.color = "#27187E",
                  matrix.color = "#0081a7", sets.bar.color = "#2a9d8f", text.scale = c(1.7,1.5,1.7,1.5,1.5,1.2))
BA.upset

#original way of doing gene clusters
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

#find gene clusters by cutting tree at 30% height 
find.gene.clust <- function(DEG.matrix, hclust.dat){
  gene.clust <- cutree(tree = as.dendrogram(hclust.dat), h = (0.18*max(hclust.dat$height))) 
  print(max(gene.clust))
  return(gene.clust)
}

BN.gene.clust.dat <- find.gene.clust(degBN.log.center, BN.DEG.gene.hclust)
BA.gene.clust.dat <- find.gene.clust(degBA.log.center, BA.DEG.gene.hclust)

#make cluster info into dataframe for annotations
gene.clust.BN <- data.frame(Cluster = BN.gene.clust.dat)
gene.clust.BN$Cluster <- as.factor(gene.clust.BN$Cluster)
gene.clust.BN$Cluster <- gsub("\\<1\\>", "Cluster1", gene.clust.BN$Cluster)
gene.clust.BN$Cluster <- gsub("\\<2\\>", "Cluster2", gene.clust.BN$Cluster)
gene.clust.BN$Cluster <- gsub("\\<3\\>", "Cluster3", gene.clust.BN$Cluster)
gene.clust.BN$Cluster <- gsub("\\<4\\>", "Cluster4", gene.clust.BN$Cluster)
gene.clust.BN$Cluster <- gsub("\\<5\\>", "Cluster5", gene.clust.BN$Cluster)
gene.clust.BN$Cluster <- gsub("\\<6\\>", "Cluster6", gene.clust.BN$Cluster)
gene.clust.BN$Cluster <- gsub("\\<7\\>", "Cluster7", gene.clust.BN$Cluster)
gene.clust.BN$Cluster <- gsub("\\<8\\>", "Cluster8", gene.clust.BN$Cluster)
gene.clust.BN$Cluster <- gsub("\\<9\\>", "Cluster9", gene.clust.BN$Cluster)
gene.clust.BN$Cluster <- gsub("\\<10\\>", "Cluster10", gene.clust.BN$Cluster)
gene.clust.BN$Cluster <- gsub("\\<11\\>", "Cluster11", gene.clust.BN$Cluster)
gene.clust.BN$Cluster <- gsub("\\<12\\>", "Cluster12", gene.clust.BN$Cluster)
gene.clust.BN$Cluster <- gsub("\\<13\\>", "Cluster13", gene.clust.BN$Cluster)
gene.clust.BN$Cluster <- gsub("\\<14\\>", "Cluster14", gene.clust.BN$Cluster)
gene.clust.BN$Cluster <- gsub("\\<15\\>", "Cluster15", gene.clust.BN$Cluster)

gene.clust.BA <- data.frame(Cluster = BA.gene.clust.dat)
gene.clust.BA$Cluster <- as.factor(gene.clust.BA$Cluster)
gene.clust.BA$Cluster <- gsub("\\<1\\>", "Cluster1", gene.clust.BA$Cluster)
gene.clust.BA$Cluster <- gsub("\\<2\\>", "Cluster2", gene.clust.BA$Cluster)
gene.clust.BA$Cluster <- gsub("\\<3\\>", "Cluster3", gene.clust.BA$Cluster)
gene.clust.BA$Cluster <- gsub("\\<4\\>", "Cluster4", gene.clust.BA$Cluster)
gene.clust.BA$Cluster <- gsub("\\<5\\>", "Cluster5", gene.clust.BA$Cluster)
gene.clust.BA$Cluster <- gsub("\\<6\\>", "Cluster6", gene.clust.BA$Cluster)
gene.clust.BA$Cluster <- gsub("\\<7\\>", "Cluster7", gene.clust.BA$Cluster)
gene.clust.BA$Cluster <- gsub("\\<8\\>", "Cluster8", gene.clust.BA$Cluster)
gene.clust.BA$Cluster <- gsub("\\<9\\>", "Cluster9", gene.clust.BA$Cluster)
gene.clust.BA$Cluster <- gsub("\\<10\\>", "Cluster10", gene.clust.BA$Cluster)
gene.clust.BA$Cluster <- gsub("\\<11\\>", "Cluster11", gene.clust.BA$Cluster)

#set colors for annotations
BN.clust.anno.colors <- list(
  Group = c("H-DO,Imm." = "#2a9d8f", "H-DO,Acc." = "#e9c46a",
            "L-DO,Imm." = "#f4a261", "L-DO,Acc." = "#e76f51"),
  Cluster = c(Cluster1 = "#780000", Cluster2 = "#ea5545", Cluster3 = "#f46a9b", Cluster4 = "#fb5607", Cluster5 =  "#ef9b20",
              Cluster6 = "#ffd500", Cluster7 = "#b5e48c", Cluster8 =  "#87bc45", Cluster9 = "#006400", 
              Cluster10 = "#80ced7", Cluster11 = "#27aeef", Cluster12 = "#00509d", Cluster13 = "#b33dc6", 
              Cluster14 = "#9d4edd", Cluster15 = "#5a189a"))

BA.clust.anno.colors <- list(
  Group = c("H-DO,Imm." = "#2a9d8f", "H-DO,Acc." = "#e9c46a",
            "L-DO,Imm." = "#f4a261", "L-DO,Acc." = "#e76f51"),
  Cluster = c(Cluster1 = "#780000", Cluster2 = "#ea5545", Cluster3 = "#f46a9b", Cluster4 = "#fb5607", Cluster5 =  "#ef9b20",
              Cluster6 = "#ffd500", Cluster7 = "#b5e48c", Cluster8 =  "#87bc45", Cluster9 = "#006400", 
              Cluster10 = "#80ced7", Cluster11 = "#27aeef"))

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

tiff("BN_heatmap_DEG.tiff", units="in", width = 7, height = 8, res = 600)
BN.DEG.heatmap.TMM.log.center <- make.heatmap(BN.DEG.log.center, sample.info.BN, gene.clust.BN, BN.DEG.sample.hclust, BN.DEG.gene.hclust, BN.DEG.breaks, BN.clust.anno.colors)
BN.DEG.heatmap.TMM.log.center
dev.off()

tiff("BA_heatmap_DEG.tiff", units="in", width = 7, height = 8, res = 600)
BA.DEG.heatmap.TMM.log.center <- make.heatmap(BA.DEG.log.center, sample.info.BA, gene.clust.BA, BA.DEG.sample.hclust, BA.DEG.gene.hclust, BA.DEG.breaks, BA.clust.anno.colors)
BA.DEG.heatmap.TMM.log.center
dev.off()
