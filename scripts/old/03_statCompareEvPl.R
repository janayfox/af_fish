######################################################################
### Goal: PCA and cluster analysis on DEG and total gene expression
### Author: Janay Fox
### R script
#####################################################################

## Set up ##
library(edgeR)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

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

## Compare direction of DEGs
#prepare data to combine 
colnames(BN.StF.vs.StP.result)[3] <- "logFC_StFvsStP"
BN.StF.vs.StP.result <- BN.StF.vs.StP.result[3]

colnames(BN.StF.vs.SwF.result)[3] <- "logFC_StFvsSwF"
BN.StF.vs.SwF.result <- BN.StF.vs.SwF.result[3]

colnames(BN.StP.vs.SwP.result)[3] <- "logFC_StPvsSwP"
BN.StP.vs.SwP.result <- BN.StP.vs.SwP.result[3]

colnames(BN.SwF.vs.SwP.result)[3] <- "logFC_SwFvsSwP"
BN.SwF.vs.SwP.result <- BN.SwF.vs.SwP.result[3]

colnames(BA.StF.vs.StP.result)[3] <- "logFC_StFvsStP"
BA.StF.vs.StP.result <- BA.StF.vs.StP.result[3]

colnames(BA.StF.vs.SwF.result)[3] <- "logFC_StFvsSwF"
BA.StF.vs.SwF.result <- BA.StF.vs.SwF.result[3]

colnames(BA.StP.vs.SwP.result)[3] <- "logFC_StPvsSwP"
BA.StP.vs.SwP.result <- BA.StP.vs.SwP.result[3]

colnames(BA.SwF.vs.SwP.result)[3] <- "logFC_SwFvsSwP"
BA.SwF.vs.SwP.result <- BA.SwF.vs.SwP.result[3]

#merge dataframes 
BN.FC <- transform(merge(BN.StF.vs.StP.result, BN.StF.vs.SwF.result, by ="row.names", all.x = TRUE, all.y = TRUE), row.names = Row.names, Row.names = NULL)
BN.FC <- transform(merge(BN.FC, BN.StP.vs.SwP.result, by ="row.names", all.x = TRUE, all.y = TRUE), row.names = Row.names, Row.names = NULL)
BN.FC <- transform(merge(BN.FC, BN.SwF.vs.SwP.result, by ="row.names", all.x = TRUE, all.y = TRUE), row.names = Row.names, Row.names = NULL)

BA.FC <- transform(merge(BA.StF.vs.StP.result, BA.StF.vs.SwF.result, by ="row.names", all.x = TRUE, all.y = TRUE), row.names = Row.names, Row.names = NULL)
BA.FC <- transform(merge(BA.FC, BA.StP.vs.SwP.result, by ="row.names", all.x = TRUE, all.y = TRUE), row.names = Row.names, Row.names = NULL)
BA.FC <- transform(merge(BA.FC, BA.SwF.vs.SwP.result, by ="row.names", all.x = TRUE, all.y = TRUE), row.names = Row.names, Row.names = NULL)

#divide into lists for LDO-I vs LDO-A DEGs
BN.FC.SwFvsSwP <- BN.FC[rownames(BN.FC) %in% rownames(BN.SwF.vs.SwP.DEG),]
BA.FC.SwFvsSwP <- BA.FC[rownames(BA.FC) %in% rownames(BA.SwF.vs.SwP.DEG),]

#remove unneeded columns 
BN.FC.SwFvsSwP <- BN.FC.SwFvsSwP[,c(2,4)]
BA.FC.SwFvsSwP <- BA.FC.SwFvsSwP[,c(2,4)]

#remove rows with NAs
BN.FC.SwFvsSwP <- na.omit(BN.FC.SwFvsSwP)
BA.FC.SwFvsSwP <- na.omit(BA.FC.SwFvsSwP)

#add column for same vs different direction 
BN.FC.SwFvsSwP$direction <- ifelse(BN.FC.SwFvsSwP$logFC_StFvsSwF > 0 & BN.FC.SwFvsSwP$logFC_SwFvsSwP < 0 | BN.FC.SwFvsSwP$logFC_StFvsSwF < 0 & BN.FC.SwFvsSwP$logFC_SwFvsSwP > 0, "Same", "Opposite")

BA.FC.SwFvsSwP$direction <- ifelse(BA.FC.SwFvsSwP$logFC_StFvsSwF > 0 & BA.FC.SwFvsSwP$logFC_SwFvsSwP < 0 | BA.FC.SwFvsSwP$logFC_StFvsSwF < 0 & BA.FC.SwFvsSwP$logFC_SwFvsSwP > 0, "Same", "Opposite")

#plot
direction.plot <- function(data){
  ggplot(data, aes(direction)) + geom_bar(stat = "count", fill = "#27187E") + theme_bw() + xlab("Direction") + ylab("Number of Genes")  +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11))
}

BN.direction.plot <- direction.plot(BN.FC.SwFvsSwP)
BN.direction.plot

BA.direction.plot <- direction.plot(BA.FC.SwFvsSwP)
BA.direction.plot

ggsave("BN_direction_plot.tiff", plot = BN.direction.plot, device = "tiff", width = 4, height = 6, units = "in", dpi = 600)  
ggsave("BA_direction_plot.tiff", plot = BA.direction.plot, device = "tiff", width = 4, height = 6, units = "in", dpi = 600)  

tiff("direction_panel.tiff", units="in", width = 9, height = 9, res = 600)
direction.panel <- ggarrange(BN.direction.plot, BA.direction.plot, 
                       labels = c("A", "B"),
                       ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")

direction.panel
dev.off()

#make list of genes in same direction 
BN.sameDir.genes <- rownames(BN.FC.SwFvsSwP[BN.FC.SwFvsSwP$direction == "Same",])
BN.sameDir.genes.up <- rownames(BN.FC.SwFvsSwP[BN.FC.SwFvsSwP$direction == "Same" & BN.FC.SwFvsSwP$logFC_StFvsSwF > 0,])
BN.sameDir.genes.down <- rownames(BN.FC.SwFvsSwP[BN.FC.SwFvsSwP$direction == "Same" & BN.FC.SwFvsSwP$logFC_StFvsSwF < 0,])

BA.sameDir.genes <- rownames(BA.FC.SwFvsSwP[BA.FC.SwFvsSwP$direction == "Same",])
BA.sameDir.genes.up <- rownames(BA.FC.SwFvsSwP[BA.FC.SwFvsSwP$direction == "Same" & BA.FC.SwFvsSwP$logFC_StFvsSwF > 0,])
BA.sameDir.genes.down <- rownames(BA.FC.SwFvsSwP[BA.FC.SwFvsSwP$direction == "Same" & BA.FC.SwFvsSwP$logFC_StFvsSwF < 0,])

BN.same.genes <- data.frame(factor = "BN_sameDir", gene_id = BN.sameDir.genes)
BN.same.genes.up <- data.frame(factor = "BN_sameDir_up", gene_id = BN.sameDir.genes.up)
BN.same.genes.down <- data.frame(factor = "BN_sameDir_down", gene_id = BN.sameDir.genes.down)

BA.same.genes <- data.frame(factor = "BA_sameDir", gene_id = BA.sameDir.genes)
BA.same.genes.up <- data.frame(factor = "BA_sameDir_up", gene_id = BA.sameDir.genes.up)
BA.same.genes.down <- data.frame(factor = "BA_sameDir_down", gene_id = BA.sameDir.genes.down)

write.table(BN.same.genes, "./data/goSeqData/BN_same_direction_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BN.same.genes.up, "./data/goSeqData/BN_same_direction_up_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BN.same.genes.down, "./data/goSeqData/BN_same_direction_down_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(BA.same.genes, "./data/goSeqData/BA_same_direction_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BA.same.genes.up, "./data/goSeqData/BA_same_direction_up_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(BA.same.genes.down, "./data/goSeqData/BA_same_direction_down_factor.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

#calculate Pearson correlation between log folds
BN.corr <- cor.test(BN.FC.SwFvsSwP$logFC_StFvsSwF, BN.FC.SwFvsSwP$logFC_SwFvsSwP, method = "pearson")
BA.corr <- cor.test(BA.FC.SwFvsSwP$logFC_StFvsSwF, BA.FC.SwFvsSwP$logFC_SwFvsSwP, method = "pearson")

#run binomial exact test 
BN.binom <- binom.test(252, 334, p = 0.5, alternative = "greater")
BA.binom <- binom.test(209, 249, p = 0.5, alternative = "greater")

## Run permutation test ##
#remove NAs from list of all genes 
BN.FC <- na.omit(BN.FC)
BA.FC <- na.omit(BA.FC)

#resample 10000 times and calculate Pearson correlation 
set.seed(113)
nsim <- 10000
BN.res <- numeric(nsim)
for (i in 1:nsim) {
  BN.perm <- BN.FC[sample(nrow(BN.FC), size = 334, replace = FALSE),]
  BN.res[i] <- cor(BN.perm$logFC_StFvsSwF, BN.perm$logFC_SwFvsSwP, method = "pearson")
}
BN.res <- data.frame(correlation = BN.res)

BA.res <- numeric(nsim)
for (i in 1:nsim) {
  BA.perm <- BA.FC[sample(nrow(BA.FC), size = 249, replace = FALSE),]
  BA.res[i] <- cor(BA.perm$logFC_StFvsSwF, BA.perm$logFC_SwFvsSwP, method = "pearson")
}
BA.res <- data.frame(correlation = BA.res)

#reverse sign to match actual relationship (because of the order of sample types in comparison)
BN.res$correlation <- -1 * BN.res$correlation
BA.res$correlation <- -1 * BA.res$correlation

#plot permutation results 
hist.plot <- function(data, line.dat){
  ggplot(data, aes(x=correlation)) + geom_histogram(bins = 40, fill = "#2a9d8f", color = "white") + theme_bw() + 
    xlab("Pearson Correlation") + ylab("Count") +
    geom_vline(xintercept = line.dat, linetype = "dashed", color = "#27187E", linewidth = 1.5) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11))
}

BN.hist.plot <- hist.plot(BN.res, 0.7482083)
BN.hist.plot

BA.hist.plot <- hist.plot(BA.res, 0.6834277)
BA.hist.plot

#save plots 
ggsave("BN_histogram_pearson.tiff", plot = BN.hist.plot, device = "tiff", width = 8, height = 6, units = "in", dpi = 600)  
ggsave("BA_histogram_pearson.tiff", plot = BA.hist.plot, device = "tiff", width = 8, height = 6, units = "in", dpi = 600)  

#calculate permutation test p value 
BN.perm.pval <- sum(unlist(BN.res) > 0.7482083) / 10000
BA.perm.pval <- sum(unlist(BA.res) > 0.6834277) / 10000

# plot scatterplot 
#adjust sign on one comparison so that it matches actual relationship 
BN.FC.SwFvsSwP$logFC_SwFvsSwP <- -1 * BN.FC.SwFvsSwP$logFC_SwFvsSwP
BA.FC.SwFvsSwP$logFC_SwFvsSwP <- -1 * BA.FC.SwFvsSwP$logFC_SwFvsSwP

scat.plot <- function(data){
  ggplot(data, aes(x = logFC_SwFvsSwP, y = logFC_StFvsSwF, color = direction)) + geom_point(alpha = 0.6) + theme_bw() +
    labs(x = paste0("Plastic change in L-DO Source", "\n", "Log2 FC (L-DO, Imm. vs L-DO, Acc.)"), 
         y = paste0("Evolutionary divergence", "\n", "Log2 FC (H-DO, Imm. vs L-DO, Imm.)")) + 
    geom_vline(xintercept = 0, linewidth = 0.5) + geom_hline(yintercept = 0, linewidth = 0.5) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), legend.position = "none") + 
    scale_color_manual(values = c("Same" = "#27187E", "Opposite" = "#ad2e24"))
}

BN.scat.plot <- scat.plot(BN.FC.SwFvsSwP)
BN.scat.plot

BA.scat.plot <- scat.plot(BA.FC.SwFvsSwP)
BA.scat.plot

#save plots
ggsave("BN_scatterplot.tiff", plot = BN.scat.plot, device = "tiff", width = 6, height = 6, units = "in", dpi = 600)  
ggsave("BA_scatterplot.tiff", plot = BA.scat.plot, device = "tiff", width = 6, height = 6, units = "in", dpi = 600)  

## Compare average magnitude FC between species ##
#isolate the DEGs we want to compare 
BN.sameDir <- BN.FC.SwFvsSwP[BN.FC.SwFvsSwP$direction == "Same",]
BA.sameDir <- BA.FC.SwFvsSwP[BA.FC.SwFvsSwP$direction == "Same",]
BN.sameDir.plastic <- data.frame(FC = BN.sameDir$logFC_SwFvsSwP, species = "EN")
BA.sameDir.plastic <- data.frame(FC = BA.sameDir$logFC_SwFvsSwP, species = "EA")
sameDir.plastic <- rbind(BN.sameDir.plastic, BA.sameDir.plastic)

BN.plastic <- data.frame(FC = BN.FC.SwFvsSwP$logFC_SwFvsSwP, species = "EN")
BA.plastic <- data.frame(FC = BA.FC.SwFvsSwP$logFC_SwFvsSwP, species = "EA")
plastic <- rbind(BN.plastic, BA.plastic)

BN.evol <- data.frame(FC = BN.StF.vs.SwF.DEG$logFC, species = "EN")
BA.evol <- data.frame(FC = BA.StF.vs.SwF.DEG$logFC, species = "EA")
evol <- rbind(BN.evol, BA.evol)
  
#convert to absolute values 
sameDir.plastic$FC <- abs(sameDir.plastic$FC)
plastic$FC <- abs(plastic$FC)
evol$FC <- abs(evol$FC)

#plot boxplot 
FC.boxplot <- function(data){
  ggplot(data, aes(x = species, y = FC, fill = species)) + geom_boxplot() + theme_bw() +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), legend.position = "none", plot.title = element_text(size = 10, hjust = 0.5)) +
    scale_fill_manual(values = c("EN" = "#ff9b54", "EA" = "#2a9d8f")) + labs(x = "Species", y = "Log2 Fold Change") 
}

sameDir.plastic.boxplot <- FC.boxplot(sameDir.plastic) + ggtitle("Same Direction DEGs") + labs(x = NULL)
sameDir.plastic.boxplot

plastic.boxplot <- FC.boxplot(plastic) + labs(y = NULL) + ggtitle("L-DO, Imm. vs L-DO. Acc. DEGs")
plastic.boxplot 

evol.boxplot <- FC.boxplot(evol) + labs(x = NULL, y = NULL) + ggtitle("H-DO, Imm. vs L-DO, Imm. DEGs")
evol.boxplot 

#save plots 
ggsave("sameDir_plastic_boxplot.tiff", plot = sameDir.plastic.boxplot, device = "tiff", width = 4, height = 7, units = "in", dpi = 600)  
ggsave("plastic_boxplot.tiff", plot = plastic.boxplot, device = "tiff", width = 4, height = 7, units = "in", dpi = 600)  
ggsave("evol_boxplot.tiff", plot = evol.boxplot, device = "tiff", width = 4, height = 7, units = "in", dpi = 600)  

#run Mann Whitney U to test for difference 
#test for normality 
sameDir.plastic %>% group_by(species) %>% summarise(`W Stat` = shapiro.test(FC)$statistic,
                                                    p.value = shapiro.test(FC)$p.value)

plastic %>% group_by(species) %>% summarise(`W Stat` = shapiro.test(FC)$statistic,
                                                    p.value = shapiro.test(FC)$p.value)

evol %>% group_by(species) %>% summarise(`W Stat` = shapiro.test(FC)$statistic,
                                                    p.value = shapiro.test(FC)$p.value)

#perform Mann-Whitney U test 
sameDir.manU <- wilcox.test(FC ~ species, data = sameDir.plastic, exact = TRUE, conf.int = TRUE, paired = FALSE)
sameDir.manU

plastic.manU <- wilcox.test(FC ~ species, data = plastic, exact = TRUE, conf.int = TRUE, paired = FALSE)
plastic.manU

evol.manU <- wilcox.test(FC ~ species, data = evol, exact = TRUE, conf.int = TRUE, paired = FALSE)
evol.manU

## Run Chi square tests on proportions of DEGs ## 
#create contigency tables 
LDOIvsLDOA.data <- data.frame("DEG" = c(436, 604), "NOT_DEG" = c(47972, 30661))
rownames(LDOIvsLDOA.data) <- c("EN", "EA")
LDOIvsLDOA.data <- t(LDOIvsLDOA.data)

HDOIvsLDOI.data <- data.frame("DEG" = c(218, 185), "NOT_DEG" = c(63209, 27630))
rownames(HDOIvsLDOI.data) <- c("EN", "EA")
HDOIvsLDOI.data <- t(HDOIvsLDOI.data)

HDOIvsHDOA.data <- data.frame("DEG" = c(516, 162), "NOT_DEG" = c(67846, 26854))
rownames(HDOIvsHDOA.data) <- c("EN", "EA")
HDOIvsHDOA.data <- t(HDOIvsHDOA.data)

HDOAvsLDOA.data <- data.frame("DEG" = c(63, 320), "NOT_DEG" = c(54717, 30228))
rownames(HDOAvsLDOA.data) <- c("EN", "EA")
HDOAvsLDOA.data <- t(HDOAvsLDOA.data)

BN.chi.data <- data.frame("DEG" = c(516, 436, 218, 63), "NOT_DEG" = c(67846, 47972, 63209, 54717))
rownames(BN.chi.data) <- c("HDOIvsHDOA", "LDOIvsLDOA", "HDOIvsLDOI", "HDOAvsLDOA")
BN.chi.data <- t(BN.chi.data)

BA.chi.data <- data.frame("DEG" = c(162, 604, 185, 320), "NOT_DEG" = c(26854, 30661, 27630, 30228))
rownames(BA.chi.data) <- c("HDOIvsHDOA", "LDOIvsLDOA", "HDOIvsLDOI", "HDOAvsLDOA")
BA.chi.data <- t(BA.chi.data)

#run chi test 
LDOIvsLDOA.chi <- chisq.test(LDOIvsLDOA.data, correct = TRUE)
LDOIvsLDOA.chi

HDOIvsLDOI.chi <- chisq.test(HDOIvsLDOI.data, correct = TRUE)
HDOIvsLDOI.chi

HDOIvsHDOA.chi <- chisq.test(HDOIvsHDOA.data, correct = TRUE)
HDOIvsHDOA.chi

HDOAvsLDOA.chi <- chisq.test(HDOAvsLDOA.data, correct = TRUE)
HDOAvsLDOA.chi

BN.chi <- chisq.test(BN.chi.data, correct = TRUE)
BN.chi

BA.chi <- chisq.test(BA.chi.data, correct = TRUE)
BA.chi

#run post hoc tests
BN.chi <- chisq.test(BN.chi.data[,c(1,2)], correct = TRUE)
BN.chi

BN.chi <- chisq.test(BN.chi.data[,c(1,3)], correct = TRUE)
BN.chi

BN.chi <- chisq.test(BN.chi.data[,c(1,4)], correct = TRUE)
BN.chi

BN.chi <- chisq.test(BN.chi.data[,c(2,3)], correct = TRUE)
BN.chi

BN.chi <- chisq.test(BN.chi.data[,c(2,4)], correct = TRUE)
BN.chi

BN.chi <- chisq.test(BN.chi.data[,c(3,4)], correct = TRUE)
BN.chi


BA.chi <- chisq.test(BA.chi.data[,c(1,2)], correct = TRUE)
BA.chi

BA.chi <- chisq.test(BA.chi.data[,c(1,3)], correct = TRUE)
BA.chi

BA.chi <- chisq.test(BA.chi.data[,c(1,4)], correct = TRUE)
BA.chi

BA.chi <- chisq.test(BA.chi.data[,c(2,3)], correct = TRUE)
BA.chi

BA.chi <- chisq.test(BA.chi.data[,c(2,4)], correct = TRUE)
BA.chi

BA.chi <- chisq.test(BA.chi.data[,c(3,4)], correct = TRUE)
BA.chi

## Create plot panels ## 
BN.patch <- BN.hist.plot / BN.direction.plot 
BN.panel <- BN.scat.plot + BN.patch
BN.panel <- BN.panel + plot_annotation(tag_levels = "A")

BA.patch <- BA.hist.plot / BA.direction.plot 
BA.panel <- BA.scat.plot + BA.patch
BA.panel <- BA.panel + plot_annotation(tag_levels = "A")

ggsave("BN_sameDir_panel.tiff", plot = BN.panel, device = "tiff", width = 9, height = 7, units = "in", dpi = 600)  
ggsave("BA_sameDir_panel.tiff", plot = BA.panel, device = "tiff", width = 9, height = 7, units = "in", dpi = 600)  

mann.whitney.panel <- sameDir.plastic.boxplot + plastic.boxplot + evol.boxplot + plot_annotation(tag_levels = "A")
mann.whitney.panel

ggsave("foldChange_panel.tiff", plot = mann.whitney.panel, device = "tiff", width = 8, height = 4, units = "in", dpi = 600)


