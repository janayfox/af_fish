######################################################################
### Goal: PCA and cluster analysis on DEG and total gene expression
### Author: Janay Fox
### R script
#####################################################################

## Set up ##
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(tidyverse)
library(rstatix)

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
BN.FC.overlap <- BN.FC[rownames(BN.FC) %in% intersect(rownames(BN.SwF.vs.SwP.DEG), rownames(BN.StP.vs.SwP.DEG)),]
BA.FC.overlap <- BA.FC[rownames(BA.FC) %in% intersect(rownames(BA.SwF.vs.SwP.DEG), rownames(BA.StP.vs.SwP.DEG)),]

#remove unneeded columns 
BN.FC.overlap <- BN.FC.overlap[,c(3,4)]
BA.FC.overlap <- BA.FC.overlap[,c(3,4)]

#remove rows with NAs
BN.FC.overlap <- na.omit(BN.FC.overlap)
BA.FC.overlap <- na.omit(BA.FC.overlap)

#add column for same vs different direction 
BN.FC.overlap$direction <- ifelse(BN.FC.overlap$logFC_StPvsSwP > 0 & BN.FC.overlap$logFC_SwFvsSwP < 0 | BN.FC.overlap$logFC_StPvsSwP < 0 & BN.FC.overlap$logFC_SwFvsSwP > 0, "Same", "Opposite")
BA.FC.overlap$direction <- ifelse(BA.FC.overlap$logFC_StPvsSwP > 0 & BA.FC.overlap$logFC_SwFvsSwP < 0 | BA.FC.overlap$logFC_StPvsSwP < 0 & BA.FC.overlap$logFC_SwFvsSwP > 0, "Same", "Opposite")

#plot
BN.dir <- count(BN.FC.overlap, direction)
BN.dir <- BN.dir %>% add_row(direction = "Same", n = 0)

BA.dir <- count(BA.FC.overlap, direction)

direction.plot <- function(data){
  ggplot(data, aes(x = direction, y = n)) + geom_bar(stat = "identity", fill = "#27187E") + theme_bw() + xlab("Direction") + ylab("Number of Genes")  +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11))
}


BN.direction.plot <- direction.plot(BN.dir)
BN.direction.plot

BA.direction.plot <- direction.plot(BA.dir)
BA.direction.plot

ggsave("BN_direction_plot_v2.tiff", plot = BN.direction.plot, device = "tiff", width = 4, height = 6, units = "in", dpi = 600)  
ggsave("BA_direction_plot_v2.tiff", plot = BA.direction.plot, device = "tiff", width = 4, height = 6, units = "in", dpi = 600)  

tiff("direction_panel.tiff", units="in", width = 9, height = 9, res = 600)
direction.panel <- ggarrange(BN.direction.plot, BA.direction.plot, 
                             labels = c("A", "B"),
                             ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")

direction.panel
dev.off()


#calculate Pearson correlation between log folds
BN.corr <- cor.test(BN.FC.overlap$logFC_StPvsSwP, BN.FC.overlap$logFC_SwFvsSwP, method = "pearson")
BA.corr <- cor.test(BA.FC.overlap$logFC_StPvsSwP, BA.FC.overlap$logFC_SwFvsSwP, method = "pearson")

#run binomial exact test 
BN.binom <- binom.test(0, 40, p = 0.5, alternative = "less")
BA.binom <- binom.test(1, 274, p = 0.5, alternative = "less")

## Run permutation test ##
#remove NAs from list of all genes 
BN.FC.noNA <- na.omit(BN.FC)
BA.FC.noNA <- na.omit(BA.FC)

#resample 10000 times and calculate Pearson correlation 
set.seed(113)
nsim <- 10000
BN.res <- numeric(nsim)
for (i in 1:nsim) {
  BN.perm <- BN.FC.noNA[sample(nrow(BN.FC.noNA), size = 40, replace = FALSE),]
  BN.res[i] <- cor(BN.perm$logFC_StPvsSwP, BN.perm$logFC_SwFvsSwP, method = "pearson")
}
BN.res <- data.frame(correlation = BN.res)

BA.res <- numeric(nsim)
for (i in 1:nsim) {
  BA.perm <- BA.FC.noNA[sample(nrow(BA.FC.noNA), size = 275, replace = FALSE),]
  BA.res[i] <- cor(BA.perm$logFC_StPvsSwP, BA.perm$logFC_SwFvsSwP, method = "pearson")
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

BN.hist.plot <- hist.plot(BN.res, -0.821)
BN.hist.plot

BA.hist.plot <- hist.plot(BA.res, -0.777)
BA.hist.plot

#save plots 
ggsave("BN_histogram_pearson_v2.tiff", plot = BN.hist.plot, device = "tiff", width = 8, height = 6, units = "in", dpi = 600)  
ggsave("BA_histogram_pearson_v2.tiff", plot = BA.hist.plot, device = "tiff", width = 8, height = 6, units = "in", dpi = 600)  

#calculate permutation test p value 
BN.perm.pval <- sum(unlist(BN.res) < -0.821) / 10000
BA.perm.pval <- sum(unlist(BA.res) < -0.777) / 10000

## Compare average magnitude FC between species ##
#isolate the DEGs we want to compare 
BN.FC.plasticDEGs <- BN.FC[rownames(BN.FC) %in% rownames(BN.SwF.vs.SwP.DEG),]
BA.FC.plasticDEGs <- BA.FC[rownames(BA.FC) %in% rownames(BA.SwF.vs.SwP.DEG),]

BN.plastic <- data.frame(FC = BN.FC.plasticDEGs$logFC_SwFvsSwP, species = "EN")
BA.plastic <- data.frame(FC = BA.FC.plasticDEGs$logFC_SwFvsSwP, species = "EA")
plastic <- rbind(BN.plastic, BA.plastic)

BN.FC.evolDEGs <- BN.FC[rownames(BN.FC) %in% rownames(BN.StP.vs.SwP.DEG),]
BA.FC.evolDEGs <- BA.FC[rownames(BA.FC) %in% rownames(BA.StP.vs.SwP.DEG),]

BN.evol <- data.frame(FC = BN.FC.evolDEGs$logFC_StPvsSwP, species = "EN")
BA.evol <- data.frame(FC = BA.FC.evolDEGs$logFC_StPvsSwP, species = "EA")
evol <- rbind(BN.evol, BA.evol)

#convert to absolute values 
plastic$FC <- abs(plastic$FC)
evol$FC <- abs(evol$FC)

#plot boxplot 
FC.boxplot <- function(data){
  ggplot(data, aes(x = species, y = FC, fill = species)) + geom_boxplot() + theme_bw() +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), legend.position = "none", plot.title = element_text(size = 10, hjust = 0.5)) +
    scale_fill_manual(values = c("EN" = "#ff9b54", "EA" = "#2a9d8f")) + labs(x = "Species", y = "Log2 Fold Change") 
}

plastic.boxplot <- FC.boxplot(plastic) + ggtitle("L-I vs L-A DEGs")
plastic.boxplot 

evol.boxplot <- FC.boxplot(evol) + labs(y = NULL) + ggtitle("H-A vs L-A DEGs")
evol.boxplot 

#save plots 
ggsave("plastic_boxplot.tiff", plot = plastic.boxplot, device = "tiff", width = 4, height = 7, units = "in", dpi = 600)  
ggsave("evol_boxplot.tiff", plot = evol.boxplot, device = "tiff", width = 4, height = 7, units = "in", dpi = 600)  

#run Mann Whitney U to test for difference 
#get summary statistics
plastic %>% group_by(species) %>% get_summary_stats(FC, type = "median_iqr")
evol %>% group_by(species) %>% get_summary_stats(FC, type = "median_iqr")

#test for normality 
plastic %>% group_by(species) %>% summarise(`W Stat` = shapiro.test(FC)$statistic,
                                            p.value = shapiro.test(FC)$p.value)

evol %>% group_by(species) %>% summarise(`W Stat` = shapiro.test(FC)$statistic,
                                         p.value = shapiro.test(FC)$p.value)

#perform Mann-Whitney U test 
plastic.manU <- wilcox.test(FC ~ species, data = plastic, exact = TRUE, conf.int = TRUE, paired = FALSE)
plastic.manU

evol.manU <- wilcox.test(FC ~ species, data = evol, exact = TRUE, conf.int = TRUE, paired = FALSE)
evol.manU

#calculate effect size 
plastic %>% wilcox_effsize(FC ~ species)
evol %>% wilcox_effsize(FC ~ species)

## Run Chi square tests on proportions of DEGs ## 
#create contigency tables 
BN.chi.data <- data.frame("DEG" = c(516, 436, 218, 63), "NOT_DEG" = c(67846, 47972, 63209, 54717))
rownames(BN.chi.data) <- c("HDOIvsHDOA", "LDOIvsLDOA", "HDOIvsLDOI", "HDOAvsLDOA")
BN.chi.data <- t(BN.chi.data)

BA.chi.data <- data.frame("DEG" = c(162, 604, 185, 320), "NOT_DEG" = c(26854, 30661, 27630, 30228))
rownames(BA.chi.data) <- c("HDOIvsHDOA", "LDOIvsLDOA", "HDOIvsLDOI", "HDOAvsLDOA")
BA.chi.data <- t(BA.chi.data)

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

BN.evol.prop.data <- data.frame("evol" = c(63,40), "not_evol" = c(54717,696))
rownames(BN.evol.prop.data) <- c("all_genes", "plastic_genes")
BN.evol.prop.data <- t(BN.evol.prop.data)

BA.evol.prop.data <- data.frame("evol" = c(320,275), "not_evol" = c(30228,324))
rownames(BA.evol.prop.data) <- c("all_genes", "plastic_genes")
BA.evol.prop.data <- t(BA.evol.prop.data)

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

BN.evol.prop.chi <- chisq.test(BN.evol.prop.data, correct = TRUE)
BN.evol.prop.chi

BA.evol.prop.chi <- chisq.test(BA.evol.prop.data, correct = TRUE)
BA.evol.prop.chi

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
combined.panel <- BN.hist.plot + BA.hist.plot + plot_annotation(tag_levels = "A")
combined.panel

ggsave("correlation_combined_panel.tiff", plot = combined.panel, device = "tiff", width = 6, height = 3, units = "in", dpi = 600)

mann.whitney.panel <- plastic.boxplot + evol.boxplot + plot_annotation(tag_levels = "A")
mann.whitney.panel

ggsave("foldChange_panel.tiff", plot = mann.whitney.panel, device = "tiff", width = 6, height = 4, units = "in", dpi = 600)
