######################################################################
### Goal: Compare evolved and plastic genes and their direction V2
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
library(VennDiagram)

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

#try filtered results
degBN <- read.table("./data/DEG/filtered/BN/diffExpr.P0.01_C2.matrix")
degBA <- read.table("./data/DEG/filtered/BA/diffExpr.P0.01_C2.matrix")

allExprBN <- read.table("./data/DEG/filtered/BN/BN_bf_filtered_sal.gene.TMM.EXPR.matrix")
allExprBA <- read.table("./data/DEG/filtered/BA/BA_bf_filtered_sal.gene.TMM.EXPR.matrix")

BN.StF.vs.StP.DEG <- read.table("./data/DEG/filtered/BN/BN_bf_filtered_sal.gene.counts.matrix.stream_field_vs_stream_pond.edgeR.DE_results.P0.01_C2.DE.subset")
BN.StF.vs.SwF.DEG <- read.table("./data/DEG/filtered/BN/BN_bf_filtered_sal.gene.counts.matrix.stream_field_vs_swamp_field.edgeR.DE_results.P0.01_C2.DE.subset")
BN.StP.vs.SwP.DEG <- read.table("./data/DEG/filtered/BN/BN_bf_filtered_sal.gene.counts.matrix.stream_pond_vs_swamp_pond.edgeR.DE_results.P0.01_C2.DE.subset")
BN.SwF.vs.SwP.DEG <- read.table("./data/DEG/filtered/BN/BN_bf_filtered_sal.gene.counts.matrix.swamp_field_vs_swamp_pond.edgeR.DE_results.P0.01_C2.DE.subset")

BA.StF.vs.StP.DEG <- read.table("./data/DEG/filtered/BA/BA_bf_filtered_sal.gene.counts.matrix.stream_field_vs_stream_pond.edgeR.DE_results.P0.01_C2.DE.subset")
BA.StF.vs.SwF.DEG <- read.table("./data/DEG/filtered/BA/BA_bf_filtered_sal.gene.counts.matrix.stream_field_vs_swamp_field.edgeR.DE_results.P0.01_C2.DE.subset")
BA.StP.vs.SwP.DEG <- read.table("./data/DEG/filtered/BA/BA_bf_filtered_sal.gene.counts.matrix.stream_pond_vs_swamp_pond.edgeR.DE_results.P0.01_C2.DE.subset")
BA.SwF.vs.SwP.DEG <- read.table("./data/DEG/filtered/BA/BA_bf_filtered_sal.gene.counts.matrix.swamp_field_vs_swamp_pond.edgeR.DE_results.P0.01_C2.DE.subset")

BN.StF.vs.StP.result <- read.table("./data/DEG/filtered/BN/BN_bf_filtered_sal.gene.counts.matrix.stream_field_vs_stream_pond.edgeR.DE_results")
BN.StF.vs.SwF.result <- read.table("./data/DEG/filtered/BN/BN_bf_filtered_sal.gene.counts.matrix.stream_field_vs_swamp_field.edgeR.DE_results")
BN.StP.vs.SwP.result <- read.table("./data/DEG/filtered/BN/BN_bf_filtered_sal.gene.counts.matrix.stream_pond_vs_swamp_pond.edgeR.DE_results")
BN.SwF.vs.SwP.result <- read.table("./data/DEG/filtered/BN/BN_bf_filtered_sal.gene.counts.matrix.swamp_field_vs_swamp_pond.edgeR.DE_results")

BA.StF.vs.StP.result <- read.table("./data/DEG/filtered/BA/BA_bf_filtered_sal.gene.counts.matrix.stream_field_vs_stream_pond.edgeR.DE_results")
BA.StF.vs.SwF.result <- read.table("./data/DEG/filtered/BA/BA_bf_filtered_sal.gene.counts.matrix.stream_field_vs_swamp_field.edgeR.DE_results")
BA.StP.vs.SwP.result <- read.table("./data/DEG/filtered/BA/BA_bf_filtered_sal.gene.counts.matrix.stream_pond_vs_swamp_pond.edgeR.DE_results")
BA.SwF.vs.SwP.result <- read.table("./data/DEG/filtered/BA/BA_bf_filtered_sal.gene.counts.matrix.swamp_field_vs_swamp_pond.edgeR.DE_results")

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

#find plastic genes 
BN.FC.plastic <- BN.FC[rownames(BN.FC) %in% setdiff(rownames(BN.SwF.vs.SwP.DEG), rownames(BN.StF.vs.StP.DEG)),]
BA.FC.plastic <- BA.FC[rownames(BA.FC) %in% setdiff(rownames(BA.SwF.vs.SwP.DEG), rownames(BA.StF.vs.StP.DEG)),]

#find evovled genes 
BN.FC.evol <- BN.FC[rownames(BN.FC) %in% rownames(BN.StP.vs.SwP.DEG),]
BA.FC.evol <- BA.FC[rownames(BA.FC) %in% rownames(BA.StP.vs.SwP.DEG),]

#find overlap genes
BN.FC.overlap <- BN.FC[rownames(BN.FC) %in% intersect(rownames(BN.FC.plastic), rownames(BN.FC.evol)),]
BA.FC.overlap <- BA.FC[rownames(BA.FC) %in% intersect(rownames(BA.FC.plastic), rownames(BA.FC.evol)),]

#remove unneeded columns 
BN.FC.overlap <- BN.FC.overlap[,c(3,4)]
BA.FC.overlap <- BA.FC.overlap[,c(3,4)]

#add column for same vs different direction 
BN.FC.overlap$direction <- ifelse(BN.FC.overlap$logFC_StPvsSwP > 0 & BN.FC.overlap$logFC_SwFvsSwP < 0 | BN.FC.overlap$logFC_StPvsSwP < 0 & BN.FC.overlap$logFC_SwFvsSwP > 0, "Same", "Opposite")
BA.FC.overlap$direction <- ifelse(BA.FC.overlap$logFC_StPvsSwP > 0 & BA.FC.overlap$logFC_SwFvsSwP < 0 | BA.FC.overlap$logFC_StPvsSwP < 0 & BA.FC.overlap$logFC_SwFvsSwP > 0, "Same", "Opposite")

#calculate Pearson correlation between log folds
BN.corr <- cor.test(BN.FC.overlap$logFC_StPvsSwP, BN.FC.overlap$logFC_SwFvsSwP, method = "pearson")
BA.corr <- cor.test(BA.FC.overlap$logFC_StPvsSwP, BA.FC.overlap$logFC_SwFvsSwP, method = "pearson")

#run binomial exact test 
BN.binom <- binom.test(0, 31, p = 0.5, alternative = "less") #need to change these numbers if i update to filtered method?? 
BA.binom <- binom.test(1, 269, p = 0.5, alternative = "less")

#save lists of evolved and plastic genes 
write.table(rownames(BN.FC.plastic), file = "./data/DEG/BN_plastic_list.tsv", row.names = FALSE, sep = "/t")
write.table(rownames(BA.FC.plastic), file = "./data/DEG/BA_plastic_list.tsv", row.names = FALSE, sep = "/t")
write.table(rownames(BN.FC.evol), file = "./data/DEG/BN_evolved_list.tsv", row.names = FALSE, sep = "/t")
write.table(rownames(BA.FC.evol), file = "./data/DEG/BA_evolved_list.tsv", row.names = FALSE, sep = "/t")

## Run permutation test ##
#remove NAs from list of all genes 
BN.FC.noNA <- na.omit(BN.FC)
BA.FC.noNA <- na.omit(BA.FC)

#resample 10000 times and calculate Pearson correlation 
set.seed(113)
nsim <- 10000
BN.res <- numeric(nsim)
for (i in 1:nsim) {
  BN.perm <- BN.FC.noNA[sample(nrow(BN.FC.noNA), size = 31, replace = FALSE),]
  BN.res[i] <- cor(BN.perm$logFC_StPvsSwP, BN.perm$logFC_SwFvsSwP, method = "pearson")
}
BN.res <- data.frame(correlation = BN.res)

BA.res <- numeric(nsim)
for (i in 1:nsim) {
  BA.perm <- BA.FC.noNA[sample(nrow(BA.FC.noNA), size = 269, replace = FALSE),]
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

BN.hist.plot <- hist.plot(BN.res, -0.838)
BN.hist.plot

BA.hist.plot <- hist.plot(BA.res, -0.792)
BA.hist.plot

#save plots 
ggsave("BN_histogram_pearson_updated.tiff", plot = BN.hist.plot, device = "tiff", width = 8, height = 6, units = "in", dpi = 600)  
ggsave("BA_histogram_pearson_updated.tiff", plot = BA.hist.plot, device = "tiff", width = 8, height = 6, units = "in", dpi = 600)  

#calculate permutation test p value 
BN.perm.pval <- sum(unlist(BN.res) < -0.838) / 10000
BA.perm.pval <- sum(unlist(BA.res) < -0.792) / 10000

# plot scatterplot 
#adjust sign on one comparison so that it matches actual relationship 
BN.FC.overlap$logFC_SwFvsSwP <- -1 * BN.FC.overlap$logFC_SwFvsSwP
BA.FC.overlap$logFC_SwFvsSwP <- -1 * BA.FC.overlap$logFC_SwFvsSwP
    
scat.plot <- function(data){
  ggplot(data, aes(x = logFC_SwFvsSwP, y = logFC_StPvsSwP)) + geom_point(alpha = 0.6, colour = "#27187E") + theme_bw() +
    labs(x = paste0("Log2 FC Plastic change in L-DO source" ), 
         y = paste0("Log2 FC Evolutionary change")) + 
    geom_vline(xintercept = 0, linewidth = 0.5) + geom_hline(yintercept = 0, linewidth = 0.5) +
    theme(axis.title = element_text(size = 9), axis.text = element_text(size = 8), legend.position = "none") 
}

BN.scat.plot <- scat.plot(BN.FC.overlap)
BN.scat.plot

BA.scat.plot <- scat.plot(BA.FC.overlap)
BA.scat.plot

#save plots 
ggsave("BN_logFC_scatterplot_updated.tiff", plot = BN.scat.plot, device = "tiff", width = 5, height = 5, units = "in", dpi = 300)  
ggsave("BA_logFC_scatterplot_updated.png", plot = BA.scat.plot, device = "png", width = 5, height = 5, units = "in", dpi = 300)  

## Compare average magnitude FC between species ##
#isolate the DEGs we want to compare 
BN.plastic <- data.frame(FC = BN.FC.plastic$logFC_SwFvsSwP, species = "EN")
BA.plastic <- data.frame(FC = BA.FC.plastic$logFC_SwFvsSwP, species = "EA")
plastic <- rbind(BN.plastic, BA.plastic)

BN.evol <- data.frame(FC = BN.FC.evol$logFC_StPvsSwP, species = "EN")
BA.evol <- data.frame(FC = BA.FC.evol$logFC_StPvsSwP, species = "EA")
evol <- rbind(BN.evol, BA.evol)

#convert to absolute values 
plastic$FC <- abs(plastic$FC)
evol$FC <- abs(evol$FC)

#plot boxplot 
FC.boxplot <- function(data){
  ggplot(data, aes(x = species, y = FC, fill = species)) + geom_boxplot() + theme_bw() +
    theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), legend.position = "none", plot.title = element_text(size = 16, hjust = 0.5)) +
    scale_fill_manual(values = c("EN" = "#ff9b54", "EA" = "#2a9d8f")) + labs(x = "Species", y = "Log2 Fold Change") 
}

plastic.boxplot <- FC.boxplot(plastic) + ggtitle("Plastic DEGs")
plastic.boxplot 

evol.boxplot <- FC.boxplot(evol) + labs(y = NULL) + ggtitle("Evolved DEGs")
evol.boxplot 

#save plots 
ggsave("plastic_boxplot_updated.tiff", plot = plastic.boxplot, device = "tiff", width = 4, height = 7, units = "in", dpi = 600)  
ggsave("evol_boxplot_updated.tiff", plot = evol.boxplot, device = "tiff", width = 4, height = 7, units = "in", dpi = 600)  

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
BN.evolvsplastic.data <- data.frame("sig" = c(63,344), "not_sig" = c(54717,48064))
rownames(BN.evolvsplastic.data) <- c("evol_genes", "plastic_genes")
BN.evolvsplastic.data <- t(BN.evolvsplastic.data)

BA.evolvsplastic.data <- data.frame("sig" = c(320,556), "not_sig" = c(30228,30709))
rownames(BA.evolvsplastic.data) <- c("evol_genes", "plastic_genes")
BA.evolvsplastic.data <- t(BA.evolvsplastic.data)

BN.evol.prop.data <- data.frame("evol" = c(63,31), "not_evol" = c(54717,313))
rownames(BN.evol.prop.data) <- c("all_genes", "plastic_genes")
BN.evol.prop.data <- t(BN.evol.prop.data)

BA.evol.prop.data <- data.frame("evol" = c(320,269), "not_evol" = c(30228,287))
rownames(BA.evol.prop.data) <- c("all_genes", "plastic_genes")
BA.evol.prop.data <- t(BA.evol.prop.data)

#run chi test 
BN.evolvsplas.chi <- chisq.test(BN.evolvsplastic.data, correct = TRUE)
BN.evolvsplas.chi

BA.evolvsplas.chi <- chisq.test(BA.evolvsplastic.data, correct = TRUE)
BA.evolvsplas.chi

BN.evol.prop.chi <- chisq.test(BN.evol.prop.data, correct = TRUE)
BN.evol.prop.chi

BA.evol.prop.chi <- chisq.test(BA.evol.prop.data, correct = TRUE)
BA.evol.prop.chi

## Create plot panels ## 
combined.panel <- BN.hist.plot + BA.hist.plot + plot_annotation(tag_levels = "A")
combined.panel

ggsave("correlation_combined_panel_updated.tiff", plot = combined.panel, device = "tiff", width = 6, height = 3, units = "in", dpi = 600)

mann.whitney.panel <- plastic.boxplot + evol.boxplot + plot_annotation(tag_levels = "A")
mann.whitney.panel

ggsave("foldChange_panel_updated.tiff", plot = mann.whitney.panel, device = "tiff", width = 6, height = 4, units = "in", dpi = 600)

scatter.panel <-  BN.scat.plot + BA.scat.plot + plot_annotation(tag_levels = "A")
scatter.panel

ggsave("scatter_panel_updated.tiff", plot = scatter.panel, device = "tiff", width = 6, height = 4, units = "in", dpi = 600)

##Plot Venn Diagram
BN.plast.list <- rownames(BN.FC.plastic)
BN.evol.list <- rownames(BN.FC.evol)

BA.plast.list <- rownames(BA.FC.plastic)
BA.evol.list <- rownames(BA.FC.evol)

venn.diagram(x = list(BN.plast.list, BN.evol.list),
             category.names = c("Plastic DEGs", "Evolved DEGs"),
             filename = "BN_venn_updated.png",
             output = TRUE,
             imagetype = "png",
             height = 800,
             width = 800,
             resolution = 600,
             lwd = 2,
             lty = "blank",
             fill = c("#ff9b54", "#2a9d8f"),
             cex = 0.6,
             fontface = "bold",
             fontfamily = "sans",
             cat.cex = 0.1,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.dist = c(0.055, 0.055),
             cat.fontfamily = "sans"
             )

venn.diagram(x = list(BA.plast.list, BA.evol.list),
             category.names = c("Plastic DEGs", "Evolved DEGs"),
             filename = "BA_venn_updated.png",
             output = TRUE,
             imagetype = "png",
             height = 800,
             width = 800,
             resolution = 600,
             lwd = 2,
             lty = "blank",
             fill = c("#ff9b54", "#2a9d8f"),
             cex = 0.6,
             fontface = "bold",
             fontfamily = "sans",
             cat.cex = 0.1,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.dist = c(0.055, 0.055),
             cat.fontfamily = "sans"
)
