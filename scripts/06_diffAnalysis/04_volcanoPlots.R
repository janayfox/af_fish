#################################
### Goal: Make volcano plots
### Author: Janay Fox
### R script
##############################

## Set up ##
library(ggplot2)

#read in data 
BN.StF.vs.StP.result <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.stream_field_vs_stream_pond.edgeR.DE_results")
BN.StF.vs.SwF.result <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.stream_field_vs_swamp_field.edgeR.DE_results")
BN.StP.vs.SwP.result <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.stream_pond_vs_swamp_pond.edgeR.DE_results")
BN.SwF.vs.SwP.result <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.swamp_field_vs_swamp_pond.edgeR.DE_results")

BA.StF.vs.StP.result <- read.table("./data/DEG/BA/salmon/BA_bf_new_sal.gene.counts.matrix.stream_field_vs_stream_pond.edgeR.DE_results")
BA.StF.vs.SwF.result <- read.table("./data/DEG/BA/salmon/BA_bf_new_sal.gene.counts.matrix.stream_field_vs_swamp_field.edgeR.DE_results")
BA.StP.vs.SwP.result <- read.table("./data/DEG/BA/salmon/BA_bf_new_sal.gene.counts.matrix.stream_pond_vs_swamp_pond.edgeR.DE_results")
BA.SwF.vs.SwP.result <- read.table("./data/DEG/BA/salmon/BA_bf_new_sal.gene.counts.matrix.swamp_field_vs_swamp_pond.edgeR.DE_results")

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

## Plot volcano plots ##
plot.volcano <- function(data){
  ggplot(data = data, aes(x = logFC, y = -log2(FDR), col = diffExpr)) + geom_point(alpha = 0.7) + 
    theme_bw() + scale_color_manual(values =  c("#e76f51", "black", "#82DCDA")) + xlab("log2(FC)") + ylab("-log2(FDR)") +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), legend.position = "none")
}

BN.StF.vs.StP.volcano <- plot.volcano(BN.StF.vs.StP.result)
BN.StF.vs.StP.volcano

BN.StF.vs.SwF.volcano <- plot.volcano(BN.StF.vs.SwF.result)
BN.StF.vs.SwF.volcano

BN.StP.vs.SwP.volcano <- plot.volcano(BN.StP.vs.SwP.result)
BN.StP.vs.SwP.volcano

BN.SwF.vs.SwP.volcano <- plot.volcano(BN.SwF.vs.SwP.result)
BN.SwF.vs.SwP.volcano

BA.StF.vs.StP.volcano <- plot.volcano(BA.StF.vs.StP.result)
BA.StF.vs.StP.volcano

BA.StF.vs.SwF.volcano <- plot.volcano(BA.StF.vs.SwF.result)
BA.StF.vs.SwF.volcano

BA.StP.vs.SwP.volcano <- plot.volcano(BA.StP.vs.SwP.result)
BA.StP.vs.SwP.volcano

BA.SwF.vs.SwP.volcano <- plot.volcano(BA.SwF.vs.SwP.result)
BA.SwF.vs.SwP.volcano
