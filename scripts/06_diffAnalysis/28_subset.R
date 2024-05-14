# # plot scatterplot 
# #adjust sign on one comparison so that it matches actual relationship 
# BN.FC.overlap$logFC_SwFvsSwP <- -1 * BN.FC.overlap$logFC_SwFvsSwP
# BA.FC.overlap$logFC_SwFvsSwP <- -1 * BA.FC.overlap$logFC_SwFvsSwP
# 
# scat.plot <- function(data){
#   ggplot(data, aes(x = logFC_SwFvsSwP, y = logFC_StFvsSwF)) + geom_point(alpha = 0.6, colour = "#27187E") + theme_bw() +
#     labs(x = paste0("Log2 FC Plastic change in L-DO Source" ), 
#          y = paste0("Log2 FC Evolutionary change", )) + 
#     geom_vline(xintercept = 0, linewidth = 0.5) + geom_hline(yintercept = 0, linewidth = 0.5) +
#     theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), legend.position = "none") 
# }

#load packages
library(edgeR)

#read in data 
data = read.table("./data/totalExpr/BN_bf_new_sal.gene.counts.matrix", header=T, row.names=1, com='')

BN.StF.vs.StP.DEG <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.stream_field_vs_stream_pond.edgeR.DE_results.P0.01_C2.DE.subset")
BN.StF.vs.SwF.DEG <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.stream_field_vs_swamp_field.edgeR.DE_results.P0.01_C2.DE.subset")
BN.StP.vs.SwP.DEG <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.stream_pond_vs_swamp_pond.edgeR.DE_results.P0.01_C2.DE.subset")
BN.SwF.vs.SwP.DEG <- read.table("./data/DEG/BN/salmon/BN_bf_new_sal.gene.counts.matrix.swamp_field_vs_swamp_pond.edgeR.DE_results.P0.01_C2.DE.subset")

#functions 
perc_overlap <- function(DEG1, DEG2) {
  # Determine the number of common rows
  common_rows <- intersect(rownames(DEG1), rownames(DEG2))
  num_common_rows <- length(common_rows)
  
  # Calculate the total number of rows in each dataframe
  total_rows_df1 <- nrow(DEG1)
  total_rows_df2 <- nrow(DEG2)
  
  # Calculate the percentage of rows that overlap
  percentage_overlap <- (num_common_rows / total_rows_df1) * 100
  
  # Print the result
  return(percentage_overlap)
}

## rerun DEG analysis by subseting 10000 times##
#stream pond and swamp pond
set.seed(113)
nsim <- 100
num.DEG.HAvsLA <- numeric(nsim)
perc.overlap.HAvsLA <- numeric(nsim)

for (i in 1:nsim) {
  #do DEG analysis
  col_ordering = c(26,12,2,1,13,24,30,22,11,20,32,31,21,8,4)
  rnaseqMatrix = data[,col_ordering]
  rnaseqMatrix = round(rnaseqMatrix)
  rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 7,]
  conditions = factor(c(rep("stream_pond", 8), rep("swamp_pond", 7)))
  
  exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
  exp_study = calcNormFactors(exp_study)
  exp_study = estimateDisp(exp_study)
  et = exactTest(exp_study, pair=c("stream_pond", "swamp_pond"))
  tTags = topTags(et,n=NULL)
  result_table = tTags$table
  result_table = data.frame(sampleA="stream_pond", sampleB="swamp_pond", result_table)
  result_table$logFC = -1 * result_table$logFC
  
  #find number of DEGs 
  DEGs_HAvsLA <- subset(result_table, abs(logFC) >= 2 & FDR <= 0.01)
  num.DEG.HAvsLA[i] <- nrow(DEGs_HAvsLA)
  
  #calculate percentage overlap 
  perc.overlap.HAvsLA[i] <- perc_overlap(BN.StP.vs.SwP.DEG, DEGs_HAvsLA)
}

#swamp field and swamp pond
num.DEG.LIvsLA <- numeric(nsim)
perc.overlap.LIvsLA <- numeric(nsim)

for (i in 1:nsim) {
  #do DEG analysis
  col_ordering = c(19,10,3,5,9,23,25,29,11,20,32,31,21,8,4)
  rnaseqMatrix = data[,col_ordering]
  rnaseqMatrix = round(rnaseqMatrix)
  rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 7,]
  conditions = factor(c(rep("swamp_field", 8), rep("swamp_pond", 7)))
  
  #randomly subset rows 
  rnaseqMatrix <- rnaseqMatrix[sample(nrow(rnaseqMatrix), 30000), ]
  
  exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
  exp_study = calcNormFactors(exp_study)
  exp_study = estimateDisp(exp_study)
  et = exactTest(exp_study, pair=c("swamp_field", "swamp_pond"))
  tTags = topTags(et,n=NULL)
  result_table = tTags$table
  result_table = data.frame(sampleA="swamp_field", sampleB="swamp_pond", result_table)
  result_table$logFC = -1 * result_table$logFC
  
  #find number of DEGs 
  DEGs_LIvsLA <- subset(result_table, abs(logFC) >= 2 & FDR <= 0.01)
  num.DEG.LIvsLA[i] <- nrow(DEGs_LIvsLA)
  
  #calculate percentage overlap 
  perc.overlap.LIvsLA[i] <- perc_overlap(BN.SwF.vs.SwP.DEG, DEGs_LIvsLA)
}

#stream field and stream pond
num.DEG.HIvsHA <- numeric(nsim)
perc.overlap.HIvsHA <- numeric(nsim)

for (i in 1:nsim) {
  #do DEG analysis
  col_ordering = c(18,14,7,6,15,17,27,16,28,26,12,2,1,13,24,30,22)
  rnaseqMatrix = data[,col_ordering]
  rnaseqMatrix = round(rnaseqMatrix)
  rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 7,]
  conditions = factor(c(rep("stream_field", 9), rep("stream_pond", 8)))
  
  #randomly subset rows 
  rnaseqMatrix <- rnaseqMatrix[sample(nrow(rnaseqMatrix), 30000), ]
  
  exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
  exp_study = calcNormFactors(exp_study)
  exp_study = estimateDisp(exp_study)
  et = exactTest(exp_study, pair=c("stream_field", "stream_pond"))
  tTags = topTags(et,n=NULL)
  result_table = tTags$table
  result_table = data.frame(sampleA="stream_field", sampleB="stream_pond", result_table)
  result_table$logFC = -1 * result_table$logFC
  
  #find number of DEGs 
  DEGs_HIvsHA <- subset(result_table, abs(logFC) >= 2 & FDR <= 0.01)
  num.DEG.HIvsHA[1] <- nrow(DEGs_HIvsHA)
  
  #calculate percentage overlap 
  perc.overlap.HIvsHA[i] <- perc_overlap(BN.StF.vs.StP.DEG, DEGs_HIvsHA)
}

missing_HAvsLA <- setdiff(rownames(DEGs_HAvsLA), rownames(BN.StP.vs.SwP.DEG))
missing_HIvsHA <- setdiff(rownames(DEGs_HIvsHA), rownames(BN.StF.vs.StP.DEG))
missing_LIvsLA <- setdiff(rownames(DEGs_LIvsLA), rownames(BN.SwF.vs.SwP.DEG))

