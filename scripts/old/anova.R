library(lme4)
BA.matrix <- readRDS("./melt_BA.RDS")
BN.matrix <- readRDS("./melt_BN.RDS")

#take random sample of 10000
sample_BA <- BA.matrix[sample(nrow(BA.matrix), 1000), ]
sample_BN <- BN.matrix[sample(nrow(BN.matrix), 1000), ]

#fit model 
fit.BA <- lm(gene_expression ~ gene_ID + pond , data = sample_BA)
fit.BN <- lm(gene_expression ~ gene_ID + pond, data = sample_BN)

an.res.BA <- anova(fit.BA)
an.res.BN <- anova(fit.BN)

an.res.BA
an.res.BN
