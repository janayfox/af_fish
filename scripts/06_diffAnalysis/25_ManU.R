#to use on cluster if analysis cant run on my computer

library(rstatix)

plastic <- read.csv("./plastManU.csv")
evol <- read.csv("./evolManU.csv")

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
