library(tidyverse)
library(ggplot2)

setwd("~/Desktop/neuro/output/gene_lengths")

csv <- "tissue_specific_lncRNAs_neural_vs_non.csv"    #SET THIS TO THE FILE YOU WISH TO PLOT FROM
title <- ""


#opening data
df <- read.csv(csv)
df$Tissue_Type <- as.factor(df$Tissue_Type)

#filtering
neuro <- df %>%
  filter(Tissue_Type == "Neural")

nonneuro <- df %>%
  filter(Tissue_Type == "Non-neural")

neurols <- neuro$Length
nonneurols <- nonneuro$Length

###Statistics
hist(neurols)
hist(nonneurols)

#t-test
t.test(neurols, nonneurols)

#wilcox
wilcox.test(neurols, nonneurols)




#plotting
ggplot(df, aes(x = Tissue_Type, y = Length, fill = Tissue_Type)) +
  stat_summary(fun = "mean", geom = "bar", width = 0.7) + 
  stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), geom = "errorbar", width = 0.2) +
  scale_fill_manual(values = c("Neural" = "skyblue3", "Non-neural" = "gray")) +
  xlab("Tissue Type") + ylab("Mean Length (bp)") + ggtitle(title) +
  theme_bw() + theme(legend.position = "none")










