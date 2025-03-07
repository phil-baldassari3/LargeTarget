library(tidyverse)
library(ggplot2)

setwd("~/Desktop/neuro/GO_coverage/output")


df <- read.csv("GO_coverage.csv")

ggplot(df, aes(x = reorder(GO_term, -CDS_coverage), y = CDS_coverage)) +
  geom_bar(stat="identity", aes(fill = ifelse(GO_term == "Neural", "Neural", "Other"))) +
  scale_fill_manual(values = c("Neural" = "cadetblue", "Other" = "gray")) +
  xlab("GO Class") + ylab("Total Coverage (bp)") + ggtitle("CDS Coverage") +
  scale_y_continuous(labels=scales::comma) + theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), legend.position = "none")










############################################################################################################################################


csv <- "neural_vs_non_gene_CDS_lengths.csv"
title <- "Average CDS Length"


#opening data
df <- read.csv(csv)
df$Ontology <- as.factor(df$Ontology)

#filtering
neuro <- df %>%
  filter(Ontology == "Neural")

nonneuro <- df %>%
  filter(Ontology == "Non-neural")



###
neurols <- neuro$Gene_length
nonneurols <- nonneuro$Gene_length
###
neurols <- neuro$CDS_length
nonneurols <- nonneuro$CDS_length
###



###Statistics
hist(neurols)
hist(nonneurols)

#t-test
t.test(neurols, nonneurols)

#wilcox
wilcox.test(neurols, nonneurols)




#plotting
ggplot(df, aes(x = Ontology, y = Gene_length, fill = Ontology)) +
  stat_summary(fun = "mean", geom = "bar", width = 0.7) + 
  stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), geom = "errorbar", width = 0.2) +
  scale_fill_manual(values = c("Neural" = "cadetblue", "Non-neural" = "gray")) +
  xlab("Tissue Type") + ylab("Mean Length (bp)") + ggtitle(title) +
  theme_bw() + theme(legend.position = "none")










