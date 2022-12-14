#loading packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(grid)
library(BSDA)


#set directory
setwd("/Users/philipbaldassari/Desktop/Biol_5403_Genomics/Final_project/site_counting")


gene_lens <- read.csv("GO_genes/gene_lengths_GO.csv")

all_gene_lens <- read.csv("GO_genes/all_total_avg_gene_lengths.csv")
head(gene_lens)

comined_gene_lens <- read.csv("GO_genes/combine_total_avg_gene_lengths.csv")

neuralnonneural_gene_lens <- read.csv("GO_genes/neuralnonneural_total_avg_gene_lengths.csv")



p1 <- ggplot(all_gene_lens, aes(x = reorder(Class, -Total), y = Total)) + geom_bar(stat="identity", fill="cadetblue") + xlab("GO Class") + ylab("Total Coverage (bp)") + scale_y_continuous(labels=scales::comma) + theme_bw()
p1 + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15))


p2 <- ggplot(comined_gene_lens, aes(x = reorder(Class, -Total), y = Total)) + geom_bar(stat="identity", fill="cadetblue") + xlab("GO Class") + ylab("Total Coverage (bp)") + scale_y_continuous(labels=scales::comma) + theme_bw()
p2 + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15))




test_results <- wilcox.test(gene_lens$Neurogenesis, gene_lens$Nonneural, alternative = "greater")
test_results





p3 <- ggplot(comined_gene_lens, aes(x = reorder(Class, -Mean), y = Mean)) + geom_bar(stat="identity", fill="cadetblue") + xlab("GO Class") + ylab("Mean Gene Size (bp)") + scale_y_continuous(labels=scales::comma) + theme_bw()
p3 + theme(axis.text.x = element_text(angle=90, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15))


p4 <- ggplot(neuralnonneural_gene_lens, aes(x = Class, y = Mean)) + geom_bar(stat="identity", fill="cadetblue") + xlab("GO Class") + ylab("Mean Gene Size (bp)") + scale_y_continuous(labels=scales::comma) + ggtitle("***") + theme_bw()
p4 + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30))





############################################################################################################################################################################################################################################################




v1 <- ggplot(exon_GOclasses, aes(x=reorder(GO, -Length, FUN = mean), y=Length)) + geom_violin(fill="cadetblue") + xlab("GO Class") + ylab("Length of Exons (bp)") + theme_bw()
v1 + stat_summary(fun.y=mean, geom="point", size=2, color="red") + theme(axis.text.x = element_text(angle=90, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15))



v2 <- ggplot(intro_GOclasses, aes(x=reorder(GO, -Length, FUN = mean), y=Length)) + geom_violin(fill="cadetblue") + xlab("GO Class") + ylab("Length of Introns (bp)") + theme_bw()
v2 + stat_summary(fun.y=mean, geom="point", size=2, color="red") + theme(axis.text.x = element_text(angle=90, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15))



v3 <- ggplot(UTR5_GOclasses, aes(x=reorder(GO, -Length, FUN = mean), y=Length)) + geom_violin(fill="cadetblue") + xlab("GO Class") + ylab("Length of 5' UTRs (bp)") + theme_bw()
v3 + stat_summary(fun.y=mean, geom="point", size=2, color="red") + theme(axis.text.x = element_text(angle=90, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15))



v4 <- ggplot(UTR3_GOclasses, aes(x=reorder(GO, -Length, FUN = mean), y=Length)) + geom_violin(fill="cadetblue") + xlab("GO Class") + ylab("Length of 3' UTRs (bp)") + theme_bw()
v4 + stat_summary(fun.y=mean, geom="point", size=2, color="red") + theme(axis.text.x = element_text(angle=90, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15))








wilcox.test(exon_neuro$Length, exon_non$Length, alternative = "greater")

b1 <- ggplot(exons, aes(x=reorder(Neuro, -Length, FUN = mean), y=Length)) + xlab("GO Class") + ylab("Mean Length of Exons (bp)") + ggtitle("p = 0.4879") + theme_bw()
b1 + stat_summary(fun.y=mean, geom="bar", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=12))


vb1 <- ggplot(exons, aes(x=reorder(Neuro, -Length, FUN = mean), y=Length)) + geom_violin(fill="cadetblue") + xlab("GO Class") + ylab("Length of Exons (bp)") + ggtitle("p = 0.4879") + theme_bw()
vb1 + stat_summary(fun.y=mean, geom="point", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=12))



wilcox.test(intro_neuro$Length, intro_non$Length, alternative = "greater")

b2 <- ggplot(introns, aes(x=reorder(Neuro, -Length, FUN = mean), y=Length)) + xlab("GO Class") + ylab("Mean Length of Introns (bp)") + ggtitle("***") + theme_bw()
b2 + stat_summary(fun.y=mean, geom="bar", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30))


vb2 <- ggplot(introns, aes(x=reorder(Neuro, -Length, FUN = mean), y=Length)) + geom_violin(fill="cadetblue") + xlab("GO Class") + ylab("Length of Introns (bp)") + ggtitle("***") + theme_bw()
vb2 + stat_summary(fun.y=mean, geom="point", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30))



wilcox.test(UTR5_neuro$Length, UTR5_non$Length, alternative = "greater")

b3 <- ggplot(UTR5, aes(x=reorder(Neuro, -Length, FUN = mean), y=Length)) + xlab("GO Class") + ylab("Mean Length of 5' UTRs (bp)") + ggtitle("***") + theme_bw()
b3 + stat_summary(fun.y=mean, geom="bar", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30))


vb3 <- ggplot(UTR5, aes(x=reorder(Neuro, -Length, FUN = mean), y=Length)) + geom_violin(fill="cadetblue") + xlab("GO Class") + ylab("Length of 5' UTRs (bp)") + ggtitle("***") + theme_bw()
vb3 + stat_summary(fun.y=mean, geom="point", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30))



wilcox.test(UTR3_neuro$Length, UTR3_non$Length, alternative = "greater")

b4 <- ggplot(UTR3, aes(x=reorder(Neuro, -Length, FUN = mean), y=Length)) + xlab("GO Class") + ylab("Mean Length of 3' UTRs (bp)") + ggtitle("***") + theme_bw()
b4 + stat_summary(fun.y=mean, geom="bar", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30))


vb4 <- ggplot(UTR3, aes(x=reorder(Neuro, -Length, FUN = mean), y=Length)) + geom_violin(fill="cadetblue") + xlab("GO Class") + ylab("Length of 3' UTRs (bp)") + ggtitle("***") + theme_bw()
vb4 + stat_summary(fun.y=mean, geom="point", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30))



############################################################################################################################################################################################################################################################


TF_per_gene <- read.csv("GO_TFBSs/TF_counts_neurononneuro.csv")

TF_length <- read.csv("GO_TFBSs/TF_lengths_neurononneuro.csv")


all_TFBS_length <- read.csv("GO_TFBSs/all.csv")

neurononneuro_TFBS_length <- read.csv("GO_TFBSs/neurononneuro.csv")

neurononneuro_TFBS_counts <- read.csv("GO_TFBSs/TFBS_counts_neurononneuro.csv")




wilcox.test(TF_per_gene$Num_of_TF ~ TF_per_gene$GO, alternative = "greater")

x1 <- ggplot(TF_per_gene, aes(x=GO, y=Num_of_TF)) + xlab("GO Class") + ylab("Mean Number of TF per gene") + ggtitle("p = 0.1196") + theme_bw()
x1 + stat_summary(fun.y=mean, geom="bar", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=12))

vx1 <- ggplot(TF_per_gene, aes(x=GO, y=Num_of_TF)) + geom_violin(fill="cadetblue") + xlab("GO Class") + ylab("Number of TF per gene") + ggtitle("p = 0.1196") + theme_bw()
vx1 + stat_summary(fun.y=mean, geom="point", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=12))


wilcox.test(TF_length$factor_length ~ TF_length$GO, alternative = "greater")

x2 <- ggplot(TF_length, aes(x=GO, y=factor_length)) + xlab("GO Class") + ylab("Mean Length of TF (bp)") + ggtitle("p = 0.7434") + theme_bw()
x2 + stat_summary(fun.y=mean, geom="bar", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=12))


vx2 <- ggplot(TF_length, aes(x=GO, y=factor_length)) + geom_violin(fill="cadetblue") + xlab("GO Class") + ylab("Length of TF (bp)") + ggtitle("p = 0.7434") + theme_bw()
vx2 + stat_summary(fun.y=mean, geom="point", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=12))



x3 <- ggplot(all_TFBS_length, aes(x=reorder(GO, -TFBS_length, FUN = mean), y=TFBS_length)) + geom_violin(fill="cadetblue") + xlab("GO Class") + ylab("Length of TFBS (bp)") + theme_bw()
x3 + stat_summary(fun.y=mean, geom="point", size=2, color="red") + theme(axis.text.x = element_text(angle=90, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15))


wilcox.test(neurononneuro_TFBS_length$TFBS_length ~ neurononneuro_TFBS_length$GO, alternative = "less")

x4 <- ggplot(neurononneuro_TFBS_length, aes(x=reorder(GO, -TFBS_length, FUN = mean), y=TFBS_length)) + xlab("GO Class") + ylab("Mean Length of TFBS (bp)") + ggtitle("***") + theme_bw()
x4 + stat_summary(fun.y=mean, geom="bar", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30))

vx4 <- ggplot(neurononneuro_TFBS_length, aes(x=reorder(GO, -TFBS_length, FUN = mean), y=TFBS_length)) + geom_violin(fill="cadetblue") + xlab("GO Class") + ylab("Length of TFBS (bp)") + ggtitle("***") + theme_bw()
vx4 + stat_summary(fun.y=mean, geom="point", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30))



wilcox.test(neurononneuro_TFBS_counts$Num_of_TFBS ~ neurononneuro_TFBS_counts$GO, alternative = "greater")

x5 <- ggplot(neurononneuro_TFBS_counts, aes(x=reorder(GO, -Num_of_TFBS, FUN = mean), y=Num_of_TFBS)) + xlab("GO Class") + ylab("Mean Number of TFBS per Gene") + ggtitle("p = 0.1172") + theme_bw()
x5 + stat_summary(fun.y=mean, geom="bar", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=12))

vx5 <- ggplot(neurononneuro_TFBS_counts, aes(x=reorder(GO, -Num_of_TFBS, FUN = mean), y=Num_of_TFBS)) + geom_violin(fill="cadetblue") + xlab("GO Class") + ylab("Number of TFBS per Gene") + ggtitle("p = 0.1172") + theme_bw()
vx5 + stat_summary(fun.y=mean, geom="point", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=12))







CRM_per_gene <- read.csv("GO_CRMs/CRM_counts_neurononneuro.csv")

CRM_length <- read.csv("GO_CRMs/neurononneuro.csv")

all_CRM_length <- read.csv("GO_CRMs/all.csv")



wilcox.test(CRM_per_gene$Num_of_CRM ~ CRM_per_gene$GO, alternative = "greater")

a1 <- ggplot(CRM_per_gene, aes(x=GO, y=Num_of_CRM)) + xlab("GO Class") + ylab("Mean Number of CRMs per gene") + ggtitle("***") + theme_bw()
a1 + stat_summary(fun.y=mean, geom="bar", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30))

va1 <- ggplot(CRM_per_gene, aes(x=GO, y=Num_of_CRM)) + geom_violin(fill="cadetblue") + xlab("GO Class") + ylab("Number of CRMs per gene") + ggtitle("***") + theme_bw()
va1 + stat_summary(fun.y=mean, geom="point", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30))



wilcox.test(CRM_length$length ~ CRM_length$GO, alternative = "greater")

a2 <- ggplot(CRM_length, aes(x=GO, y=length)) + xlab("GO Class") + ylab("Mean Length of CRMs (bp)") + ggtitle("***") + theme_bw()
a2 + stat_summary(fun.y=mean, geom="bar", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30))

va2 <- ggplot(CRM_length, aes(x=GO, y=length)) + geom_violin(fill="cadetblue") + xlab("GO Class") + ylab("Length of CRMs (bp)") + ggtitle("***") + theme_bw()
va2 + stat_summary(fun.y=mean, geom="point", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30))



a3 <- ggplot(all_CRM_length, aes(x=reorder(GO, -length, FUN = mean), y=length)) + geom_violin(fill="cadetblue") + xlab("GO Class") + ylab("Length of CRMs (bp)") + theme_bw()
a3 + stat_summary(fun.y=mean, geom="point", size=2, color="red") + theme(axis.text.x = element_text(angle=90, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15))





############################################################################################################################################################################################################################################################


NS <- read.csv("CDS4site_counting/ann_potential_NS_site_counting.csv")

NS_all  <- NS %>%
  filter(GO != "None")

NS_neuro <- NS %>%
  filter(Neuro == "Neurogenesis")

NS_non <- NS %>%
  filter(Neuro == "Non-neural")


wilcox.test(NS$Number_NS_sites ~ NS$Neuro, alternative = "greater")

c1 <- ggplot(NS, aes(x=Neuro, y=Number_NS_sites)) + xlab("GO Class") + ylab("Mean Number of Potential NS Sites") + ggtitle("***") + theme_bw()
c1 + stat_summary(fun.y=mean, geom="bar", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30))



z.test(NS_neuro$NS_sites_per_bp, NS_non$NS_sites_per_bp, alternative = "greater", sigma.x=sd(NS_neuro$NS_sites_per_bp), sigma.y=sd(NS_non$NS_sites_per_bp),conf.level=.95)

c2 <- ggplot(NS, aes(x=Neuro, y=NS_sites_per_bp)) + xlab("GO Class") + ylab("Mean Number of Potential NS Sites per bp") + ggtitle("p = 0.9971") + theme_bw()
c2 + stat_summary(fun.y=mean, geom="bar", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=12))



c3 <- ggplot(NS, aes(x=Neuro, y=NS_sites_per_bp)) + geom_violin(fill="cadetblue") + xlab("GO Class") + ylab("Number of Potential NS Sites per bp") + ggtitle("p = 0.9971") + theme_bw()
c3 + stat_summary(fun.y=mean, geom="point", size=2, fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=12))



c4 <- ggplot(NS_all, aes(x=reorder(GO, -NS_sites_per_bp, FUN = mean), y=NS_sites_per_bp)) + geom_violin(fill="cadetblue") + xlab("GO Class") + ylab("Number of Potential NS Sites per bp") + theme_bw()
c4 + stat_summary(fun.y=mean, geom="point", size=2, color="red") + theme(axis.text.x = element_text(angle=90, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15))










