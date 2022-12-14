#loading packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(grid)
library(BSDA)
library(magrittr)
library(ggridges)
library(scales)


#set directory
setwd("/Users/philipbaldassari/Desktop/Biol_5403_Genomics/Final_project/pi_per_site")


all <- read.csv("pi_r6_maf1.ann.vcf.sites.pi", sep='\t')

genes_mods <- read.csv("pi_all_genes_andmods_r6_maf1.ann.vcf.sites.pi", sep='\t')

genes <- read.csv("pi_all_genes_r6_maf1.ann.vcf.sites.pi", sep='\t')




den_all <- density(all$PI, adjust=2) %$%
  data.frame(x = x, y = y)

xall25 <- length(all$PI) * 0.75
xall10 <- length(all$PI) * 0.90
xall5 <- length(all$PI) * 0.95
xall1 <- length(all$PI) * 0.99

den25_all <- den_all %>%
  filter(x > sort(all$PI)[xall25])

den10_all <- den_all %>%
  filter(x > sort(all$PI)[xall10])

den5_all <- den_all %>%
  filter(x > sort(all$PI)[xall5])

den1_all <- den_all %>%
  filter(x > sort(all$PI)[xall1])


ggplot(den_all, aes(x=x, y=y)) + geom_line() + xlab("Pi on all SNPs") + ylab("Density") +
  layer(data=den25_all, geom="area", stat="identity", position_dodge(width=0), mapping=aes(x=x, y=y, fill="red", alpha=0.7)) +
  layer(data=den10_all, geom="area", stat="identity", position_dodge(width=0), mapping=aes(x=x, y=y, fill="red", alpha=0.7)) +
  layer(data=den5_all, geom="area", stat="identity", position_dodge(width=0), mapping=aes(x=x, y=y, fill="red", alpha=0.7)) +
  layer(data=den1_all, geom="area", stat="identity", position_dodge(width=0), mapping=aes(x=x, y=y, fill="red", alpha=0.7)) +
  theme_classic() + theme(legend.position = "none") +
  geom_vline(xintercept = c(sort(all$PI)[xall25], sort(all$PI)[xall10], sort(all$PI)[xall5], sort(all$PI)[xall1], max(den_all$x)), linetype=c("dotted", "dotted", "dotted", "dotted", "dashed"), color = c("indianred1", "indianred2", "indianred", "indianred3", "gray"), linewidth=0.6)





den_genes_mods <- density(genes_mods$PI, adjust=2) %$%
  data.frame(x = x, y = y)


xgenes_mods25 <- length(genes_mods$PI) * 0.75
xgenes_mods10 <- length(genes_mods$PI) * 0.90
xgenes_mods5 <- length(genes_mods$PI) * 0.95
xgenes_mods1 <- length(genes_mods$PI) * 0.99

den25_genes_mods <- den_genes_mods %>%
  filter(x > sort(genes_mods$PI)[xgenes_mods25])

den10_genes_mods <- den_genes_mods %>%
  filter(x > sort(genes_mods$PI)[xgenes_mods10])

den5_genes_mods <- den_genes_mods %>%
  filter(x > sort(genes_mods$PI)[xgenes_mods5])

den1_genes_mods <- den_genes_mods %>%
  filter(x > sort(genes_mods$PI)[xgenes_mods1])


ggplot(den_genes_mods, aes(x=x, y=y)) + geom_line() + xlab("Pi on SNPs in Genes and Upstream/Downstream Modifier Regions") + ylab("Density") +
  layer(data=den25_genes_mods, geom="area", stat="identity", position_dodge(width=0), mapping=aes(x=x, y=y, fill="red", alpha=0.7)) +
  layer(data=den10_genes_mods, geom="area", stat="identity", position_dodge(width=0), mapping=aes(x=x, y=y, fill="red", alpha=0.7)) +
  layer(data=den5_genes_mods, geom="area", stat="identity", position_dodge(width=0), mapping=aes(x=x, y=y, fill="red", alpha=0.7)) +
  layer(data=den1_genes_mods, geom="area", stat="identity", position_dodge(width=0), mapping=aes(x=x, y=y, fill="red", alpha=0.7)) +
  theme_classic() + theme(legend.position = "none") +
  geom_vline(xintercept = c(sort(genes_mods$PI)[xgenes_mods25], sort(genes_mods$PI)[xgenes_mods10], sort(genes_mods$PI)[xgenes_mods5], sort(genes_mods$PI)[xgenes_mods1], max(den_genes_mods$x)), linetype=c("dotted", "dotted", "dotted", "dotted", "dashed"), color = c("indianred1", "indianred2", "indianred", "indianred3", "gray"), linewidth=0.6)







den_genes <- density(genes$PI, adjust=2) %$%
  data.frame(x = x, y = y)

xgenes25 <- length(genes$PI) * 0.75
xgenes10 <- length(genes$PI) * 0.90
xgenes5 <- length(genes$PI) * 0.95
xgenes1 <- length(genes$PI) * 0.99

den25_genes <- den_genes %>%
  filter(x > sort(genes$PI)[xgenes25])

den10_genes <- den_genes %>%
  filter(x > sort(genes$PI)[xgenes10])

den5_genes <- den_genes %>%
  filter(x > sort(genes$PI)[xgenes5])

den1_genes <- den_genes %>%
  filter(x > sort(genes$PI)[xgenes1])


ggplot(den_genes, aes(x=x, y=y)) + geom_line() + xlab("Pi on SNPs in Genes") + ylab("Density") +
  layer(data=den25_genes, geom="area", stat="identity", position_dodge(width=0), mapping=aes(x=x, y=y, fill="red", alpha=0.7)) +
  layer(data=den10_genes, geom="area", stat="identity", position_dodge(width=0), mapping=aes(x=x, y=y, fill="red", alpha=0.7)) +
  layer(data=den5_genes, geom="area", stat="identity", position_dodge(width=0), mapping=aes(x=x, y=y, fill="red", alpha=0.7)) +
  layer(data=den1_genes, geom="area", stat="identity", position_dodge(width=0), mapping=aes(x=x, y=y, fill="red", alpha=0.7)) +
  theme_classic() + theme(legend.position = "none") +
  geom_vline(xintercept = c(sort(genes$PI)[xgenes25], sort(genes$PI)[xgenes10], sort(genes$PI)[xgenes5], sort(genes$PI)[xgenes1], max(den_genes$x)), linetype=c("dotted", "dotted", "dotted", "dotted", "dashed"), color = c("indianred1", "indianred2", "indianred", "indianred3", "gray"), linewidth=0.6)



############################################################################################################################################################################################################################################################


neuromods <- read.csv("pi_neurogenesis_andmods_r6_maf1.ann.vcf.sites.pi", sep='\t')
neuromods$GO <- "Neurogenesis"

nonmods <- read.csv("pi_nonneural_andmods_r6_maf1.ann.vcf.sites.pi", sep='\t')
nonmods$GO <- "Non-neural"

neurononenruomods <- rbind(neuromods, nonmods)
head(neurononenruomods)



wilcox.test(neurononenruomods$PI ~ neurononenruomods$GO, alternative = "greater")

v1 <- ggplot(neurononenruomods, aes(x=GO, y=PI)) + geom_violin(fill="cadetblue",  draw_quantiles=c(0.25,0.5,0.75)) + xlab("GO Class") + ylab("Pi") + ggtitle("p = 0.999") + theme_linedraw()
v1 + stat_summary(fun=mean, geom="point", size=3, color="black") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=12)) +
  scale_y_continuous(minor_breaks = seq(0, 0.55, 0.01))











neuro <- read.csv("pi_neurogenesis_r6_maf1.ann.vcf.sites.pi", sep='\t')
neuro$GO <- "Neurogenesis"

non <- read.csv("pi_nonneural_r6_maf1.ann.vcf.sites.pi", sep='\t')
non$GO <- "Non-neural"

neurononenruo <- rbind(neuro, non)
head(neurononenruo)



wilcox.test(neurononenruo$PI ~ neurononenruo$GO, alternative = "greater")

v2 <- ggplot(neurononenruo, aes(x=GO, y=PI)) + geom_violin(fill="cadetblue",  draw_quantiles=c(0.25,0.5,0.75)) + xlab("GO Class") + ylab("Pi") + ggtitle("p = 0.999") + theme_bw()
v2 + stat_summary(fun=mean, geom="point", size=3, color="black") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=12)) +
  scale_y_continuous(minor_breaks = seq(0, 0.55, 0.01))
  















neuroTFmods <- read.csv("pi_neurogenesis_TF_andmods_r6_maf1.ann.vcf.sites.pi", sep='\t')
neuroTFmods$GO <- "Neurogenesis TF"

nonTFmods <- read.csv("pi_nonneural_TF_andmods_r6_maf1.ann.vcf.sites.pi", sep='\t')
nonTFmods$GO <- "Non-neural TF"

neurononenruoTFmods <- rbind(neuroTFmods, nonTFmods)
head(neurononenruomods)



wilcox.test(neurononenruoTFmods$PI ~ neurononenruoTFmods$GO, alternative = "greater")

v3 <- ggplot(neurononenruoTFmods, aes(x=GO, y=PI)) + geom_violin(fill="cadetblue", draw_quantiles=c(0.25,0.5,0.75)) + xlab("GO Class") + ylab("Pi") + ggtitle("**") + theme_linedraw()
v3 + stat_summary(fun=mean, geom="point", size=3, color="black") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30)) +
  scale_y_continuous(minor_breaks = seq(0, 0.55, 0.01))



bv3 <- ggplot(neurononenruoTFmods, aes(x=GO, y=PI)) + xlab("GO Class") + ylab("Pi") + ggtitle("**") + theme_linedraw()
bv3 + stat_summary(fun=mean, geom="bar", fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30)) +
  geom_violin(aes(x=GO, y=PI), fill = NA, size=0.5, width = 0.5, draw_quantiles=c(0.25,0.5,0.75)) +
  scale_y_continuous(minor_breaks = seq(0, 0.55, 0.01))



bv4 <- ggplot(neurononenruoTFmods, aes(x=GO, y=PI)) + xlab("GO Class") + ylab("Pi") + ggtitle("**") + theme_linedraw()
bv4 + stat_summary(fun=mean, geom="bar", fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30)) +
  scale_y_continuous(minor_breaks = seq(0, 0.2, 0.005))








neuroTF <- read.csv("pi_neurogenesis_TF_r6_maf1.ann.vcf.sites.pi", sep='\t')
neuroTF$GO <- "Neurogenesis TF"

nonTF <- read.csv("pi_nonneural_TF_r6_maf1.ann.vcf.sites.pi", sep='\t')
nonTF$GO <- "Non-neural TF"

neurononenruoTF <- rbind(neuroTF, nonTF)
head(neurononenruoTF)



wilcox.test(neurononenruoTF$PI ~ neurononenruoTF$GO, alternative = "greater")

v4 <- ggplot(neurononenruoTF, aes(x=GO, y=PI)) + geom_violin(fill="cadetblue", draw_quantiles=c(0.25,0.5,0.75)) + xlab("GO Class") + ylab("Pi") + ggtitle("***") + theme_bw()
v4 + stat_summary(fun=mean, geom="point", size=3, color="black") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30)) +
  scale_y_continuous(minor_breaks = seq(0, 0.55, 0.01))




bv5 <- ggplot(neurononenruoTF, aes(x=GO, y=PI)) + xlab("GO Class") + ylab("Pi") + ggtitle("***") + theme_linedraw()
bv5 + stat_summary(fun=mean, geom="bar", fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30)) +
  geom_violin(aes(x=GO, y=PI), fill = NA, size=0.5, width = 0.5, draw_quantiles=c(0.25,0.5,0.75)) +
  scale_y_continuous(minor_breaks = seq(0, 0.55, 0.01))



bv6 <- ggplot(neurononenruoTF, aes(x=GO, y=PI)) + xlab("GO Class") + ylab("Pi") + ggtitle("***") + theme_linedraw()
bv6 + stat_summary(fun=mean, geom="bar", fill="cadetblue") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30)) +
  scale_y_continuous(minor_breaks = seq(0, 0.2, 0.005))





############################################################################################################################################################################################################################################################



neuroSS <- read.csv("pi_neurogenesis_SS_r6_maf1.ann.vcf.sites.pi", sep='\t')
neuroSS$GO <- "Neurogenesis"

nonSS <- read.csv("pi_nonneural_SS_r6_maf1.ann.vcf.sites.pi", sep='\t')
nonSS$GO <- "Non-neural"

neurononeuroSS <- rbind(neuroSS, nonSS)
head(neurononeuroSS)



wilcox.test(neurononeuroSS$PI ~ neurononeuroSS$GO, alternative = "greater")

v5 <- ggplot(neurononeuroSS, aes(x=GO, y=PI)) + geom_violin(fill="cadetblue", draw_quantiles=c(0.25,0.5,0.75)) + xlab("GO Class") + ylab("Ps") + ggtitle("p = 0.9942") + theme_linedraw() 
v5 + stat_summary(fun=mean, geom="point", size=3, color="black") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=12)) +
  scale_y_continuous(minor_breaks = seq(0, 0.55, 0.01))
  



logv5 <- ggplot(neurononeuroSS, aes(x=GO, y=PI)) + geom_violin(fill="cadetblue", draw_quantiles=c(0.25,0.5,0.75)) + xlab("GO Class") + ylab("log(Ps)") + ggtitle("p = 0.9942") + theme_bw() + scale_y_continuous(trans='log10')
logv5 + stat_summary(fun=mean, geom="point", size=3, color="black") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=12))
  






reproductionSS <- read.csv("pi_reproduction_SS_r6_maf1.ann.vcf.sites.pi", sep='\t')
reproductionSS$GO <- "Reproduction"

immuneSS <- read.csv("pi_immune_SS_r6_maf1.ann.vcf.sites.pi", sep='\t')
immuneSS$GO <- "Immune"


allSS <- rbind(neurononeuroSS, reproductionSS, immuneSS)
head(allSS)



ggplot(allSS, aes(x=PI, y=reorder(GO, PI, FUN = mean))) + ylab("GO Class") + xlab("Ps") + theme_bw() +
  geom_density_ridges(scale = 2, alpha = 0.7, fill="cadetblue") + 
  stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.25,0.5,0.75), scale = 2, alpha = 0.7, fill="cadetblue") +
  stat_summary(fun=mean, geom="point", size=3, color="red")




ggplot(allSS, aes(x=PI, y=reorder(GO, PI, FUN = mean))) + ylab("GO Class") + xlab("log(Ps)") + theme_bw() + scale_x_log10(oob = scales::squish_infinite) +
  geom_density_ridges(scale = 2, alpha = 0.7, fill="cadetblue") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.25,0.5,0.75), scale = 2, alpha = 0.7, fill="cadetblue") +
  stat_summary(fun=mean, geom="point", size=3, color="red")




ggplot(allSS, aes(x=reorder(GO, -PI, FUN = mean), y=PI)) + geom_violin(fill="cadetblue", draw_quantiles=c(0.25,0.5,0.75)) + xlab("GO Class") + ylab("Ps") + theme_bw() +
  stat_summary(fun=mean, geom="point", size=2, color="red")




############################################################################################################################################################################################################################################################


pnps <- read.csv("PnPs_final_all_genes_effects.csv")

wilcox.test(pnps$PnPs ~ pnps$Neuro, alternative = "less")
ggplot(pnps, aes(x=Neuro, y=PnPs)) + geom_violin(fill="cadetblue") + xlab("GO Class") + ylab("Pn/Ps") + ggtitle("***") + theme_bw() +
  stat_summary(fun=mean, geom="point", size=2) + theme(axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30)) +
  geom_hline(yintercept=1, linetype="dotted", color="gray", linewidth=1)



ggplot(pnps, aes(x=Neuro, y=PnPs)) + geom_boxplot(fill="cadetblue") + xlab("GO Class") + ylab("Pn/Ps (log scale)") + ggtitle("***") + theme_bw() + scale_y_log10(oob = scales::squish_infinite) +
  stat_summary(fun=mean, geom="point", size=2, color="red") + theme(axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30)) +
  geom_hline(yintercept=1, linetype="dotted", color="gray", linewidth=1)



ggplot(pnps, aes(x=Neuro, y=PnPs)) + xlab("GO Class") + ylab("Mean Pn/Ps") + ggtitle("***") + theme_bw() +
  stat_summary(fun=mean, geom="bar", fill="cadetblue") + theme(axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30))




wilcox.test(pnps$Num_of_SNPs ~ pnps$Neuro, alternative = "greater")
ggplot(pnps, aes(x=Neuro, y=Num_of_SNPs)) + geom_violin(fill="cadetblue", draw_quantiles=c(0.25,0.5,0.75)) + xlab("GO Class") + ylab("CDS SNPs per Gene (log scale)") + ggtitle("***") + theme_bw() + scale_y_log10(oob = scales::squish_infinite) +
  stat_summary(fun=mean, geom="point", size=2) + theme(axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30))


ggplot(pnps, aes(x=Neuro, y=Num_of_SNPs)) + xlab("GO Class") + ylab("Mean CDS SNPs per Gene") + ggtitle("***") + theme_bw() +
  stat_summary(fun=mean, geom="bar", fill="cadetblue") + theme(axis.title = element_text(size=15), plot.title = element_text(hjust = 0.5, size=30))





ggplot(pnps, aes(x=PnPs, y=reorder(GO, PnPs, FUN = mean))) + ylab("GO Class") + xlab("log(Pn/Ps)") + theme_bw()  + scale_x_log10(oob = scales::squish_infinite) +
  geom_density_ridges(scale = 1, alpha = 0.7, fill="cadetblue") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.25,0.5,0.75), scale = 1, alpha = 0.7, fill="cadetblue") +
  stat_summary(fun=mean, geom="point", size=2, color="red") +
  geom_vline(xintercept=1, linetype="dotted", color="gray", linewidth=1)


ggplot(pnps, aes(x=PnPs, y=reorder(GO, PnPs, FUN = mean))) + geom_boxplot(fill="cadetblue") + ylab("GO Class") + xlab("Pn/Ps (log scale)") + theme_bw() + scale_x_log10(oob = scales::squish_infinite) +
  stat_summary(fun=mean, geom="point", size=2, color="red") +
  geom_vline(xintercept=1, linetype="dotted", color="gray", linewidth=1)









































