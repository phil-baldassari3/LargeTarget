library(ggplot2)

setwd("/Users/philipbaldassari/Desktop/neuro/TFBS")

tfs <- read.csv("ChIP_seq_TF.csv")

tfbytissue <- read.csv("ChIP_seq_TF_by_cell_type.csv")


head(tfs)
head(tfbytissue)


ggplot(tfs, aes(x=reorder(Cell_type, -Binding_site_length), y=Binding_site_length)) + geom_boxplot(fill="cadetblue") + 
  xlab("Cell Type") + ylab("chIP-seq peak length") + theme_bw() + scale_y_log10(oob = scales::squish_infinite) +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15))


plot(hist(tfs$Binding_site_length, breaks=100), xlab="length (bp)", main="Distribution of chIP-seq TF Peak lengths")




p1 <- ggplot(tfbytissue, aes(x = reorder(Cell_Type, -Number_of_binding_sites), y = Number_of_binding_sites)) + geom_bar(stat="identity", fill="cadetblue") + xlab("Cell Type") + ylab("Number of Peaks (binding sites)") + scale_y_continuous(labels=scales::comma) + theme_bw()
p1 + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15))



p2 <- ggplot(tfbytissue, aes(x = reorder(Cell_Type, -Total_length_of_binding_sites), y = Total_length_of_binding_sites)) + geom_bar(stat="identity", fill="cadetblue") + xlab("Cell Type") + ylab("Total Length of Peaks (binding sites)") + scale_y_continuous(labels=scales::comma) + theme_bw()
p2 + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, size=12), axis.title = element_text(size=15))


