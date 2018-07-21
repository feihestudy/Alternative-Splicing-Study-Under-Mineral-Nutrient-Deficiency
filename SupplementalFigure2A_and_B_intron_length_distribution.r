#### all gene expression
gene_feature <- read.table("gene_intron_distribution1.xls", sep="\t", header=TRUE)
intron_number <- gene_feature$intron_number
intron_number <- intron_number[intron_number!=0]
a <- sum(intron_number < 10)
b <- length(intron_number)
stat1 <- c("Total_intron_genes", b, b)
stat2 <- c("Few_than_10", a, round(a/b, 3))
a <- sum(intron_number > 1)
stat3 <- c("More_than_1", a, round(a/b, 3))
#write.table(rbind(stat1,stat2,stat3), file="fewer_than_10_percentage.xls", sep="\t", col.names=FALSE, row.names=FALSE)
intron_number <- table(intron_number)
#write.table(intron_number, file="intron_vs_gene_number.xls", sep="\t", row.names=FALSE)
intron_number <- intron_number[1:30]
gene_feature <- gene_feature[gene_feature[,"intron_number"] > 0,]
#write.table(data.frame(gene_id=rownames(gene_feature), gene_feature), row.names=FALSE, file="gene_intron_containing_distribution.xls", sep="\t")
pdf(file="SupplementalFigure2A.pdf", width=9)
library(gplots)
par(mar=c(6,8,2,2))
barplot2(intron_number, yaxt="n",col="red", cex.lab=1.8, space=0.5,cex.axis=2,cex.names=1.5,las=2, ylim=c(0, max(intron_number)*1.3), ylab="Number of genes", xlab="Number of introns", names.arg=names(intron_number), las=2)
axis(2, seq(0,10000, by=1000),paste(seq(0, 10, by=1),"K",sep=""),las=2, cex.axis=1.5)
gene_feature <- read.table("gene_intron_distribution2.xls", sep="\t", header=TRUE)
gene_feature <- gene_feature[gene_feature$intron_length != 0,]
intron_size_dis <- NULL
for(i in 1:500)
{
   intron_size_dis <- c(intron_size_dis, length(unique(as.vector(gene_feature[gene_feature$intron_length == i,"gene_id"]))))
}
dev.off()
pdf(file="SupplementalFigure2B.pdf")
par(mar=c(6,8,2,2))
plot(1:500, intron_size_dis, type="l",ylab="Number of introns", cex.axis=1.4, xlab="Intron size (bp)", col="red",cex.lab=1.8)
intron_length <- gene_feature$intron_length
intron_length <- intron_length[intron_length!=0]
md <- median(intron_length)
mn <- mean(intron_length)
mn <- round(mn, 1)
tb <- table(intron_length)
maxnum <- as.numeric(names(tb)[tb==max(tb)])
ttintron <- length(intron_length)
text(300, 1200, paste("Mean intron size: ", mn, sep=""), cex=1.5)
text(300, 1100, paste("Median intron size: ", md, sep=""), cex=1.5)
text(300, 1000, paste("Mode intron size:      ", maxnum, sep=""), cex=1.5)
text(300, 900, paste("Total of introns: ", ttintron, sep=""), cex=1.5)
dev.off()

