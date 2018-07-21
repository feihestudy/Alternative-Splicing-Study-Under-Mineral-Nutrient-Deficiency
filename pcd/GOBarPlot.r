argsx <- commandArgs(TRUE)
inputfile <- argsx[1]
library(plyr)
library(ggplot2)
library(scales)
library(gplots)
data <- read.table(inputfile,sep="\t", header=TRUE, row.names=1)
data_deg <- data[data$DEG >= -log10(0.05), ]
data_deg <- data_deg[order(data_deg$DEG, decreasing=FALSE),]
data_dasg <- data[!rownames(data) %in% rownames(data_deg),]
data_dasg <- data_dasg[order(data_dasg$DASG, decreasing=FALSE),]
data <- rbind(data_deg, data_dasg)
dd <- as.matrix(data)
ymin <- min(dd)
ymax <- max(dd)
dd <- t(dd)
pdf(file="Plot.pdf")
par(mar=c(5, 22, 2, 2))
bp <- barplot2(dd,plot.grid=FALSE, col=c("blue3", "red"), cex.axis=1.5, names.arg=colnames(dd), beside=TRUE, las=2, cex.names=1, horiz = TRUE, xlim=c(ymin, ymax*1.5), main="", cex.main=1.8, xlab="-log10(Pvalue)", cex.lab=1.6)
legend("topright",cex=1.5, legend=c('DEG','DASG'),pch=c(15,15),box.col="NA", col=c("blue3", "red"))
dev.off()
