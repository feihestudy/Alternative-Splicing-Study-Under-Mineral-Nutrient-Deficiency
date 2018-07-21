argsx <- commandArgs(TRUE)
#pdf(file="ASG_condition.pdf", width=10)
ttintron <- argsx[1]
ttintron <- as.numeric(ttintron)
pdf("FigS1_ASG_condition.pdf")
library(gplots)
data <-read.table("ASG_condition.xls", sep="\t", header=FALSE)
data <- data[,2]
per <- data/ttintron
per <- paste("Total:", data, "(",round(per*100,2), "%", ")",sep="")
#per <- paste("", "",round(per*100,2), "%", "",sep="")
par(mar=c(10,8,2,5))
ymax <- max(data)*1.5
bp <- barplot2(data,plot.grid=FALSE, col="darkblue", xaxt="n",cex.axis=1.3, las=2, names.arg=c("Control", "-Fe", "-Zn", "-Cu", "-Mn", "+Pi Root", "-Pi Root", "+Pi Shoot", "-Pi Shoot", "Root","Shoot","Total"), beside=TRUE, ylim=c(0, ymax),  cex.names=1.2, horiz = F, main="", cex.lab=1.5, cex.main=1.8, ylab="")
axis(1, pos=0, at=bp[,1],labels=c("Control", "-Fe", "-Zn", "-Cu", "-Mn", "+Pi Root", "-Pi Root", "+Pi Shoot", "-Pi Shoot", "Root","Shoot","Total"),cex.axis=1.6, las=2)
mtext("Number of detected AS genes",side=2,line=5, cex=1.6)
cexv <- 1.5
#text(bp[1,1], data[1]+1500, per[1], cex=cexv)
#text(bp[2,1], data[2]+1900, per[2], cex=cexv)
#text(bp[3,1], data[3]+2800, per[3], cex=cexv)
#text(bp[4,1], data[4]+3500, per[4], cex=cexv)
#text(bp[5,1], data[5]+2800, per[5], cex=cexv)
#text(bp[6,1], data[6]+2000, per[6], cex=cexv)
#text(bp[7,1], data[7]+3500, per[7], cex=cexv)
#text(bp[8,1], data[8]+2800, per[8], cex=cexv)
#text(bp[9,1], data[9]+3500, per[9], cex=cexv)
#text(bp[10,1], data[10]+3200, per[10], cex=cexv)
#text(bp[11,1], data[11]+3800, per[11], cex=cexv)
#text(bp[12,1], data[12]+4500, per[12], cex=cexv)
text(bp[6,1], data[6]+4500, per[12], cex=cexv)
#lines(c(bp[12,1], bp[12,1]), c(data[12], data[12]+4500), lwd=2)
dev.off()

