library(gplots)
data <- read.table("high_AS_genes_domain_enrichment_significant.xls",sep="\t", header=TRUE)
data <- data[order(data$classic_fisher),]
data <- data[1:20,]
data2 <- data[,c("Significant", "Expected")]
data2 <- as.matrix(data2)
data2 <- t(data2)
pv_log <- -log10(data$classic_fisher)
pv_log[is.infinite(pv_log)] <- 25
cls <- rep("red", length(pv_log))
goids <- as.vector(data$Domain_id)
pname <- as.vector(data$Domain_name)
pname <- gsub("'","",pname)
one <- pname[goids == "IPR000504"]
one <- sub("\\.\\s+\\(.*\\)", "", one)
pname[goids == "IPR000504"] <- one
one <- pname[goids == "IPR001680"]
one <- sub("\\,\\s+G\\-beta\\s+repeat", "", one)
pname[goids == "IPR001680"] <- one
one <- pname[goids == "IPR001327"]
pname[goids == "IPR001327"] <- paste("oxidoreductase", "(IPR001327)", sep="")
one <- pname[goids == "IPR023753"]
pname[goids == "IPR023753"] <- paste("oxidoreductase", "(IPR023753)", sep="")
one <- pname[goids == "IPR003594"]
pname[goids == "IPR003594"] <- paste("ATPase", "(IPR003594)", sep="")
one <- pname[goids == "IPR003959"]
pname[goids == "IPR003959"] <- paste("ATPase", "(IPR003959)", sep="")
data3 <- data2
cls <- c(rep(c("red4", "red"), length(pname)))
jpeg(filename="Figure4A.jpg",width=500,height=500, quality = 100)
par(xpd=T, mar=par()$mar+c(17,3,0,8))
bp <- barplot2(data3, beside=TRUE, col=cls, cex.names=1.2, cex.axis=1.5,las=2, ylim=c(0, max(data3)*1.2), names.arg=pname, las=2)
legend1 <- NULL
legend2 <- NULL
legend3 <- NULL
legend1 <- c(legend1, c('Observed','Expected'))
legend2 <- c(legend2, c(15, 15))
legend3 <- c(legend3, c("red4","red"))
legend("topright",cex=1.4,legend=legend1,pch=legend2,box.col="NA", col=legend3)
mtext("Number of genes",side=2,line=3, cex=1.5)
par(new=TRUE)
pvs <- c(pv_log)
cls <- c(rep("red", length(pname)))
plot(bp[1,], pvs, type="o", col=NA, axes=F, xlab="",cex.axis=1.2, ylim=c(0, max(pvs)*1.2),ylab="", pch=16, cex=0.8)
points(bp[1,][grep("red", cls)], pvs[grep("red", cls)], type="o", col="blue", pch=16, cex=0.6)
axis(4, cex.axis=1.5, las=2)
mtext("-log10(Pvalue)",side=4,line=3, cex=1.5)
dev.off()
