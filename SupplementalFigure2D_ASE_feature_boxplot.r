tiff(filename="SupplementalFigure2D.tiff", width=650)
layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow=TRUE))
color4 <- c("darkred", "pink", "lightblue", "darkblue")
boxplotdata <- list()
data1 <- read.table("high_AS_genes_seq_feature.xls", sep="\t", row.names=1, header=TRUE)
data1 <- as.vector(data1[,"exon_num"])
data1 <- data1[data1 != "---"]
data1 <- as.numeric(data1)
boxplotdata[["High"]] <- data1
data2 <- read.table("medium_AS_genes_seq_feature.xls", sep="\t", row.names=1, header=TRUE)
data2 <- as.vector(data2[,"exon_num"])
data2 <- data2[data2 != "---"]
data2 <- as.numeric(data2)
boxplotdata[["Medium"]] <- data2
data3 <- read.table("low_AS_genes_seq_feature.xls", sep="\t", row.names=1, header=TRUE)
data3 <- as.vector(data3[,"exon_num"])
data3 <- data3[data3 != "---"]
data3 <- as.numeric(data3)
boxplotdata[["Low"]] <- data3
data4 <- read.table("no_AS_genes_seq_feature.xls", sep="\t", row.names=1, header=TRUE)
data4 <- as.vector(data4[,"exon_num"])
data4 <- data4[data4 != "---"]
data4 <- as.numeric(data4)
boxplotdata[["No"]] <- data4
par(mar=c(3,5,2,2))
ylim1 <- min(c(data1,data2,data3,data4)) 
ylim1 <- 0
ylim2 <- max(c(data1,data2,data3,data4)) 
ylim2 <- 40
boxplot(boxplotdata,outline = FALSE, ylim=c(ylim1,ylim2), cex=1, cex.axis=2, pch=16, border=color4, las=1, boxwex = 0.3, frame.plot=FALSE, col=color4, ylab="Exon number", cex.lab=2.2)
result <- wilcox.test(data1, data2)$p.value
result
result <- wilcox.test(c(data1, data2), data3)$p.value
result
lines(c(1,1), c(29,32),col="black")
lines(c(2,2), c(29,32),col="black")
lines(c(1,2), c(32,32),col="black")
text(1.5, 31, "**",cex=2)

lines(c(1.5,1.5), c(32,35),col="black")
lines(c(1.5,3), c(35,35),col="black")
lines(c(3,3), c(30,35),col="black")
text(2.3, 33, "**",cex=2)
lines(c(2.3,2.3), c(35,40),col="black")
lines(c(4,4), c(35,40),col="black")
lines(c(2.3,4), c(40,40),col="black")
result <- wilcox.test(c(data1, data2, data3), data4)$p.value
result
text(3.3, 39, "**",cex=2)



boxplotdata <- list()
data1 <- read.table("high_AS_genes_seq_feature.xls", sep="\t", row.names=1, header=TRUE)
data1 <- as.vector(data1[,"intron_num"])
data1 <- data1[data1 != "---"]
data1 <- as.numeric(data1)
boxplotdata[["High"]] <- data1

data2 <- read.table("medium_AS_genes_seq_feature.xls", sep="\t", row.names=1, header=TRUE)
data2 <- as.vector(data2[,"intron_num"])
data2 <- data2[data2 != "---"]
data2 <- as.numeric(data2)
boxplotdata[["Medium"]] <- data2

data3 <- read.table("low_AS_genes_seq_feature.xls", sep="\t", row.names=1, header=TRUE)
data3 <- as.vector(data3[,"intron_num"])
data3 <- data3[data3 != "---"]
data3 <- as.numeric(data3)
boxplotdata[["Low"]] <- data3
data4 <- read.table("no_AS_genes_seq_feature.xls", sep="\t", row.names=1, header=TRUE)
data4 <- as.vector(data4[,"intron_num"])
data4 <- data4[data4 != "---"]
data4 <- as.numeric(data4)
boxplotdata[["No"]] <- data4
par(mar=c(3,5,2,2))
ylim1 <- min(c(data1,data2,data3,data4)) 
ylim2 <- max(c(data1,data2,data3,data4))
ylim2 <- 40
boxplot(boxplotdata,outline = FALSE, ylim=c(ylim1,ylim2), cex=1, cex.axis=2, pch=16, border=color4, boxwex = 0.3, frame.plot=FALSE, col=color4, ylab="Intron number", cex.lab=2.2)


result <- wilcox.test(data1, data2)$p.value
result
result <- wilcox.test(c(data1, data2), data3)$p.value
result
lines(c(1,1), c(28,31),col="black")
lines(c(2,2), c(28,31),col="black")
lines(c(1,2), c(31,31),col="black")
text(1.5, 30, "**",cex=2)

lines(c(1.5,1.5), c(31,35),col="black")
lines(c(1.5,3), c(35,35),col="black")
lines(c(3,3), c(30,35),col="black")
text(2.3, 33, "**",cex=2)
lines(c(2.3,2.3), c(35,40),col="black")
lines(c(4,4), c(35,40),col="black")
lines(c(2.3,4), c(40,40),col="black")
result <- wilcox.test(c(data1, data2, data3), data4)$p.value
result
text(3.3, 39, "**",cex=2)






boxplotdata <- list()
data1 <- read.table("high_AS_genes_seq_feature.xls", sep="\t", row.names=1, header=TRUE)
data1 <- as.vector(data1[,"exon_len"])
data1 <- data1[data1 != "---"]
data1 <- as.numeric(data1)
data1 <- data1/1000
boxplotdata[["High"]] <- data1

data2 <- read.table("medium_AS_genes_seq_feature.xls", sep="\t", row.names=1, header=TRUE)
data2 <- as.vector(data2[,"exon_len"])
data2 <- data2[data2 != "---"]
data2 <- as.numeric(data2)
data2 <- data2/1000
boxplotdata[["Medium"]] <- data2

data3 <- read.table("low_AS_genes_seq_feature.xls", sep="\t", row.names=1, header=TRUE)
data3 <- as.vector(data3[,"exon_len"])
data3 <- data3[data3 != "---"]
data3 <- as.numeric(data3)
data3 <- data3/1000
#data3 <- data.frame(Group=rep("Low AS", length(data3)), Value=data3)
boxplotdata[["Low"]] <- data3
data4 <- read.table("no_AS_genes_seq_feature.xls", sep="\t", row.names=1, header=TRUE)
data4 <- as.vector(data4[,"exon_len"])
data4 <- data4[data4 != "---"]
data4 <- as.numeric(data4)
data4 <- data4/1000
boxplotdata[["No"]] <- data4
par(mar=c(3,5,2,2))
ylim1 <- min(c(data1,data2,data3,data4)) 
ylim2 <- max(c(data1,data2,data3,data4))
ylim2 <- 10
boxplot(boxplotdata,outline = FALSE, ylim=c(ylim1,ylim2), cex=1, cex.axis=2, pch=16, border=color4, boxwex = 0.3, las=1, frame.plot=FALSE, col=color4, ylab="Exon length", cex.lab=2.2)


result <- wilcox.test(data1, data2)$p.value
result
result <- wilcox.test(c(data1, data2), data3)$p.value
result
lines(c(1,1), c(7.5,8),col="black")
lines(c(2,2), c(7.5,8),col="black")
lines(c(1,2), c(8,8),col="black")
#text(2.3, 7.8, "**",cex=2)
text(1.5, 7.5, "**",cex=2)

lines(c(1.5,1.5), c(8,9),col="black")
lines(c(1.5,3), c(9,9),col="black")
lines(c(3,3), c(8,9),col="black")

text(2.3, 8.5, "**",cex=2)


lines(c(2.3,2.3), c(9,10),col="black")
lines(c(4,4), c(7.5,10),col="black")
lines(c(2.3,4), c(10,10),col="black")
result <- wilcox.test(c(data1, data2, data3), data4)$p.value
result
text(3.3, 9.5, "**",cex=2)










boxplotdata <- list()
data1 <- read.table("high_AS_genes_seq_feature.xls", sep="\t", row.names=1, header=TRUE)
data1 <- as.vector(data1[,"intron_len"])
data1 <- data1[data1 != "---"]
data1 <- as.numeric(data1)
data1 <- data1/1000
boxplotdata[["High"]] <- data1

data2 <- read.table("medium_AS_genes_seq_feature.xls", sep="\t", row.names=1, header=TRUE)
data2 <- as.vector(data2[,"intron_len"])
data2 <- data2[data2 != "---"]
data2 <- as.numeric(data2)
data2 <- data2/1000
boxplotdata[["Medium"]] <- data2

data3 <- read.table("low_AS_genes_seq_feature.xls", sep="\t", row.names=1, header=TRUE)
data3 <- as.vector(data3[,"intron_len"])
data3 <- data3[data3 != "---"]
data3 <- as.numeric(data3)
data3 <- data3/1000
#data3 <- data.frame(Group=rep("Low AS", length(data3)), Value=data3)
boxplotdata[["Low"]] <- data3
data4 <- read.table("no_AS_genes_seq_feature.xls", sep="\t", row.names=1, header=TRUE)
data4 <- as.vector(data4[,"intron_len"])
data4 <- data4[data4 != "---"]
data4 <- as.numeric(data4)
data4 <- data4/1000
#data4 <- data.frame(Group=rep("No AS", length(data4)), Value=data4)
boxplotdata[["No"]] <- data4
par(mar=c(3,5,2,2))
ylim1 <- min(c(data1,data2,data3,data4)) 
ylim2 <- max(c(data1,data2,data3,data4))*0.6
ylim2 <- 15
boxplot(boxplotdata,outline = FALSE, ylim=c(ylim1,ylim2), cex=1, cex.axis=2, pch=16, border=color4, boxwex = 0.3, frame.plot=FALSE, col=color4, las=1, ylab="Intron length", cex.lab=2.2)


result <- wilcox.test(data1, data2)$p.value
result
result <- wilcox.test(c(data1, data2), data3)$p.value
result
lines(c(1,1), c(9.5,10),col="black")
lines(c(2,2), c(9.5,10),col="black")
lines(c(1,2), c(10,10),col="black")
text(1.5, 9.7, "**",cex=2)

lines(c(1.5,1.5), c(10,12),col="black")
lines(c(1.5,3), c(12,12),col="black")
lines(c(3,3), c(9.5,12),col="black")
text(2.3, 11.5, "**",cex=2)

lines(c(2.3,2.3), c(12,13),col="black")
lines(c(4,4), c(9.5,13),col="black")
lines(c(2.3,4), c(13,13),col="black")
text(3.3, 12.5, "**",cex=2)


result <- wilcox.test(c(data1, data2, data3), data4)$p.value
result 
#text(3.3, 11.9, "**",cex=2)





boxplotdata <- list()
data1 <- read.table("high_AS_genes_exprs.xls", sep="\t", row.names=1, header=TRUE)
data1 <- as.vector(data1[,"expression"])
data1 <- data1[data1 != "---"]
data1 <- as.numeric(data1)
boxplotdata[["High"]] <- data1

data2 <- read.table("medium_AS_genes_exprs.xls", sep="\t", row.names=1, header=TRUE)
data2 <- as.vector(data2[,"expression"])
data2 <- data2[data2 != "---"]
data2 <- as.numeric(data2)
boxplotdata[["Medium"]] <- data2

data3 <- read.table("low_AS_genes_exprs.xls", sep="\t", row.names=1, header=TRUE)
data3 <- as.vector(data3[,"expression"])
data3 <- data3[data3 != "---"]
data3 <- as.numeric(data3)
boxplotdata[["Low"]] <- data3
data4 <- read.table("no_AS_genes_exprs.xls", sep="\t", row.names=1, header=TRUE)
data4 <- as.vector(data4[,"expression"])
data4 <- data4[data4 != "---"]
data4 <- as.numeric(data4)
boxplotdata[["No"]] <- data4
par(mar=c(3,5,2,2))
ylim1 <- min(c(data1,data2,data3,data4)) 
ylim2 <- max(c(data1,data2,data3,data4))
ylim2 <- 15
boxplot(boxplotdata, cex=1, pch=16, cex.axis=2,outline = FALSE, ylim=c(ylim1,ylim2),boxwex = 0.3, frame.plot=FALSE, border=color4, col=color4, ylab="Log2(TPM+1)",las=1, cex.lab=2.2)

result <- wilcox.test(data1, data2)$p.value
result
result <- wilcox.test(c(data1, data2), data3)$p.value
result
lines(c(1,1), c(8,9),col="black")
lines(c(2,2), c(8,9),col="black")
lines(c(1,2), c(9,9),col="black")
text(1.5, 8.5, "**",cex=2)
lines(c(1.5,1.5), c(9,10),col="black")
lines(c(1.5,3), c(10,10),col="black")
lines(c(3,3), c(9,10),col="black")
text(2.3, 9.5, "**",cex=2)
lines(c(2.3,2.3), c(9.5,12),col="black")
lines(c(4,4), c(9.5,12),col="black")
lines(c(2.3,4), c(12,12),col="black")
result <- wilcox.test(c(data1, data2, data3), data4)$p.value
result
text(3.3, 11.5, "**",cex=2)




boxplotdata <- list()
data1 <- read.table("high_AS_genes_GC.xls", sep="\t", row.names=1, header=TRUE)
data1 <- as.vector(data1[,"GC"])
data1 <- data1[data1 != "---"]
data1 <- as.numeric(data1)
boxplotdata[["High"]] <- data1

data2 <- read.table("medium_AS_genes_GC.xls", sep="\t", row.names=1, header=TRUE)
data2 <- as.vector(data2[,"GC"])
data2 <- data2[data2 != "---"]
data2 <- as.numeric(data2)
boxplotdata[["Medium"]] <- data2

data3 <- read.table("low_AS_genes_GC.xls", sep="\t", row.names=1, header=TRUE)
data3 <- as.vector(data3[,"GC"])
data3 <- data3[data3 != "---"]
data3 <- as.numeric(data3)
#data3 <- data.frame(Group=rep("Low AS", length(data3)), Value=data3)
boxplotdata[["Low"]] <- data3

data4 <- read.table("no_AS_genes_GC.xls", sep="\t", row.names=1, header=TRUE)
data4 <- as.vector(data4[,"GC"])
data4 <- data4[data4 != "---"]
data4 <- as.numeric(data4)
#data4 <- data.frame(Group=rep("No AS", length(data4)), Value=data4)
boxplotdata[["No"]] <- data4
par(mar=c(3,5,2,2))
ylim1 <- min(c(data1,data2,data3,data4)) 
ylim2 <- max(c(data1,data2,data3,data4))
boxplot(boxplotdata, cex=1, cex.axis=2, pch=16, ylim=c(ylim1, ylim2), outline = FALSE,  boxwex = 0.3, border=color4, frame.plot=FALSE, col=color4, ylab="GC(%)", las=1, cex.lab=2.2)


result <- wilcox.test(data1, data2)$p.value
result
result <- wilcox.test(c(data1, data2), data3)$p.value
result
lines(c(1,1), c(55,60),col="black")
lines(c(2,2), c(55,60),col="black")
lines(c(1,2), c(60,60),col="black")
text(1.5, 58, "**",cex=2)

lines(c(1.5,1.5), c(60,65),col="black")
lines(c(1.5,3), c(65,65),col="black")
lines(c(3,3), c(55,65),col="black")
text(2.3, 64, "**",cex=2)

lines(c(2.3,2.3), c(65,72),col="black")
lines(c(4,4), c(70,72),col="black")
lines(c(2.3,4), c(72,72),col="black")
result <- wilcox.test(c(data1, data2, data3), data4)$p.value
text(3.3, 71, "**",cex=2)




dev.off()


