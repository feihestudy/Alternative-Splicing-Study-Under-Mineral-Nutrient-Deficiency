argsx <- commandArgs(TRUE)
inputfile <- argsx[1]
intron_list_expressed <- argsx[2]
data <- read.table(inputfile, sep="\t", row.names=1, header=TRUE)
sums <- data[,1]
numcutoff_high <- mean(sums)+sd(sums)
#numcutoff <- mean(sums)
numcutoff_high <- round(numcutoff_high, 0)
numcutoff_high

numcutoff_low <- mean(sums)-sd(sums)
#numcutoff <- mean(sums)
numcutoff_low <- round(numcutoff_low, 0)
numcutoff_low

sum(sums<=numcutoff_low)
sum(sums>=numcutoff_high)
sum(sums>numcutoff_low & sums<numcutoff_high)

highASgene <- rownames(data)[sums>=numcutoff_high]
write.table(highASgene, sep="\t", file="high_AS_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

mediumASgene <- rownames(data)[sums>numcutoff_low & sums<numcutoff_high]
write.table(mediumASgene, sep="\t", file="medium_AS_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

lowASgene <- rownames(data)[sums<=numcutoff_low]
write.table(lowASgene, sep="\t", file="low_AS_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


noas <- as.vector(read.table(intron_list_expressed, sep="\t",header=FALSE)[,1])
noas <- noas[!noas %in% rownames(data)]
length(noas)
write.table(noas, sep="\t", file="no_AS_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

