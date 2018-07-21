args<-commandArgs(T)
inputfile1 <- args[1]
d <- read.table(inputfile1, sep="\t", header=T, row.names=1)
d <- log2(d+1)
write.table(data.frame(geneid=rownames(d), d), file="result.xls", row.names=FALSE, sep="\t", quote=FALSE)

