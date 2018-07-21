args<-commandArgs(T)
inputfile1 <- args[1]
inputfile2 <- args[2]
d <- read.table(inputfile1, sep="\t", header=T, row.names=1)
sample2con <- read.table(inputfile2, sep="\t", header=FALSE)
samplenames <- as.vector(sample2con[,1])
connames <- as.vector(sample2con[,2])
d <- d[,samplenames]
colnames(d) <- connames
write.table(data.frame(geneid=rownames(d), d), file="result.xls", row.names=FALSE, sep="\t", quote=FALSE)

