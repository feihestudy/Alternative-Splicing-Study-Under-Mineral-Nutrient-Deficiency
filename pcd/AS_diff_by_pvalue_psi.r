data <- read.table("AS.xls", sep="\t", header=TRUE, row.names=1)
n <- length(colnames(data))
pvaluename <- colnames(data)[n]
fdrs <- p.adjust(data[,n], method ="fdr")
result <- data.frame(ASE_id=rownames(data), data, FDR=fdrs)
n <- length(colnames(result))
colnames(result)[n] <- sub("Pvalue","FDR",pvaluename)
write.table(result, file="AS_FDR.xls", sep="\t", row.names=FALSE, quote=FALSE)
result <- result[result[,n] <= 0.05 & abs(result[,(n-2)]) >= 0.1,]
write.table(result, file="AS_diff.xls", sep="\t", row.names=FALSE, quote=FALSE)

