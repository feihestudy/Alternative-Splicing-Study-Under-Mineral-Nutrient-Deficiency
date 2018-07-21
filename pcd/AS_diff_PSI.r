data <- read.table("AS.xls", sep="\t", header=TRUE, row.names=1)
index <- grep("_PSI", colnames(data))
baseline <- data[,index[1:3]]
mean2 <- function(one)
{
   mean(one, na.rm=TRUE)
}
baseline <- apply(baseline, 1, mean2)
data[,index[4]] <- data[,index[4]]-baseline
data[,index[5]] <- data[,index[5]]-baseline
data[,index[6]] <- data[,index[6]]-baseline
t123 <- data.frame(ASE_id=rownames(data), data[,index[4:6]])
write.table(t123, file="AS_diff_filter_PSI.xls", sep="\t", row.names=FALSE, quote=FALSE)

