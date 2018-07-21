data <- read.table("PiFeZnCuMn_psi.xls", sep="\t", header=TRUE, row.names=1)
data2 <- (1-data)/data
write.table(data.frame(ASE=rownames(data2), data2), file="PiFeZnCuMn_splicing_efficiency.xls", sep="\t", row.names=FALSE, quote=FALSE)

