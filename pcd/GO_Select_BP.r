argsx <- commandArgs(TRUE)
inputfile <- argsx[1]

outputfile <- inputfile

outputfile <- sub("\\.xls", "", outputfile)

outputfile <- paste(outputfile, "_BP.xls", sep="")

dt <- read.table(inputfile,sep="\t", header=TRUE)

dt <- dt[as.vector(dt$GO_category) == "biological_process",]

dt <- dt[order(dt$Pvalue_by_Hypergeometric_test),]

write.table(dt, outputfile, quote=FALSE, sep="\t", row.names=FALSE)



