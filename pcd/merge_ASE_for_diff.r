argsx <- commandArgs(TRUE)
paste2 <- function(one)
{
    re <- paste(one, collapse=",")
	return(re)
}
g1s <- strsplit(argsx[1], ",")[[1]]
g2s <- strsplit(argsx[2], ",")[[1]]
IC_SAMPLE_1 <- NULL
SC_SAMPLE_1 <- NULL
IC_SAMPLE_2 <- NULL
SC_SAMPLE_2 <- NULL
ID <- NULL
for(i in 1:length(g1s))
{
    one <- read.table(paste("../ASE/", g1s[i], "/JCEC.RNASeq.RI.MATS.input.txt", sep=""), sep="\t", header=TRUE)
    IC_SAMPLE_1 <- cbind(IC_SAMPLE_1, one[,2])
	SC_SAMPLE_1 <- cbind(SC_SAMPLE_1, one[,3])
    ID <- one[,1]
	IncFormLen <- one[,6]
	SkipFormLen <- one[,7]
}
#IC_SAMPLE_1 <- paste(IC_SAMPLE_1[,1:2], sep=",")
IC_SAMPLE_1 <- apply(IC_SAMPLE_1, 1, paste2)
SC_SAMPLE_1 <- apply(SC_SAMPLE_1, 1, paste2)
for(i in 1:length(g2s))
{
    one <- read.table(paste("../ASE/", g2s[i], "/JCEC.RNASeq.RI.MATS.input.txt", sep=""), sep="\t", header=TRUE)
    IC_SAMPLE_2 <- cbind(IC_SAMPLE_2, one[,2])
	SC_SAMPLE_2 <- cbind(SC_SAMPLE_2, one[,3])

}
#IC_SAMPLE_1 <- paste(IC_SAMPLE_1[,1:2], sep=",")
IC_SAMPLE_2 <- apply(IC_SAMPLE_2, 1, paste2)
SC_SAMPLE_2 <- apply(SC_SAMPLE_2, 1, paste2)
result <- data.frame(ID=ID, IC_SAMPLE_1=IC_SAMPLE_1, SC_SAMPLE_1=SC_SAMPLE_1, IC_SAMPLE_2=IC_SAMPLE_2, SC_SAMPLE_2=SC_SAMPLE_2, IncFormLen=IncFormLen, SkipFormLen=SkipFormLen)
write.table(result, file="JCEC.RNASeq.RI.MATS.input.txt", sep="\t", row.names=FALSE, quote=FALSE)






####A3SS
IC_SAMPLE_1 <- NULL
SC_SAMPLE_1 <- NULL
IC_SAMPLE_2 <- NULL
SC_SAMPLE_2 <- NULL
ID <- NULL
for(i in 1:length(g1s))
{
    one <- read.table(paste("../ASE/", g1s[i], "/JCEC.RNASeq.A3SS.MATS.input.txt", sep=""), sep="\t", header=TRUE)
    IC_SAMPLE_1 <- cbind(IC_SAMPLE_1, one[,2])
	SC_SAMPLE_1 <- cbind(SC_SAMPLE_1, one[,3])
    ID <- one[,1]
	IncFormLen <- one[,6]
	SkipFormLen <- one[,7]
}
#IC_SAMPLE_1 <- paste(IC_SAMPLE_1[,1:2], sep=",")
IC_SAMPLE_1 <- apply(IC_SAMPLE_1, 1, paste2)
SC_SAMPLE_1 <- apply(SC_SAMPLE_1, 1, paste2)
for(i in 1:length(g2s))
{
    one <- read.table(paste("../ASE/", g2s[i], "/JCEC.RNASeq.A3SS.MATS.input.txt", sep=""), sep="\t", header=TRUE)
    IC_SAMPLE_2 <- cbind(IC_SAMPLE_2, one[,2])
	SC_SAMPLE_2 <- cbind(SC_SAMPLE_2, one[,3])

}
#IC_SAMPLE_1 <- paste(IC_SAMPLE_1[,1:2], sep=",")
IC_SAMPLE_2 <- apply(IC_SAMPLE_2, 1, paste2)
SC_SAMPLE_2 <- apply(SC_SAMPLE_2, 1, paste2)
result <- data.frame(ID=ID, IC_SAMPLE_1=IC_SAMPLE_1, SC_SAMPLE_1=SC_SAMPLE_1, IC_SAMPLE_2=IC_SAMPLE_2, SC_SAMPLE_2=SC_SAMPLE_2, IncFormLen=IncFormLen, SkipFormLen=SkipFormLen)
write.table(result, file="JCEC.RNASeq.A3SS.MATS.input.txt", sep="\t", row.names=FALSE, quote=FALSE)





####A5SS
IC_SAMPLE_1 <- NULL
SC_SAMPLE_1 <- NULL
IC_SAMPLE_2 <- NULL
SC_SAMPLE_2 <- NULL
ID <- NULL
for(i in 1:length(g1s))
{
    one <- read.table(paste("../ASE/", g1s[i], "/JCEC.RNASeq.A5SS.MATS.input.txt", sep=""), sep="\t", header=TRUE)
    IC_SAMPLE_1 <- cbind(IC_SAMPLE_1, one[,2])
	SC_SAMPLE_1 <- cbind(SC_SAMPLE_1, one[,3])
    ID <- one[,1]
	IncFormLen <- one[,6]
	SkipFormLen <- one[,7]
}
#IC_SAMPLE_1 <- paste(IC_SAMPLE_1[,1:2], sep=",")
IC_SAMPLE_1 <- apply(IC_SAMPLE_1, 1, paste2)
SC_SAMPLE_1 <- apply(SC_SAMPLE_1, 1, paste2)
for(i in 1:length(g2s))
{
    one <- read.table(paste("../ASE/", g2s[i], "/JCEC.RNASeq.A5SS.MATS.input.txt", sep=""), sep="\t", header=TRUE)
    IC_SAMPLE_2 <- cbind(IC_SAMPLE_2, one[,2])
	SC_SAMPLE_2 <- cbind(SC_SAMPLE_2, one[,3])

}
#IC_SAMPLE_1 <- paste(IC_SAMPLE_1[,1:2], sep=",")
IC_SAMPLE_2 <- apply(IC_SAMPLE_2, 1, paste2)
SC_SAMPLE_2 <- apply(SC_SAMPLE_2, 1, paste2)
result <- data.frame(ID=ID, IC_SAMPLE_1=IC_SAMPLE_1, SC_SAMPLE_1=SC_SAMPLE_1, IC_SAMPLE_2=IC_SAMPLE_2, SC_SAMPLE_2=SC_SAMPLE_2, IncFormLen=IncFormLen, SkipFormLen=SkipFormLen)
write.table(result, file="JCEC.RNASeq.A5SS.MATS.input.txt", sep="\t", row.names=FALSE, quote=FALSE)






####SE
IC_SAMPLE_1 <- NULL
SC_SAMPLE_1 <- NULL
IC_SAMPLE_2 <- NULL
SC_SAMPLE_2 <- NULL
ID <- NULL
for(i in 1:length(g1s))
{
    one <- read.table(paste("../ASE/", g1s[i], "/JCEC.RNASeq.SE.MATS.input.txt", sep=""), sep="\t", header=TRUE)
    IC_SAMPLE_1 <- cbind(IC_SAMPLE_1, one[,2])
	SC_SAMPLE_1 <- cbind(SC_SAMPLE_1, one[,3])
    ID <- one[,1]
	IncFormLen <- one[,6]
	SkipFormLen <- one[,7]
}
#IC_SAMPLE_1 <- paste(IC_SAMPLE_1[,1:2], sep=",")
IC_SAMPLE_1 <- apply(IC_SAMPLE_1, 1, paste2)
SC_SAMPLE_1 <- apply(SC_SAMPLE_1, 1, paste2)
for(i in 1:length(g2s))
{
    one <- read.table(paste("../ASE/", g2s[i], "/JCEC.RNASeq.SE.MATS.input.txt", sep=""), sep="\t", header=TRUE)
    IC_SAMPLE_2 <- cbind(IC_SAMPLE_2, one[,2])
	SC_SAMPLE_2 <- cbind(SC_SAMPLE_2, one[,3])

}
#IC_SAMPLE_1 <- paste(IC_SAMPLE_1[,1:2], sep=",")
IC_SAMPLE_2 <- apply(IC_SAMPLE_2, 1, paste2)
SC_SAMPLE_2 <- apply(SC_SAMPLE_2, 1, paste2)
result <- data.frame(ID=ID, IC_SAMPLE_1=IC_SAMPLE_1, SC_SAMPLE_1=SC_SAMPLE_1, IC_SAMPLE_2=IC_SAMPLE_2, SC_SAMPLE_2=SC_SAMPLE_2, IncFormLen=IncFormLen, SkipFormLen=SkipFormLen)
write.table(result, file="JCEC.RNASeq.SE.MATS.input.txt", sep="\t", row.names=FALSE, quote=FALSE)






####MXE
IC_SAMPLE_1 <- NULL
SC_SAMPLE_1 <- NULL
IC_SAMPLE_2 <- NULL
SC_SAMPLE_2 <- NULL
ID <- NULL
for(i in 1:length(g1s))
{
    one <- read.table(paste("../ASE/", g1s[i], "/JCEC.RNASeq.MXE.MATS.input.txt", sep=""), sep="\t", header=TRUE)
    IC_SAMPLE_1 <- cbind(IC_SAMPLE_1, one[,2])
	SC_SAMPLE_1 <- cbind(SC_SAMPLE_1, one[,3])
    ID <- one[,1]
	IncFormLen <- one[,6]
	SkipFormLen <- one[,7]
}
#IC_SAMPLE_1 <- paste(IC_SAMPLE_1[,1:2], sep=",")
IC_SAMPLE_1 <- apply(IC_SAMPLE_1, 1, paste2)
SC_SAMPLE_1 <- apply(SC_SAMPLE_1, 1, paste2)
for(i in 1:length(g2s))
{
    one <- read.table(paste("../ASE/", g2s[i], "/JCEC.RNASeq.MXE.MATS.input.txt", sep=""), sep="\t", header=TRUE)
    IC_SAMPLE_2 <- cbind(IC_SAMPLE_2, one[,2])
	SC_SAMPLE_2 <- cbind(SC_SAMPLE_2, one[,3])

}
#IC_SAMPLE_1 <- paste(IC_SAMPLE_1[,1:2], sep=",")
IC_SAMPLE_2 <- apply(IC_SAMPLE_2, 1, paste2)
SC_SAMPLE_2 <- apply(SC_SAMPLE_2, 1, paste2)
result <- data.frame(ID=ID, IC_SAMPLE_1=IC_SAMPLE_1, SC_SAMPLE_1=SC_SAMPLE_1, IC_SAMPLE_2=IC_SAMPLE_2, SC_SAMPLE_2=SC_SAMPLE_2, IncFormLen=IncFormLen, SkipFormLen=SkipFormLen)
write.table(result, file="JCEC.RNASeq.MXE.MATS.input.txt", sep="\t", row.names=FALSE, quote=FALSE)

