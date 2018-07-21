ASEdata <- read.table("PiFeZnCuMn_ASE_type.xls", sep="\t", header=TRUE)
RIdata <- as.vector(ASEdata[1,2])
A3SSdata <- as.vector(ASEdata[2,2])
A5SSdata <- as.vector(ASEdata[3,2])
SEdata <- as.vector(ASEdata[4,2])
MXEdata <- as.vector(ASEdata[5,2])
tiff(filename="Figure1A.tiff", width=500)
plot(NA, NA,xlim=c(1, 12), ylim=c(0, 7),type="n",col="grey", main="", cex.main=0.9, ylab="", xlab="", frame.plot=FALSE, yaxt="n",xaxt="n")
i=5
rect(1, i-0.1, 2, i+0.1, col="black")
rect(2, i-0.1, 5, i+0.1, col="white")
rect(5, i-0.1, 6, i+0.1, col="black")
lines(c(2,3.5),c(i-0.1,i-0.3),col = "black")
lines(c(3.5,5),c(i-0.3,i-0.1),col = "black")
#### A3SS
i=4
rect(1, i-0.1, 2, i+0.1, col="black")
rect(5, i-0.1, 6, i+0.1, col="black")
rect(4, i-0.1, 5, i+0.1, col="white")
lines(c(2,3),c(i,i+0.3),col = "black")
lines(c(3,4),c(i+0.3,i),col = "black")
lines(c(2,3.5),c(i,i-0.3),col = "black")
lines(c(5,3.5),c(i-0.1,i-0.3),col = "black")
#### A5SS
i=3
rect(1, i-0.1, 2, i+0.1, col="black")
rect(2, i-0.1, 3, i+0.1, col="white")
rect(5, i-0.1, 6, i+0.1, col="black")
lines(c(2+0.1,3.5),c(i+0.1,i+0.3),col = "black")
lines(c(3.5,5),c(i+0.3,i),col = "black")
lines(c(3,4),c(i,i-0.3),col = "black")
lines(c(5,4),c(i,i-0.3),col = "black")
#### SE
i=2
rect(1, i-0.1, 2, i+0.1, col="black")
rect(3, i-0.1, 4, i+0.1, col="white")
rect(5, i-0.1, 6, i+0.1, col="black")
lines(c(2,2.5),c(i,i+0.3),col = "black")
lines(c(2.5,3),c(i+0.3,i),col = "black")
lines(c(4,4.5),c(i,i+0.3),col = "black")
lines(c(4.5,5),c(i+0.3,i),col = "black")
#lines(c(3,4),c(i,i-0.3),col = "black")
lines(c(2,3.5),c(i,i-0.3),col = "black")
lines(c(5,3.5),c(i,i-0.3),col = "black")
#### MXE
i=1
rect(1, i-0.1, 2, i+0.1, col="black")
rect(2.5, i-0.1, 3, i+0.1, col="white")
rect(4, i-0.1, 4.5, i+0.1, col="white")
rect(5, i-0.1, 6, i+0.1, col="black")
lines(c(2,2.25),c(i,i+0.3),col = "black")
lines(c(2.25,2.5),c(i+0.3,i),col = "black")
lines(c(3,4),c(i,i+0.3),col = "black")
lines(c(4,5),c(i+0.3,i),col = "black")
lines(c(4.5,4.75),c(i,i-0.3),col = "black")
lines(c(4.75,5),c(i-0.3,i),col = "black")
lines(c(2,3),c(i,i-0.3),col = "black")
lines(c(3,4),c(i-0.3,i),col = "black")
text(7, 5.05, "RI", cex=1.5)
text(10, 5.05, RIdata, cex=1.5)
text(7.3, 4.05, "A3SS", cex=1.5)
text(10, 4.05, A3SSdata, cex=1.5)
text(7.3, 3.05, "A5SS", cex=1.5)
text(10, 3.05, A5SSdata, cex=1.5)
text(7, 2.05, "SE", cex=1.5)
text(10, 2.05, SEdata, cex=1.5)
text(7.2, 1.05, "MXE", cex=1.5)
text(10, 1.05, MXEdata, cex=1.5)
text(7.6, 6.05, "Event type", cex=1.5)
text(10.6, 6.05, "Percentage(%)", cex=1.5)
text(2, 6.05, "Event plot", cex=1.5)
lines(c(1,12),c(5.5,5.5),col = "black", lwd=2)
lines(c(1,12),c(0.5,0.5),col = "black", lwd=2)
lines(c(1,12),c(6.5,6.5),col = "black", lwd=2)
dev.off()
