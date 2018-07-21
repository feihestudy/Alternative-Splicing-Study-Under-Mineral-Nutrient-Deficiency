argsx <- commandArgs(TRUE)
ASE <- argsx[1]
pcd <- argsx[2]
samplelist <- strsplit(argsx[3], ",")[[1]]
fun1 <- paste(pcd, "/plot_gene_exon_structure_AS_Fun5.r", sep="")
tmp <- strsplit(ASE, "\\&")[[1]]
target_st <- NULL
target_en <- NULL
upstream_st <- NULL
upstream_en <- NULL
downstream_st <- NULL
downstream_en <- NULL
strand <- tmp[4]


if(tmp[1] == "SE")
{
target_st <- tmp[5]
target_en <- tmp[6]
target_st <- as.numeric(target_st)
target_en <- as.numeric(target_en)
upstream_st <- as.numeric(tmp[7])
upstream_en <- as.numeric(tmp[8])
target_st <- as.numeric(target_st)
target_en <- as.numeric(target_en)
downstream_st <- as.numeric(tmp[9])
downstream_en <- as.numeric(tmp[10])
}

if(tmp[1] == "RI")
{
target_st <- tmp[8]
target_en <- tmp[9]
target_st <- as.numeric(target_st)
target_en <- as.numeric(target_en)
}

if(tmp[1] == "A3SS" && strand == "+")
{
	target_st <- tmp[5]
	target_en <- tmp[7]
	upstream_st <- tmp[9]
	upstream_en <- tmp[10]
	downstream_st <- tmp[7]
	downstream_en <- tmp[8]	
}
if(tmp[1] == "A3SS" && strand == "-")
{
	target_st <- tmp[8]
	target_en <- tmp[6]
	upstream_st <- tmp[7]
	upstream_en <- tmp[8]
	downstream_st <- tmp[9]
	downstream_en <- tmp[10]
}


if(tmp[1] == "A5SS" && strand == "+")
{
	target_st <- tmp[8]
	target_en <- tmp[6]
	upstream_st <- tmp[7]
	upstream_en <- tmp[8]
	downstream_st <- tmp[9]
	downstream_en <- tmp[10]	
}
if(tmp[1] == "A5SS" && strand == "-")
{
	target_st <- tmp[5]
	target_en <- tmp[7]
	upstream_st <- tmp[9]
	upstream_en <- tmp[10]
	downstream_st <- tmp[7]
	downstream_en <- tmp[8]
}


target_st <- as.numeric(target_st)
target_en <- as.numeric(target_en)
upstream_st <- as.numeric(upstream_st)
upstream_en <- as.numeric(upstream_en)
downstream_st <- as.numeric(downstream_st)
downstream_en <- as.numeric(downstream_en)


n <- length(tmp)
gst <- min(as.numeric(tmp[5:n]))
gen <- max(as.numeric(tmp[5:n]))

otfile <- paste(ASE, ".tiff", sep="")
tiff(filename=otfile)

##########################  1


colorss <- c("lightblue", "red")
layout(matrix(c(1,2), 2, 1, byrow=TRUE), height=c(2.5,1))
par(mar=c(2,13,0,3))

plot(NA, NA,xlim=c(gst, gen), ylim=c(0, 2.5),type="n",col="grey", main="", cex.main=0.9, ylab="", xlab="Pos(bp)", frame.plot=FALSE, yaxt="n",xaxt="n")

distance <- seq(gst, gen, by=(gen-gst)/2)

theta <- NULL
sss_n <- length(samplelist)






	#Control1
	Control <- read.table("Control.txt", sep="\t", header=FALSE)
	n <- dim(Control)[1]
	depth <- Control[,3]
	depth <- (depth-min(depth))/(max(depth)-min(depth)+0.1) 
	for(i in 1:n)
	{
		lines(c(Control[i,2], Control[i,2]), c(1.2,depth[i]+1.2), col="lightblue")
	}
	#filename <- paste(samplelist[2], ".txt", sep="")
	#Control1
	Case <- read.table("Case.txt", sep="\t", header=FALSE)
	n <- dim(Case)[1]
	depth <- Case[,3]
	depth <- (depth-min(depth))/(max(depth)-min(depth)+0.1) 
	for(i in 1:n)
	{
		lines(c(Case[i,2], Case[i,2]), c(0.2,depth[i]+0.2), col="red")
	}



if(tmp[1] == "A3SS")
{
	#filename <- paste(samplelist[1], "_junction.txt", sep="")
	Control <- read.table("Control_junction.txt", sep="\t", header=FALSE)
	Control <- as.vector(Control[,1])
	tmparr <- strsplit(Control, ",")[[1]]
	tmparr <- as.numeric(tmparr)
	tmparr[1] <- round(tmparr[1], 0)
	tmparr[2] <- round(tmparr[2], 0)
	#theta <- c(theta, tmparr[3])
	if(strand == "+")
	{
	    ii <- 1.2
		mid <- (upstream_en+target_st)/2
		#lines(c(target_st, target_en), c(3.3, 3.3), lty=3)
		lines(c(upstream_en, mid), c(ii, ii+0.1), lwd=1.5)
		lines(c(mid, target_st), c(ii+0.1, ii), lwd=1.5)
		#points(upstream_en, ii, pch=10, cex=0.5)
		#points(target_st, ii, pch=10, cex=0.5)
		text(mid, ii+0.15, tmparr[1], cex=2)
		mid <- (upstream_en+downstream_st)/2
		lines(c(upstream_en, mid), c(ii, ii-0.1), lwd=1.5)
		lines(c(mid, downstream_st), c(ii-0.1, ii), lwd=1.5)
		#lines(c(mid, target_en), c(ii-0.3, ii), lwd=1)
		#points(upstream_en, ii, pch=10, cex=0.5)
		#points(downstream_st, ii, pch=10, cex=0.5)
		text(mid, ii-0.15, tmparr[2], cex=2)
	}
	if(strand == "-")
	{
	    ii <- 1.2
		mid <- (target_en+downstream_st)/2
		#lines(c(target_st, target_en), c(3.3, 3.3), lty=3)
		lines(c(target_en, mid), c(ii, ii+0.1), lwd=1.5)
		lines(c(mid, downstream_st), c(ii+0.1, ii), lwd=1.5)
		#points(target_en, ii, pch=10, cex=0.5)
		#points(downstream_st, ii, pch=10, cex=0.5)
		text(mid, ii+0.15, tmparr[1], cex=2)
		mid <- (upstream_en+downstream_st)/2
		lines(c(upstream_en, mid), c(ii, ii-0.1), lwd=1.5)
		lines(c(mid, downstream_st), c(ii-0.1, ii), lwd=1.5)
		#lines(c(mid, target_en), c(ii-0.3, ii), lwd=1)
		#points(upstream_en, ii, pch=10, cex=0.5)
		#points(downstream_st, ii, pch=10, cex=0.5)
		text(mid, ii-0.15, tmparr[2], cex=2)
	}	

	
	
	#filename <- paste(samplelist[2], "_junction.txt", sep="")
	Case <- read.table("Case_junction.txt", sep="\t", header=FALSE)
	Case <- as.vector(Case[,1])
	tmparr <- strsplit(Case, ",")[[1]]
	tmparr <- as.numeric(tmparr)
	tmparr[1] <- round(tmparr[1], 0)
	tmparr[2] <- round(tmparr[2], 0)
	theta <- c(theta, tmparr[3])
	if(strand == "+")
	{
	    ii <- 0.2
		mid <- (upstream_en+target_st)/2
		#lines(c(target_st, target_en), c(3.3, 3.3), lty=3)
		lines(c(upstream_en, mid), c(ii, ii+0.1), lwd=1.5)
		lines(c(mid, target_st), c(ii+0.1, ii), lwd=1.5)
		#points(upstream_en, ii, pch=10, cex=0.5)
		#points(target_st, ii, pch=10, cex=0.5)
		text(mid, ii+0.15, tmparr[1], cex=2)
		mid <- (upstream_en+downstream_st)/2
		lines(c(upstream_en, mid), c(ii, ii-0.1), lwd=1.5)
		lines(c(mid, downstream_st), c(ii-0.1, ii), lwd=1.5)
		#lines(c(mid, target_en), c(ii-0.3, ii), lwd=1)
		#points(upstream_en, ii, pch=10, cex=0.5)
		#points(downstream_st, ii, pch=10, cex=0.5)
		text(mid, ii-0.15, tmparr[2], cex=2)
	}
	if(strand == "-")
	{
	    ii <- 0.2
		mid <- (target_en+downstream_st)/2
		#lines(c(target_st, target_en), c(3.3, 3.3), lty=3)
		lines(c(target_en, mid), c(ii, ii+0.1), lwd=1.5)
		lines(c(mid, downstream_st), c(ii+0.1, ii), lwd=1.5)
		#points(target_en, ii, pch=10, cex=0.5)
		#points(downstream_st, ii, pch=10, cex=0.5)
		text(mid, ii+0.1, tmparr[1], cex=2)
		mid <- (upstream_en+downstream_st)/2
		lines(c(upstream_en, mid), c(ii, ii-0.1), lwd=1.5)
		lines(c(mid, downstream_st), c(ii-0.1, ii), lwd=1.5)
		#lines(c(mid, target_en), c(ii-0.3, ii), lwd=1)
		#points(upstream_en, ii, pch=10, cex=0.5)
		#points(downstream_st, ii, pch=10, cex=0.5)
		text(mid, ii-0.15, tmparr[2], cex=2)
	}	
}



if(tmp[1] == "A5SS")
{

	Control <- read.table("Control_junction.txt", sep="\t", header=FALSE)
	Control <- as.vector(Control[,1])
	tmparr <- strsplit(Control, ",")[[1]]
	tmparr <- as.numeric(tmparr)
	tmparr[1] <- round(tmparr[1], 0)
	tmparr[2] <- round(tmparr[2], 0)
	if(strand == "+")
	{
		mid <- (target_en+downstream_st)/2
		lines(c(target_en, mid), c(1.2, 1.2+0.1), lwd=1.5)
		lines(c(mid, downstream_st), c(1.2+0.1, ii), lwd=1.5)
		text(mid, 1.2+0.15, tmparr[1], cex=2)
		mid <- (upstream_en+downstream_st)/2
		lines(c(upstream_en, mid), c(1.2, 1.2-0.1), lwd=1.5)
		lines(c(mid, downstream_st), c(1.2-0.1, 1.2), lwd=1.5)
		text(mid, 1.2-0.15, tmparr[2], cex=2)
	}
	if(strand == "-")
	{
		mid <- (upstream_en+target_st)/2
		lines(c(upstream_en, mid), c(1.2, 1.2+0.1), lwd=1.5)
		lines(c(mid, target_st), c(1.2+0.1, 1.2), lwd=1.5)
		text(mid, 1.2+0.15, tmparr[1], cex=2)
		mid <- (upstream_en+downstream_st)/2
		lines(c(upstream_en, mid), c(1.2, 1.2-0.1), lwd=1.5)
		lines(c(mid, downstream_st), c(1.2-0.1, 1.2), lwd=1.5)
		text(mid, 1.2-0.15, tmparr[2], cex=2)
	}




	Control <- read.table("Case_junction.txt", sep="\t", header=FALSE)
	Control <- as.vector(Control[,1])
	tmparr <- strsplit(Control, ",")[[1]]
	tmparr <- as.numeric(tmparr)
	tmparr[1] <- round(tmparr[1], 0)
	tmparr[2] <- round(tmparr[2], 0)
	if(strand == "+")
	{
		mid <- (target_en+downstream_st)/2
		lines(c(target_en, mid), c(0.2, 0.2+0.1), lwd=1.5)
		lines(c(mid, downstream_st), c(0.2+0.1, ii), lwd=1.5)
		text(mid, 0.2+0.15, tmparr[1], cex=2)
		mid <- (upstream_en+downstream_st)/2
		lines(c(upstream_en, mid), c(0.2, 0.2-0.1), lwd=1.5)
		lines(c(mid, downstream_st), c(0.2-0.1, 0.2), lwd=1.5)
		text(mid, 0.2-0.15, tmparr[2], cex=2)
	}
	if(strand == "-")
	{
		mid <- (upstream_en+target_st)/2
		#lines(c(target_st, target_en), c(3.3, 3.3), lty=3)
		lines(c(upstream_en, mid), c(0.2, 0.2+0.1), lwd=1.5)
		lines(c(mid, target_st), c(0.2+0.1, 0.2), lwd=1.5)
		#points(upstream_en, 1.5, pch=10, cex=0.5)
		#points(target_st, 1.5, pch=10, cex=0.5)
		text(mid, 0.2+0.15, tmparr[1], cex=2)
		mid <- (upstream_en+downstream_st)/2
		lines(c(upstream_en, mid), c(0.2, 0.2-0.1), lwd=1)
		lines(c(mid, downstream_st), c(0.2-0.1, 0.2), lwd=1)
		text(mid, 0.2-0.15, tmparr[2], cex=1.8)
	}

}	
	


if(tmp[1] == "SE")
{

	Control1 <- read.table("Control_junction.txt", sep="\t", header=FALSE)
	Control1 <- as.vector(Control1[,1])
	tmparr <- strsplit(Control1, ",")[[1]]
	tmparr <- as.numeric(tmparr)
	tmparr[1] <- round(tmparr[1], 0)
	tmparr[2] <- round(tmparr[2], 0)
	ii=1.2
	mid <- (upstream_en+target_st)/2
	lines(c(upstream_en, mid), c(ii, ii+0.1), lwd=1.5)
	lines(c(mid, target_st), c(ii+0.1, ii), lwd=1.5)
	mid <- (target_en+downstream_st)/2
	lines(c(target_en, mid), c(ii, ii+0.1), lwd=1.5)
	lines(c(mid, downstream_st), c(ii+0.1, ii), lwd=1.5)
	mid <- (target_st+target_en)/2
	text(mid, ii+0.15, paste("Total:", tmparr[1], sep=""), cex=2)
	mid <- (upstream_en+downstream_st)/2
	lines(c(upstream_en, mid), c(ii, ii-0.1), lwd=1.5)
	lines(c(mid, downstream_st), c(ii-0.1, ii), lwd=1.5)
	text(mid, ii-0.15, tmparr[2], cex=2)
	


	Control1 <- read.table("Case_junction.txt", sep="\t", header=FALSE)
	Control1 <- as.vector(Control1[,1])
	tmparr <- strsplit(Control1, ",")[[1]]
	tmparr <- as.numeric(tmparr)
	tmparr[1] <- round(tmparr[1], 0)
	tmparr[2] <- round(tmparr[2], 0)
	ii=0.2
	mid <- (upstream_en+target_st)/2
	lines(c(upstream_en, mid), c(ii, ii+0.1), lwd=1.5)
	lines(c(mid, target_st), c(ii+0.1, ii), lwd=1.5)
	mid <- (target_en+downstream_st)/2
	lines(c(target_en, mid), c(ii, ii+0.1), lwd=1.5)
	lines(c(mid, downstream_st), c(ii+0.1, ii), lwd=1.5)
	mid <- (target_st+target_en)/2
	text(mid, ii+0.15, paste("Total:", tmparr[1], sep=""), cex=2)	
	mid <- (upstream_en+downstream_st)/2
	lines(c(upstream_en, mid), c(ii, ii-0.1), lwd=1.5)
	lines(c(mid, downstream_st), c(ii-0.1, ii), lwd=1.5)
	points(upstream_en, ii, pch=10, cex=0.5)
	points(downstream_st, ii, pch=10, cex=0.5)
	text(mid, ii-0.15, tmparr[2], cex=2)	

}



if(tmp[1] == "RI")
{
	Control <- read.table("Control_junction.txt", sep="\t", header=FALSE)
	Control <- as.vector(Control[,1])
	tmparr <- strsplit(Control, ",")[[1]]
	tmparr <- as.numeric(tmparr)
	tmparr[1] <- round(tmparr[1], 0)
	tmparr[2] <- round(tmparr[2], 0)
	mid <- (target_en+target_st)/2
	text(mid, 1.2+0.15, tmparr[1], cex=2)
	lines(c(target_st, mid), c(1.2, 1.2-0.15), lwd=1)
	lines(c(mid, target_en), c(1.2-0.15, 1.2), lwd=1)
	text((target_en+target_st)/2, 1.2-0.15, tmparr[2], cex=2)


	Case <- read.table("Case_junction.txt", sep="\t", header=FALSE)
	Case <- as.vector(Case[,1])
	tmparr <- strsplit(Case, ",")[[1]]
	tmparr <- as.numeric(tmparr)
	tmparr[1] <- round(tmparr[1], 0)
	tmparr[2] <- round(tmparr[2], 0)
	mid <- (target_en+target_st)/2
	text(mid, 0.2+0.15, tmparr[1], cex=1.8)
	lines(c(target_st, mid), c(0.2, 0.2-0.15), lwd=1)
	lines(c(mid, target_en), c(0.2-0.15, 0.2), lwd=1)
	#points(target_st, ii, pch=10, cex=0.5)
	#points(target_en, ii, pch=10, cex=0.5)
	text((target_en+target_st)/2, 0.2-0.15, tmparr[2], cex=2)
	


}
	
samplelist <- sub("\\_.*", "", samplelist)
#mtext(rev(samplelist), side=2, at=c(0.5, 1.5), las=2, cex=1.5)
axis(2, at=c(0.7, 1.7),labels=c(samplelist[2], samplelist[1]), cex.axis=2.2, tick = FALSE, las=1)

##### whole gene structure
par(mar=c(4,13,0,3))
gtf1 <- read.table("include.gtf", sep="\t", header=FALSE)
tmp1 <-  strsplit(unique(as.vector(gtf1$V9)), "\\;\\s")
tmp1 <- tmp1[[1]][2]
tmp1 <- sub("transcript_id ", "", tmp1)
tmp1 <- sub("\\;", "", tmp1)

gtf2 <- read.table("exclude.gtf", sep="\t", header=FALSE)
tmp2 <-  strsplit(unique(as.vector(gtf2$V9)), "\\;\\s")
tmp2 <- tmp2[[1]][2]
tmp2 <- sub("transcript_id ", "", tmp2)
tmp2 <- sub("\\;", "", tmp2)

tids <- c(tmp1, tmp2)

pos <- rbind(gtf1[,c("V4", "V5")],gtf2[,c("V4", "V5")])
gst <- min(pos)
gen <- max(pos)
plot(NA,NA, col="red",xlim=c(gst, gen),ylim=c(0,3),frame.plot=FALSE, yaxt="n", ylab="",xaxt="n", xlab="")
distance <- seq(gst, gen, by=(gen-gst)/2)
#axis(1, at=distance,labels=round(distance,0),cex.axis=1.8, las=1)
axis(2, at=1:2,labels=rev(tids), cex.axis=1.5, tick = FALSE, las=1)

#dev.off()
gtf <- read.table("include.gtf", sep="\t", header=FALSE)
i=2
source(fun1)
i=1
gtf <- read.table("exclude.gtf", sep="\t", header=FALSE)
source(fun1)

tmp <- strsplit(ASE, "\\&")[[1]]
n <- length(tmp)
target_st <- tmp[8]
target_en <- tmp[9]
n <- length(tmp)
ASEgst <- min(as.numeric(tmp[5:n]))
ASEgen <- max(as.numeric(tmp[5:n]))
#arrows(ASEgst, 2.2, gst, 3, col = "black", length=0.1, lwd=1,lty=2)
#arrows(ASEgen, 2.2, gen, 3, col = "black", length=0.1, lwd=1,lty=2)
lines(c(ASEgst,gst), c(2.5, 3), col="black", lwd=2.5)
lines(c(ASEgen,gen), c(2.5, 3), col="black", lwd=2.5)






equal10 <- (gen-gst)/50	
arrowheight <- 0.01
if(strand == "+")
{

	arrows(gst+equal10*24, arrowheight, gst+equal10*26, arrowheight, col = "black", length=0.1, lwd=2)
	arrows(gst+equal10*9, arrowheight, gst+equal10*11, arrowheight, col = "black", length=0.1, lwd=2)
	arrows(gst+equal10*39, arrowheight, gst+equal10*41, arrowheight, col = "black", length=0.1, lwd=2)
}else
{
	arrows(gst+equal10*26, arrowheight, gst+equal10*24, arrowheight, col = "black", length=0.1, lwd=2)
	arrows(gst+equal10*11, arrowheight, gst+equal10*9, arrowheight, col = "black", length=0.1, lwd=2)
	arrows(gst+equal10*41, arrowheight, gst+equal10*39, arrowheight, col = "black", length=0.1, lwd=2)
}
dev.off()


