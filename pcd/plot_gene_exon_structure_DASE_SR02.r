argsx <- commandArgs(TRUE)
ASE <- argsx[1]
pcd <- argsx[2]
samplelist <- strsplit(argsx[3], ",")[[1]]
fun1 <- paste(pcd, "/plot_gene_exon_structure_AS_FunG3.r", sep="")
tmp <- strsplit(ASE, "\\&")[[1]]
target_st <- NULL
target_en <- NULL
upstream_st <- NULL
upstream_en <- NULL
downstream_st <- NULL
downstream_en <- NULL
strand <- tmp[4]
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
target_st <- as.numeric(target_st)
target_en <- as.numeric(target_en)
upstream_st <- as.numeric(upstream_st)
upstream_en <- as.numeric(upstream_en)
downstream_st <- as.numeric(downstream_st)
downstream_en <- as.numeric(downstream_en)


n <- length(tmp)
gst <- min(as.numeric(tmp[5:n]))
gen <- max(as.numeric(tmp[5:n]))
otfile <- paste(ASE, ".jpeg", sep="")
otfile <- paste(ASE, ".tiff", sep="")
otfile <- paste(ASE, ".pdf", sep="")
pdf(file=otfile)
colorss <- c("lightblue", "red")
layout(matrix(c(1,2), 2, 1, byrow=TRUE), height=c(3.5,3.5))
par(mar=c(0,12,0,3))
plot(NA, NA,xlim=c(gst, gen), ylim=c(0, 2.5),type="n",col="grey", main="", cex.main=0.9, ylab="", xlab="Pos(bp)", frame.plot=FALSE, yaxt="n",xaxt="n")
distance <- seq(gst, gen, by=(gen-gst)/2)
equal10 <- (gen-gst)/50	
if(strand == "+")
{
	arrows(gst+equal10*24, 0, gst+equal10*26, 0, col = "black", length=0.1, lwd=1)
	arrows(gst+equal10*9, 0, gst+equal10*11, 0, col = "black", length=0.1, lwd=1)
	arrows(gst+equal10*39, 0, gst+equal10*41, 0, col = "black", length=0.1, lwd=1)
}else
{
	arrows(gst+equal10*26, 0, gst+equal10*24, 0, col = "black", length=0.1, lwd=1)
	arrows(gst+equal10*11, 0, gst+equal10*9, 0, col = "black", length=0.1, lwd=1)
	arrows(gst+equal10*41, 0, gst+equal10*39, 0, col = "black", length=0.1, lwd=1)
}
theta <- NULL
sss_n <- length(samplelist)



	filename <- paste(samplelist[1], ".txt", sep="")
	#Control1
	Control <- read.table(filename, sep="\t", header=FALSE)
	n <- dim(Control)[1]
	depth <- Control[,3]
	depth <- (depth-min(depth))/(max(depth)-min(depth)+0.1) 
	for(i in 1:n)
	{
		lines(c(Control[i,2], Control[i,2]), c(1.5,depth[i]+1.5), col="lightblue")
	}
	filename <- paste(samplelist[2], ".txt", sep="")
	#Control1
	Case <- read.table(filename, sep="\t", header=FALSE)
	n <- dim(Case)[1]
	depth <- Case[,3]
	depth <- (depth-min(depth))/(max(depth)-min(depth)+0.1) 
	for(i in 1:n)
	{
		lines(c(Case[i,2], Case[i,2]), c(0.5,depth[i]+0.5), col="red")
	}




	filename <- paste(samplelist[1], "_junction.txt", sep="")
	Control <- read.table(filename, sep="\t", header=FALSE)
	Control <- as.vector(Control[,1])
	tmparr <- strsplit(Control, ",")[[1]]
	tmparr <- as.numeric(tmparr)
	tmparr[1] <- round(tmparr[1], 0)
	tmparr[2] <- round(tmparr[2], 0)
	#theta <- c(theta, tmparr[3])
	if(strand == "+")
	{
	    ii <- 1.5
		mid <- (upstream_en+target_st)/2
		#lines(c(target_st, target_en), c(3.3, 3.3), lty=3)
		lines(c(upstream_en, mid), c(ii, ii+0.1), lwd=1)
		lines(c(mid, target_st), c(ii+0.1, ii), lwd=1)
		#points(upstream_en, ii, pch=10, cex=0.5)
		#points(target_st, ii, pch=10, cex=0.5)
		text(mid, ii+0.1, tmparr[1], cex=1.2)
		mid <- (upstream_en+downstream_st)/2
		lines(c(upstream_en, mid), c(ii, ii-0.1), lwd=1)
		lines(c(mid, downstream_st), c(ii-0.1, ii), lwd=1)
		#lines(c(mid, target_en), c(ii-0.3, ii), lwd=1)
		#points(upstream_en, ii, pch=10, cex=0.5)
		#points(downstream_st, ii, pch=10, cex=0.5)
		text(mid, ii-0.1, tmparr[2], cex=1.2)
	}
	if(strand == "-")
	{
	    ii <- 1.5
		mid <- (target_en+downstream_st)/2
		#lines(c(target_st, target_en), c(3.3, 3.3), lty=3)
		lines(c(target_en, mid), c(ii, ii+0.1), lwd=1)
		lines(c(mid, downstream_st), c(ii+0.1, ii), lwd=1)
		#points(target_en, ii, pch=10, cex=0.5)
		#points(downstream_st, ii, pch=10, cex=0.5)
		text(mid, ii+0.1, tmparr[1], cex=1.2)
		mid <- (upstream_en+downstream_st)/2
		lines(c(upstream_en, mid), c(ii, ii-0.1), lwd=1)
		lines(c(mid, downstream_st), c(ii-0.1, ii), lwd=1)
		#lines(c(mid, target_en), c(ii-0.3, ii), lwd=1)
		#points(upstream_en, ii, pch=10, cex=0.5)
		#points(downstream_st, ii, pch=10, cex=0.5)
		text(mid, ii-0.1, tmparr[2], cex=1.2)
	}	
	
	
	
	

	
	
	filename <- paste(samplelist[2], "_junction.txt", sep="")
	Case <- read.table(filename, sep="\t", header=FALSE)
	Case <- as.vector(Case[,1])
	tmparr <- strsplit(Case, ",")[[1]]
	tmparr <- as.numeric(tmparr)
	tmparr[1] <- round(tmparr[1], 0)
	tmparr[2] <- round(tmparr[2], 0)
	theta <- c(theta, tmparr[3])
	if(strand == "+")
	{
	    ii <- 0.5
		mid <- (upstream_en+target_st)/2
		#lines(c(target_st, target_en), c(3.3, 3.3), lty=3)
		lines(c(upstream_en, mid), c(ii, ii+0.1), lwd=1)
		lines(c(mid, target_st), c(ii+0.1, ii), lwd=1)
		#points(upstream_en, ii, pch=10, cex=0.5)
		#points(target_st, ii, pch=10, cex=0.5)
		text(mid, ii+0.1, tmparr[1], cex=1.2)
		mid <- (upstream_en+downstream_st)/2
		lines(c(upstream_en, mid), c(ii, ii-0.1), lwd=1)
		lines(c(mid, downstream_st), c(ii-0.1, ii), lwd=1)
		#lines(c(mid, target_en), c(ii-0.3, ii), lwd=1)
		#points(upstream_en, ii, pch=10, cex=0.5)
		#points(downstream_st, ii, pch=10, cex=0.5)
		text(mid, ii-0.1, tmparr[2], cex=1.2)
	}
	if(strand == "-")
	{
	    ii <- 0.5
		mid <- (target_en+downstream_st)/2
		#lines(c(target_st, target_en), c(3.3, 3.3), lty=3)
		lines(c(target_en, mid), c(ii, ii+0.1), lwd=1)
		lines(c(mid, downstream_st), c(ii+0.1, ii), lwd=1)
		#points(target_en, ii, pch=10, cex=0.5)
		#points(downstream_st, ii, pch=10, cex=0.5)
		text(mid, ii+0.1, tmparr[1], cex=1.2)
		mid <- (upstream_en+downstream_st)/2
		lines(c(upstream_en, mid), c(ii, ii-0.1), lwd=1)
		lines(c(mid, downstream_st), c(ii-0.1, ii), lwd=1)
		#lines(c(mid, target_en), c(ii-0.3, ii), lwd=1)
		#points(upstream_en, ii, pch=10, cex=0.5)
		#points(downstream_st, ii, pch=10, cex=0.5)
		text(mid, ii-0.1, tmparr[2], cex=1.2)
	}	

	
samplelist <- sub("\\_.*", "", samplelist)
#mtext(rev(samplelist), side=2, at=c(0.5, 1.5), las=2, cex=1.5)
axis(2, at=c(0.7, 1.7),labels=c("-Pi", "+Pi"), cex.axis=2, tick = FALSE, las=1)


##### whole gene structure
par(mar=c(4,12,0,3))
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
plot(NA,NA, col="red",xlim=c(gst, gen),ylim=c(0,6),frame.plot=FALSE,yaxt="n", ylab="",xaxt="n", xlab="")
distance <- seq(gst, gen, by=(gen-gst)/2)
axis(1, at=distance,labels=round(distance,0),cex.axis=1.8, las=1)
#axis(2, at=2:4,labels=c("Domain",rev(tids)), cex.axis=1.4, tick = FALSE, las=1)
#dev.off()
gtf <- read.table("include.gtf", sep="\t", header=FALSE)
gtf <- gtf[order(gtf$V4), ]
i=4
source(fun1)
i=2
gtf <- read.table("exclude.gtf", sep="\t", header=FALSE)
gtf <- gtf[order(gtf$V4), ]
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
lines(c(ASEgst,gst), c(5.5, 6), col="black", lwd=1)
lines(c(ASEgen,gen), c(5.5, 6), col="black", lwd=1)

#### domain
gtf <- read.table("Domain.gtf", sep="\t", header=FALSE)
gidtid <- as.vector(gtf[,"V9"])
allgids <- NULL
tmp <- strsplit(gidtid, "\\;\\s")
for(i in 1:length(tmp))
{
	allgids <- c(allgids, tmp[[i]][1])
}
allgids <- sub("gene_id ", "", allgids)
exon_type <- as.vector(gtf[,"V3"])
exon_type_uniq <- unique(exon_type)
for(i in 1:length(exon_type_uniq))
{

	exsts <- gtf[exon_type == exon_type_uniq[i], "V4"]
	exens <- gtf[exon_type == exon_type_uniq[i], "V5"]
	n_ex  <- length(exsts)
	if(n_ex > 0)
	{
		for(j in 1:n_ex)
		{
		    if(exon_type_uniq[i] == "RRM")
			{
				rect(exsts[j], 1-0.3, exens[j], 1+0.3, col="purple",border="purple")
			}
		    if(exon_type_uniq[i] == "Zinc_knuckle")
			{
				rect(exsts[j], 1-0.3, exens[j], 1+0.3, col="yellow",border="yellow")
			}
			
		}
	}
}
legend1 <- NULL
legend2 <- NULL
legend3 <- NULL
#if(is.element("RRM", exon_type_uniq) && is.element("Zinc_knuckle", exon_type_uniq))
#{
#legend1 <- c(legend1, c('UTR','CDS','Exon','RRM', 'Zinc_knuckle'))
if(is.element('RRM', exon_type_uniq))
{
#legend1 <- c('RRM', 'Zinc_knuckle')
legend1 <- c(legend1, 'RRM')
legend2 <- c(legend2, 15)
legend3 <- c(legend3, "purple")
}
if(is.element('Zinc_knuckle', exon_type_uniq))
{
#legend1 <- c('RRM', 'Zinc_knuckle')
legend1 <- c(legend1, 'Zinc_knuckle')
legend2 <- c(legend2, 15)
legend3 <- c(legend3, "yellow")
}
#legend2 <- c(legend2, c(15, 15, 15,15,15))
#legend2 <- c(15,15)
#legend3 <- c(legend3, c("darkred","blue", "lightblue","purple", "yellow"))
#legend("topright",cex=1.2,legend=legend1,pch=legend2,box.col="NA", col=legend3)
#legend("bottom",inset=0, yjust=-2, pt.cex=1.2, x.intersp=0.5, legend=legend1, pch=legend2 ,box.col="NA", col=legend3, horiz=TRUE, cex=1)
dev.off()


