argsx <- commandArgs(TRUE)
GENE <- argsx[1]
pcd <- argsx[2]
samplelist <- strsplit(argsx[3], ",")[[1]]
fun1 <- paste(pcd, "/plot_gene_exon_structure_AS_FunG3.r", sep="")
genefile <- read.table(GENE, header=FALSE, sep="\t")
strand <- as.vector(genefile[,6])
gst <- genefile[,2]
gen <- genefile[,3]
otfile <- paste("Gene", ".pdf", sep="")
pdf(file=otfile, width=15, height=8)
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

	filename <- paste(samplelist[1], "_junction.txt", sep="")
	Control <- read.table(filename, sep="\t", header=FALSE)
	for(i in 1:dim(Control)[1])
	{
	tmp <- strsplit(as.vector(Control[i,1]), "\\&")[[1]]
    target_st <- tmp[8]
    target_en <- tmp[9]
    target_st <- as.numeric(target_st)
    target_en <- as.numeric(target_en)
	tmpin <- as.numeric(Control[i,2])
	tmpex <- as.numeric(Control[i,3])
	tmpin <- round(tmpin, 0)
	tmpex <- round(tmpex, 0)
	mid <- (target_en+target_st)/2
	text(mid, 1.5+0.2, tmpin, cex=1.8)
	lines(c(target_st, mid), c(1.5, 1.5-0.1), lwd=1)
	lines(c(mid, target_en), c(1.5-0.1, 1.5), lwd=1)
	text((target_en+target_st)/2, 1.5-0.1, tmpex, cex=1.8)
	}

	
	
	#Case
	filename <- paste(samplelist[2], ".txt", sep="")
	Case <- read.table(filename, sep="\t", header=FALSE)
	n <- dim(Case)[1]
	depth <- Case[,3]
	depth <- (depth-min(depth))/(max(depth)-min(depth)+0.1) 
	for(i in 1:n)
	{
		lines(c(Case[i,2], Case[i,2]), c(0.5,depth[i]+0.5), col="red")
	}

	filename <- paste(samplelist[2], "_junction.txt", sep="")
	Control <- read.table(filename, sep="\t", header=FALSE)
	for(i in 1:dim(Control)[1])
	{
	tmp <- strsplit(as.vector(Control[i,1]), "\\&")[[1]]
    target_st <- tmp[8]
    target_en <- tmp[9]
    target_st <- as.numeric(target_st)
    target_en <- as.numeric(target_en)
	tmpin <- as.numeric(Control[i,2])
	tmpex <- as.numeric(Control[i,3])
	tmpin <- round(tmpin, 0)
	tmpex <- round(tmpex, 0)
	#theta <- c(theta, tmparr[3])
	#ii=1-index+2
	mid <- (target_en+target_st)/2
	text(mid, 0.5+0.2, tmpin, cex=1.8)
	lines(c(target_st, mid), c(0.5, 0.5-0.1), lwd=1)
	lines(c(mid, target_en), c(0.5-0.1, 0.5), lwd=1)
	text((target_en+target_st)/2, 0.5-0.1, tmpex, cex=1.8)
	}
	
	
samplelist <- sub("\\_.*", "", samplelist)
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

plot(NA,NA, col="red",xlim=c(gst, gen),ylim=c(0,5),frame.plot=FALSE,yaxt="n", ylab="",xaxt="n", xlab="")
distance <- seq(gst, gen, by=(gen-gst)/2)

axis(1, at=distance,labels=round(distance,0),cex.axis=1.8, las=1)

#axis(2, at=c(1.5, 2, 4),labels=c("Domain",rev(tids)), cex.axis=1.4, tick = FALSE, las=1)

#dev.off()
gtf <- read.table("include.gtf", sep="\t", header=FALSE)
gtf <- gtf[order(gtf$V4), ]
i=4
source(fun1)
i=2
gtf <- read.table("exclude.gtf", sep="\t", header=FALSE)
gtf <- gtf[order(gtf$V4), ]
source(fun1)

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
		    if(exon_type_uniq[i] == "PF06925")
			{
				rect(exsts[j], 1-0.3, exens[j], 1+0.3, col="purple",border="purple")
			}
		    if(exon_type_uniq[i] == "PF04101")
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
legend1 <- c(legend1, c('5UTR','CDS','3UTR'))
legend2 <- c(legend2, c(15,15,15))
legend3 <- c(legend3, c("green","blue","red"))
if(is.element('PF06925', exon_type_uniq))
{
legend1 <- c(legend1, 'MGDG synthase')
legend2 <- c(legend2, 15)
legend3 <- c(legend3, "purple")
}
if(is.element('PF04101', exon_type_uniq))
{
legend1 <- c(legend1, 'Glycosyltransferase')
legend2 <- c(legend2, 15)
legend3 <- c(legend3, "yellow")
}
legend("bottomright",inset=0, yjust=-2, pt.cex=1.2, x.intersp=0.5, legend=legend1, pch=legend2 ,box.col="NA", col=legend3, horiz=TRUE, cex=1)

dev.off()

