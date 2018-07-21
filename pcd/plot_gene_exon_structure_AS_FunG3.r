minEXON <- NULL
maxEXON <- NULL
exon_type <- as.vector(gtf$V3)
exsts <- gtf[exon_type == "exon", "V4"]
exens <- gtf[exon_type == "exon", "V5"]
if(length(exsts) > 0)
{
	strand <- gtf[, "V7"]
	strand <- as.vector(strand)
	strand <- unique(strand)
	n_ex  <- length(exsts)
	for(j in 1:n_ex)
	{
		rect(exsts[j], i-0.2, exens[j], i+0.2, col="blue")
		if(j < length(exsts))
		{
			lines(c(exens[j],exsts[j+1]),c(i, i), col = "black",lwd=1.8)
			
			if(strand == "+")
			{
				mid <- (exens[j]+exsts[j+1])/2
			}
			if(strand == "-")
			{
				mid <- (exens[j]+exsts[j+1])/2
				#lines(c(t_exens[j],t_exsts[j+1]),c(i, i), col = "grey",lwd=1.8)
				#arrows(exsts[j+1], i, mid, i, col = "black", length=0.1, lwd=1)
			}
		}
	}
	minEXON <- min(exsts)
	maxEXON <- max(exens)
}




minCDS <- NULL
maxCDS <- NULL
exon_type <- as.vector(gtf$V3)
exsts <- gtf[exon_type == "CDS", "V4"]
exens <- gtf[exon_type == "CDS", "V5"]
if(length(exsts) > 0)
{
	strand <- gtf[, "V7"]
	strand <- as.vector(strand)
	strand <- unique(strand)
	n_ex  <- length(exsts)
	for(j in 1:n_ex)
	{
		rect(exsts[j], i-0.3, exens[j], i+0.3, col="blue")
		if(j < length(exsts))
		{
			lines(c(exens[j],exsts[j+1]),c(i, i), col = "black",lwd=1.8)
			
			if(strand == "+")
			{
				mid <- (exens[j]+exsts[j+1])/2
			}
			if(strand == "-")
			{
				mid <- (exens[j]+exsts[j+1])/2
				#lines(c(t_exens[j],t_exsts[j+1]),c(i, i), col = "grey",lwd=1.8)
				#arrows(exsts[j+1], i, mid, i, col = "black", length=0.1, lwd=1)
			}
		}
	}
	minCDS <- min(exsts)
	maxCDS <- max(exens)
	if(strand == "+")
	{
	     text(minCDS, i+0.5, "ATG", cex=0.8)
	     text(maxCDS, i+0.5, "TGA", cex=0.8)
		 
	}
	if(strand == "-")
	{
	     text(maxCDS, i+0.5, "GTA", cex=0.8)
	     text(minCDS, i+0.5, "AGT", cex=0.8)
		 
	}
}



min5UTR <- NULL
max5UTR <- NULL
exsts <- gtf[exon_type == "five_prime_UTR", "V4"]
exens <- gtf[exon_type == "five_prime_UTR", "V5"]
if(length(exsts) > 0)
{
	strand <- gtf[, "V7"]
	strand <- as.vector(strand)
	strand <- unique(strand)
	n_ex  <- length(exsts)
	for(j in 1:n_ex)
	{
		rect(exsts[j], i-0.2, exens[j], i+0.2, col="green")
		if(j < length(exsts))
		{
			lines(c(exens[j],exsts[j+1]),c(i, i), col = "black",lwd=1.8)
			if(strand == "+")
			{
				mid <- (exens[j]+exsts[j+1])/2
			}
			if(strand == "-")
			{
				mid <- (exens[j]+exsts[j+1])/2
				#lines(c(t_exens[j],t_exsts[j+1]),c(i, i), col = "grey",lwd=1.8)
				#arrows(exsts[j+1], i, mid, i, col = "black", length=0.1, lwd=1)
			}
		}
	}
	min5UTR <- min(exsts)
	max5UTR <- max(exens)
}


min3UTR <- NULL
max3UTR <- NULL
exsts <- gtf[exon_type == "three_prime_UTR", "V4"]
exens <- gtf[exon_type == "three_prime_UTR", "V5"]
if(length(exsts) > 0)
{
	strand <- gtf[, "V7"]
	strand <- as.vector(strand)
	strand <- unique(strand)
	n_ex  <- length(exsts)
	for(j in 1:n_ex)
	{
		rect(exsts[j], i-0.2, exens[j], i+0.2, col="red")
		if(j < length(exsts))
		{
			lines(c(exens[j],exsts[j+1]),c(i, i), col = "black",lwd=1.8)
			if(strand == "+")
			{
				mid <- (exens[j]+exsts[j+1])/2
				#lines(c(t_exens[j],t_exsts[j+1]),c(i, i), col = "grey",lwd=1.8)
				#arrows(exens[j], i, mid, i, col = "black", length=0.1, lwd=1)				
			}
			if(strand == "-")
			{
				mid <- (exens[j]+exsts[j+1])/2
				#lines(c(t_exens[j],t_exsts[j+1]),c(i, i), col = "grey",lwd=1.8)
				#arrows(exsts[j+1], i, mid, i, col = "black", length=0.1, lwd=1)
			}
		}
	}
	min3UTR <- min(exsts)
	max3UTR <- max(exens)
}


if(!is.null(min5UTR) && !is.null(maxCDS))
{
	if(min5UTR > maxCDS && !is.null(min5UTR))
	{
		lines(c(maxCDS,min5UTR),c(i, i), col = "black",lwd=1.8)  
	}
}

if(!is.null(min3UTR) && !is.null(maxCDS))
{


if(min3UTR > maxCDS)
{
    lines(c(maxCDS,min3UTR),c(i, i), col = "black",lwd=1.8)  
}

}


if(!is.null(max5UTR) && !is.null(minCDS))
{



if(max5UTR < minCDS  && !is.null(max5UTR))
{
    lines(c(max5UTR,minCDS),c(i, i), col = "black",lwd=1.8)  
}

}

if(!is.null(max3UTR) && !is.null(minCDS))
{
	if(max3UTR < minCDS)
	{
		lines(c(max3UTR,minCDS),c(i, i), col = "black",lwd=1.8)  
	}
}