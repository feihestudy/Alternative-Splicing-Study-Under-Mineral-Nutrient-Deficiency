library(gplots)
library(pheatmap)
FeZnCuMn <- read.table("SRs_FeZnCuMn_DASG_overlap.xls", sep="\t", header=TRUE)
rownames(FeZnCuMn) <- as.vector(FeZnCuMn[,1])
FeZnCuMn[,1] <- NULL
FeZnCuMn <- as.matrix(FeZnCuMn)
root <- read.table("SRs_PiRoot_DASG_overlap.xls", sep="\t", header=TRUE)
rownames(root) <- as.vector(root[,1])
root[,1] <- NULL
root <- as.matrix(root)
colnames(root) <- NULL
root2 <- NULL
for(i in 1:dim(root)[1])
{
   one <- as.numeric(root[i,])
   if(sum(one==1) >= 1)
   {
      root2 <- c(root2, 1)
	  next
   }
   if(sum(one==-1) >= 1)
   {
      root2 <- c(root2, -1)
	  next
   }
   root2 <- c(root2, 0)
}
shoot <- read.table("SRs_PiShoot_DASG_overlap.xls", sep="\t", header=TRUE)
rownames(shoot) <- as.vector(shoot[,1])
shoot[,1] <- NULL
shoot <- as.matrix(shoot)
colnames(shoot) <- NULL
shoot2 <- NULL
for(i in 1:dim(shoot)[1])
{
   one <- as.numeric(shoot[i,])
   if(sum(one==1) >= 1)
   {
      shoot2 <- c(shoot2, 1)
	  next
   }
   if(sum(one==-1) >= 1)
   {
      shoot2 <- c(shoot2, -1)
	  next
   }
   shoot2 <- c(shoot2, 0)
}
rootshoot <- cbind(FeZnCuMn, root2, shoot2)
colnames(rootshoot) <- c("Fe","Zn", "Cu", "Mn","PiR","PiS")
annotation <- data.frame(Tissue = c("Root", "Root", "Root", "Root","Root","Shoot"), Condition=c("Fe","Zn", "Cu", "Mn","Pi", "Pi"))
rownames(annotation) <- colnames(rootshoot) # check out the row names of annotation
Tissue <- c("grey", "green")
names(Tissue) <- c("Root", "Shoot")
Condition <- c("red", "grey26","green4", "hotpink", "gold")
names(Condition) <- c("Fe", "Zn", "Cu", "Mn","Pi")
anno_colors <- list(Tissue = Tissue,Condition = Condition)
cols <- colorRampPalette(c("white", "red"))(12)
pheatmap(rootshoot, fontsize = 15, filename="Figure5B.pdf", annotation_colors = anno_colors, col=cols, annotation = annotation, scale="none", show_rownames = TRUE, show_colnames = FALSE, cluster_rows= FALSE, cluster_cols = FALSE, fontsize_row=13, fontsize_col=12, border_color="black", hclustfun=function(c) (hclust(c, "ave")))


