library(DESeq2)
fpkm <- read.table("FeZnCuMn_gene_sample_TPM_filter.xls", sep="\t", header=T, row.names=1)
d <- read.table("FeZnCuMn_gene_sample_count_filter.xls", sep="\t", header=T, row.names=1)
cons1 <- c(rep("ControlR",3), rep("FeR", 3), rep("ZnR", 3), rep("CuR", 3), rep("MnR", 3))
design <- data.frame(row.names = colnames(d), condition=c(cons1))
conds <- design$condition
rsums <- apply(d, 1, sum)
d <- d[rsums > 0,]
dds <- DESeqDataSetFromMatrix(countData = d, colData = design, design=~condition)
dds <- DESeq(dds, fitType="local")
normalizedcount <- counts(dds,normalized=TRUE)
rawcount <- counts(dds,normalized=FALSE)
write.table(data.frame(geneid=rownames(normalizedcount),normalizedcount), file="FeZnCuMn_gene_sample_normalized_count.xls", sep="\t", quote=FALSE, row.names=FALSE)
normalizedcount <- log2(normalizedcount+1)
write.table(data.frame(geneid=rownames(normalizedcount),normalizedcount), file="FeZnCuMn_gene_sample_normalized_count_log2.xls", sep="\t", quote=FALSE,row.names=FALSE)



##### Fe_vs_Control
d_one <- d[,cons1 %in% c("ControlR", "FeR")]
cons_one <- c(rep("ControlR",3), rep("FeR", 3))
design_one <- data.frame(row.names = colnames(d_one), condition=c(cons_one))
rsum_one1 <- apply(d_one[,1:3], 1, sum)
rsum_one2 <- apply(d_one[,4:6], 1, sum)
d_one <- d_one[rsum_one1 >= 30 | rsum_one2 >= 30,]
dds_one <- DESeqDataSetFromMatrix(countData = d_one, colData = design_one, design=~condition)
dds_one <- DESeq(dds_one, fitType="local")
res <- results(dds_one, contrast=c("condition","FeR","ControlR"), cooksCutoff=FALSE)
res2 <- as.data.frame(res)
gids <- rownames(res2)
length(gids)
subsetraw <- rawcount[gids,design$condition %in% c("ControlR", "FeR")]
colnames(subsetraw) <- paste(colnames(subsetraw), "_raw_count", sep="")
subsetnorm <- normalizedcount[gids,design$condition %in% c("ControlR", "FeR")]
colnames(subsetnorm) <- paste(colnames(subsetnorm), "_normalized_count", sep="")
subsetFPKM <- fpkm[rownames(subsetnorm),design$condition %in% c("ControlR", "FeR")]
colnames(subsetFPKM) <- paste(colnames(subsetFPKM), "_TPM", sep="")
#res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("Fe_vs_Control_normalized_count_log2FoldChange", "Pvalue")])
res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("log2FoldChange", "pvalue")])
n <- length(colnames(res3))
colnames(res3)[(n-1):n] <- c("Fe_vs_Control_Log2FoldChange", "Pvalue")
fdrs <- p.adjust(res3$Pvalue , method ="fdr")
res3 <- data.frame(res3, FDR=fdrs)
res3 <- res3[order(res3$Pvalue, decreasing=FALSE), ]
diff_res3 <- res3[res3$FDR <= 0.05 & abs(res3$Fe_vs_Control_Log2FoldChange) >= 1,]
write.table(data.frame(Gene_locus_id=rownames(diff_res3), diff_res3), "Fe_vs_Control_gene_exp_significant.xls", sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(Gene_locus_id=rownames(res3), res3), "Fe_vs_Control_gene_exp.xls", sep="\t", row.names=FALSE, quote=FALSE)


##### Zn_vs_Control
d_one <- d[,cons1 %in% c("ControlR", "ZnR")]
cons_one <- c(rep("ControlR",3), rep("ZnR", 3))
design_one <- data.frame(row.names = colnames(d_one), condition=c(cons_one))
rsum_one1 <- apply(d_one[,1:3], 1, sum)
rsum_one2 <- apply(d_one[,4:6], 1, sum)
d_one <- d_one[rsum_one1 >= 30 | rsum_one2 >= 30,]
dds_one <- DESeqDataSetFromMatrix(countData = d_one, colData = design_one, design=~condition)
dds_one <- DESeq(dds_one, fitType="local")
res <- results(dds_one, contrast=c("condition","ZnR","ControlR"), cooksCutoff=FALSE)
res2 <- as.data.frame(res)
gids <- rownames(res2)
length(gids)
#alldata <- read.table(normcoutfile, sep="\t", header=TRUE, row.names=1)
#alldata <- alldata[gids,c(g1s,g2s)]
subsetraw <- rawcount[gids,design$condition %in% c("ControlR", "ZnR")]
colnames(subsetraw) <- paste(colnames(subsetraw), "_raw_count", sep="")
subsetnorm <- normalizedcount[gids,design$condition %in% c("ControlR", "ZnR")]
colnames(subsetnorm) <- paste(colnames(subsetnorm), "_normalized_count", sep="")
subsetFPKM <- fpkm[rownames(subsetnorm),design$condition %in% c("ControlR", "ZnR")]
colnames(subsetFPKM) <- paste(colnames(subsetFPKM), "_TPM", sep="")
#res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("Fe_vs_Control_normalized_count_log2FoldChange", "Pvalue")])
res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("log2FoldChange", "pvalue")])
n <- length(colnames(res3))
colnames(res3)[(n-1):n] <- c("Zn_vs_Control_Log2FoldChange", "Pvalue")
fdrs <- p.adjust(res3$Pvalue , method ="fdr")
res3 <- data.frame(res3, FDR=fdrs)
res3 <- res3[order(res3$Pvalue, decreasing=FALSE), ]
diff_res3 <- res3[res3$FDR <= 0.05 & abs(res3$Zn_vs_Control_normalized_count_log2FoldChange) >= 1,]
write.table(data.frame(Gene_locus_id=rownames(diff_res3), diff_res3), "Zn_vs_Control_gene_exp_significant.xls", sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(Gene_locus_id=rownames(res3), res3), "Zn_vs_Control_gene_exp.xls", sep="\t", row.names=FALSE, quote=FALSE)



##### Cu_vs_Control
d_one <- d[,cons1 %in% c("ControlR", "CuR")]
cons_one <- c(rep("ControlR",3), rep("CuR", 3))
design_one <- data.frame(row.names = colnames(d_one), condition=c(cons_one))
rsum_one1 <- apply(d_one[,1:3], 1, sum)
rsum_one2 <- apply(d_one[,4:6], 1, sum)
d_one <- d_one[rsum_one1 >= 30 | rsum_one2 >= 30,]
dds_one <- DESeqDataSetFromMatrix(countData = d_one, colData = design_one, design=~condition)
dds_one <- DESeq(dds_one, fitType="local")
res <- results(dds_one, contrast=c("condition","CuR","ControlR"), cooksCutoff=FALSE)
res2 <- as.data.frame(res)
gids <- rownames(res2)
length(gids)
#alldata <- read.table(normcoutfile, sep="\t", header=TRUE, row.names=1)
#alldata <- alldata[gids,c(g1s,g2s)]
subsetraw <- rawcount[gids,design$condition %in% c("ControlR", "CuR")]
colnames(subsetraw) <- paste(colnames(subsetraw), "_raw_count", sep="")
subsetnorm <- normalizedcount[gids,design$condition %in% c("ControlR", "CuR")]
colnames(subsetnorm) <- paste(colnames(subsetnorm), "_normalized_count", sep="")
subsetFPKM <- fpkm[rownames(subsetnorm),design$condition %in% c("ControlR", "CuR")]
colnames(subsetFPKM) <- paste(colnames(subsetFPKM), "_TPM", sep="")
#res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("Fe_vs_Control_normalized_count_log2FoldChange", "Pvalue")])
res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("log2FoldChange", "pvalue")])
n <- length(colnames(res3))
colnames(res3)[(n-1):n] <- c("Cu_vs_Control_Log2FoldChange", "Pvalue")
fdrs <- p.adjust(res3$Pvalue , method ="fdr")
res3 <- data.frame(res3, FDR=fdrs)
res3 <- res3[order(res3$Pvalue, decreasing=FALSE), ]
diff_res3 <- res3[res3$FDR <= 0.05 & abs(res3$Cu_vs_Control_normalized_count_log2FoldChange) >= 1,]
write.table(data.frame(Gene_locus_id=rownames(diff_res3), diff_res3), "Cu_vs_Control_gene_exp_significant.xls", sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(Gene_locus_id=rownames(res3), res3), "Cu_vs_Control_gene_exp.xls", sep="\t", row.names=FALSE, quote=FALSE)





##### Mn_vs_Control
d_one <- d[,cons1 %in% c("ControlR", "MnR")]
cons_one <- c(rep("ControlR",3), rep("MnR", 3))
design_one <- data.frame(row.names = colnames(d_one), condition=c(cons_one))
rsum_one1 <- apply(d_one[,1:3], 1, sum)
rsum_one2 <- apply(d_one[,4:6], 1, sum)
d_one <- d_one[rsum_one1 >= 30 | rsum_one2 >= 30,]
dds_one <- DESeqDataSetFromMatrix(countData = d_one, colData = design_one, design=~condition)
dds_one <- DESeq(dds_one, fitType="local")
res <- results(dds_one, contrast=c("condition","MnR","ControlR"), cooksCutoff=FALSE)
res2 <- as.data.frame(res)
gids <- rownames(res2)
length(gids)
#alldata <- read.table(normcoutfile, sep="\t", header=TRUE, row.names=1)
#alldata <- alldata[gids,c(g1s,g2s)]
subsetraw <- rawcount[gids,design$condition %in% c("ControlR", "MnR")]
colnames(subsetraw) <- paste(colnames(subsetraw), "_raw_count", sep="")
subsetnorm <- normalizedcount[gids,design$condition %in% c("ControlR", "MnR")]
colnames(subsetnorm) <- paste(colnames(subsetnorm), "_normalized_count", sep="")
subsetFPKM <- fpkm[rownames(subsetnorm),design$condition %in% c("ControlR", "MnR")]
colnames(subsetFPKM) <- paste(colnames(subsetFPKM), "_TPM", sep="")
#res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("Fe_vs_Control_normalized_count_log2FoldChange", "Pvalue")])
res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("log2FoldChange", "pvalue")])
n <- length(colnames(res3))
colnames(res3)[(n-1):n] <- c("Mn_vs_Control_Log2FoldChange", "Pvalue")
fdrs <- p.adjust(res3$Pvalue , method ="fdr")
res3 <- data.frame(res3, FDR=fdrs)
res3 <- res3[order(res3$Pvalue, decreasing=FALSE), ]
diff_res3 <- res3[res3$FDR <= 0.05 & abs(res3$Mn_vs_Control_normalized_count_log2FoldChange) >= 1,]
write.table(data.frame(Gene_locus_id=rownames(diff_res3), diff_res3), "Mn_vs_Control_gene_exp_significant.xls", sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(Gene_locus_id=rownames(res3), res3), "Mn_vs_Control_gene_exp.xls", sep="\t", row.names=FALSE, quote=FALSE)

