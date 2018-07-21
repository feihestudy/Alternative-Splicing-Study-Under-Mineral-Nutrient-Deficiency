library(DESeq2)
fpkm <- read.table("PiRoot_gene_sample_TPM_filter.xls", sep="\t", header=T, row.names=1)
d <- read.table("PiRoot_gene_sample_count_filter.xls", sep="\t", header=T, row.names=1)
cons1 <- c(rep("H1piPlus_R",3), rep("H1piMinus_R", 3), rep("H6piPlus_R", 3), rep("H6piMinus_R", 3), rep("H24piPlus_R", 3), rep("H24piMinus_R", 3), rep("D3piPlus_R", 3), rep("D3piMinus_R", 3), rep("D7piPlus_R", 3), rep("D7piMinus_R", 3), rep("D21piPlus_R", 3), rep("D21piMinus_R", 3), rep("D21H1piPlus_R", 3), rep("D21H1piMinus_R", 3), rep("D21H1piPlusRec_R", 3), rep("D21H6piPlus_R", 3), rep("D21H6piMinus_R", 3), rep("D21H6piPlusRec_R", 3), rep("D21H24piPlus_R", 3), rep("D21H24piMinus_R", 3), rep("D21H24piPlusRec_R", 3))
#design <- data.frame(row.names = colnames(d), condition=c(cons1,cons2,cons3))
design <- data.frame(row.names = colnames(d), condition=c(cons1))
#design <- data.frame(row.names = colnames(d), condition=c("O", "O", "O", "YC", "YC", "YC"))
conds <- design$condition
rsums <- apply(d, 1, sum)
d <- d[rsums > 0,]
dds <- DESeqDataSetFromMatrix(countData = d, colData = design, design=~condition)
dds <- DESeq(dds, fitType="local")
normalizedcount <- counts(dds,normalized=TRUE)
rawcount <- counts(dds,normalized=FALSE)
write.table(data.frame(geneid=rownames(normalizedcount),normalizedcount), file="PiRoot_gene_sample_normalized_count.xls", sep="\t", quote=FALSE, row.names=FALSE)
normalizedcount <- log2(normalizedcount+1)
write.table(data.frame(geneid=rownames(normalizedcount),normalizedcount), file="PiRoot_gene_sample_normalized_count_log2.xls", sep="\t", quote=FALSE, row.names=FALSE)



##### H1piMinus_vs_H1piPlus
d_one <- d[,cons1 %in% c("H1piPlus_R", "H1piMinus_R")]
cons_one <- c(rep("H1piPlus_R",3), rep("H1piMinus_R", 3))
design_one <- data.frame(row.names = colnames(d_one), condition=c(cons_one))
rsum_one1 <- apply(d_one[,1:3], 1, sum)
rsum_one2 <- apply(d_one[,4:6], 1, sum)
d_one <- d_one[rsum_one1 >= 30 | rsum_one2 >= 30,]
dds_one <- DESeqDataSetFromMatrix(countData = d_one, colData = design_one, design=~condition)
dds_one <- DESeq(dds_one, fitType="local")
res <- results(dds_one, contrast=c("condition","H1piMinus_R","H1piPlus_R"), cooksCutoff=FALSE)
res2 <- as.data.frame(res)
gids <- rownames(res2)
length(gids)
#alldata <- read.table(normcoutfile, sep="\t", header=TRUE, row.names=1)
#alldata <- alldata[gids,c(g1s,g2s)]
subsetraw <- rawcount[gids,design$condition %in% c("H1piPlus_R", "H1piMinus_R")]
colnames(subsetraw) <- paste(colnames(subsetraw), "_raw_count", sep="")
subsetnorm <- normalizedcount[gids,design$condition %in% c("H1piPlus_R", "H1piMinus_R")]
colnames(subsetnorm) <- paste(colnames(subsetnorm), "_normalized_count", sep="")
subsetFPKM <- fpkm[rownames(subsetnorm),design$condition %in% c("H1piPlus_R", "H1piMinus_R")]
colnames(subsetFPKM) <- paste(colnames(subsetFPKM), "_TPM", sep="")
#res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("Fe_vs_Control_normalized_count_log2FoldChange", "Pvalue")])
res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("log2FoldChange", "pvalue")])
n <- length(colnames(res3))
colnames(res3)[(n-1):n] <- c("H1piMinus_vs_H1piPlus_Log2FoldChange", "Pvalue")
fdrs <- p.adjust(res3$Pvalue , method ="fdr")
res3 <- data.frame(res3, FDR=fdrs)
res3 <- res3[order(res3$Pvalue, decreasing=FALSE), ]
diff_res3 <- res3[res3$FDR <= 0.05 & abs(res3$H1piMinus_vs_H1piPlus_normalized_count_log2FoldChange) >= 1,]
write.table(data.frame(Gene_locus_id=rownames(diff_res3), diff_res3), "H1piMinus_vs_H1piPlus_R_gene_exp_significant.xls", sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(Gene_locus_id=rownames(res3), res3), "H1piMinus_vs_H1piPlus_R_gene_exp.xls", sep="\t", row.names=FALSE, quote=FALSE)



##### H6piMinus_vs_H6piPlus
d_one <- d[,cons1 %in% c("H6piPlus_R", "H6piMinus_R")]
cons_one <- c(rep("H6piPlus_R",3), rep("H6piMinus_R", 3))
design_one <- data.frame(row.names = colnames(d_one), condition=c(cons_one))
rsum_one1 <- apply(d_one[,1:3], 1, sum)
rsum_one2 <- apply(d_one[,4:6], 1, sum)
d_one <- d_one[rsum_one1 >= 30 | rsum_one2 >= 30,]
dds_one <- DESeqDataSetFromMatrix(countData = d_one, colData = design_one, design=~condition)
dds_one <- DESeq(dds_one, fitType="local")
res <- results(dds_one, contrast=c("condition","H6piMinus_R","H6piPlus_R"), cooksCutoff=FALSE)
res2 <- as.data.frame(res)
gids <- rownames(res2)
length(gids)
#alldata <- read.table(normcoutfile, sep="\t", header=TRUE, row.names=1)
#alldata <- alldata[gids,c(g1s,g2s)]
subsetraw <- rawcount[gids,design$condition %in% c("H6piPlus_R", "H6piMinus_R")]
colnames(subsetraw) <- paste(colnames(subsetraw), "_raw_count", sep="")
subsetnorm <- normalizedcount[gids,design$condition %in% c("H6piPlus_R", "H6piMinus_R")]
colnames(subsetnorm) <- paste(colnames(subsetnorm), "_normalized_count", sep="")
subsetFPKM <- fpkm[rownames(subsetnorm),design$condition %in% c("H6piPlus_R", "H6piMinus_R")]
colnames(subsetFPKM) <- paste(colnames(subsetFPKM), "_TPM", sep="")
#res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("Fe_vs_Control_normalized_count_log2FoldChange", "Pvalue")])
res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("log2FoldChange", "pvalue")])
n <- length(colnames(res3))
colnames(res3)[(n-1):n] <- c("H6piMinus_vs_H6piPlus_Log2FoldChange", "Pvalue")
fdrs <- p.adjust(res3$Pvalue , method ="fdr")
res3 <- data.frame(res3, FDR=fdrs)
res3 <- res3[order(res3$Pvalue, decreasing=FALSE), ]
diff_res3 <- res3[res3$FDR <= 0.05 & abs(res3$H6piMinus_vs_H6piPlus_normalized_count_log2FoldChange) >= 1,]
write.table(data.frame(Gene_locus_id=rownames(diff_res3), diff_res3), "H6piMinus_vs_H6piPlus_R_gene_exp_significant.xls", sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(Gene_locus_id=rownames(res3), res3), "H6piMinus_vs_H6piPlus_R_gene_exp.xls", sep="\t", row.names=FALSE, quote=FALSE)


##### H24piMinus_vs_H24piPlus
d_one <- d[,cons1 %in% c("H24piPlus_R", "H24piMinus_R")]
cons_one <- c(rep("H24piPlus_R",3), rep("H24piMinus_R", 3))
design_one <- data.frame(row.names = colnames(d_one), condition=c(cons_one))
rsum_one1 <- apply(d_one[,1:3], 1, sum)
rsum_one2 <- apply(d_one[,4:6], 1, sum)
d_one <- d_one[rsum_one1 >= 30 | rsum_one2 >= 30,]
dds_one <- DESeqDataSetFromMatrix(countData = d_one, colData = design_one, design=~condition)
dds_one <- DESeq(dds_one, fitType="local")
res <- results(dds_one, contrast=c("condition","H24piMinus_R","H24piPlus_R"), cooksCutoff=FALSE)
res2 <- as.data.frame(res)
gids <- rownames(res2)
length(gids)
#alldata <- read.table(normcoutfile, sep="\t", header=TRUE, row.names=1)
#alldata <- alldata[gids,c(g1s,g2s)]
subsetraw <- rawcount[gids,design$condition %in% c("H24piPlus_R", "H24piMinus_R")]
colnames(subsetraw) <- paste(colnames(subsetraw), "_raw_count", sep="")
subsetnorm <- normalizedcount[gids,design$condition %in% c("H24piPlus_R", "H24piMinus_R")]
colnames(subsetnorm) <- paste(colnames(subsetnorm), "_normalized_count", sep="")
subsetFPKM <- fpkm[rownames(subsetnorm),design$condition %in% c("H24piPlus_R", "H24piMinus_R")]
colnames(subsetFPKM) <- paste(colnames(subsetFPKM), "_TPM", sep="")
#res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("Fe_vs_Control_normalized_count_log2FoldChange", "Pvalue")])
res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("log2FoldChange", "pvalue")])
n <- length(colnames(res3))
colnames(res3)[(n-1):n] <- c("H24piMinus_vs_H24piPlus_Log2FoldChange", "Pvalue")
fdrs <- p.adjust(res3$Pvalue , method ="fdr")
res3 <- data.frame(res3, FDR=fdrs)
res3 <- res3[order(res3$Pvalue, decreasing=FALSE), ]
diff_res3 <- res3[res3$FDR <= 0.05 & abs(res3$H24piMinus_vs_H24piPlus_normalized_count_log2FoldChange) >= 1,]
write.table(data.frame(Gene_locus_id=rownames(diff_res3), diff_res3), "H24piMinus_vs_H24piPlus_R_gene_exp_significant.xls", sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(Gene_locus_id=rownames(res3), res3), "H24piMinus_vs_H24piPlus_R_gene_exp.xls", sep="\t", row.names=FALSE, quote=FALSE)



##### D3piMinus_vs_D3piPlus
d_one <- d[,cons1 %in% c("D3piPlus_R", "D3piMinus_R")]
cons_one <- c(rep("D3piPlus_R",3), rep("D3piMinus_R", 3))
design_one <- data.frame(row.names = colnames(d_one), condition=c(cons_one))
rsum_one1 <- apply(d_one[,1:3], 1, sum)
rsum_one2 <- apply(d_one[,4:6], 1, sum)
d_one <- d_one[rsum_one1 >= 30 | rsum_one2 >= 30,]
dds_one <- DESeqDataSetFromMatrix(countData = d_one, colData = design_one, design=~condition)
dds_one <- DESeq(dds_one, fitType="local")
res <- results(dds_one, contrast=c("condition","D3piMinus_R","D3piPlus_R"), cooksCutoff=FALSE)
res2 <- as.data.frame(res)
gids <- rownames(res2)
length(gids)
#alldata <- read.table(normcoutfile, sep="\t", header=TRUE, row.names=1)
#alldata <- alldata[gids,c(g1s,g2s)]
subsetraw <- rawcount[gids,design$condition %in% c("D3piPlus_R", "D3piMinus_R")]
colnames(subsetraw) <- paste(colnames(subsetraw), "_raw_count", sep="")
subsetnorm <- normalizedcount[gids,design$condition %in% c("D3piPlus_R", "D3piMinus_R")]
colnames(subsetnorm) <- paste(colnames(subsetnorm), "_normalized_count", sep="")
subsetFPKM <- fpkm[rownames(subsetnorm),design$condition %in% c("D3piPlus_R", "D3piMinus_R")]
colnames(subsetFPKM) <- paste(colnames(subsetFPKM), "_TPM", sep="")
#res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("Fe_vs_Control_normalized_count_log2FoldChange", "Pvalue")])
res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("log2FoldChange", "pvalue")])
n <- length(colnames(res3))
colnames(res3)[(n-1):n] <- c("D3piMinus_vs_D3piPlus_Log2FoldChange", "Pvalue")
fdrs <- p.adjust(res3$Pvalue , method ="fdr")
res3 <- data.frame(res3, FDR=fdrs)
res3 <- res3[order(res3$Pvalue, decreasing=FALSE), ]
diff_res3 <- res3[res3$FDR <= 0.05 & abs(res3$D3piMinus_vs_D3piPlus_normalized_count_log2FoldChange) >= 1,]
write.table(data.frame(Gene_locus_id=rownames(diff_res3), diff_res3), "D3piMinus_vs_D3piPlus_R_gene_exp_significant.xls", sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(Gene_locus_id=rownames(res3), res3), "D3piMinus_vs_D3piPlus_R_gene_exp.xls", sep="\t", row.names=FALSE, quote=FALSE)



##### D7piMinus_vs_D7piPlus
d_one <- d[,cons1 %in% c("D7piPlus_R", "D7piMinus_R")]
cons_one <- c(rep("D7piPlus_R",3), rep("D7piMinus_R", 3))
design_one <- data.frame(row.names = colnames(d_one), condition=c(cons_one))
rsum_one1 <- apply(d_one[,1:3], 1, sum)
rsum_one2 <- apply(d_one[,4:6], 1, sum)
d_one <- d_one[rsum_one1 >= 30 | rsum_one2 >= 30,]
dds_one <- DESeqDataSetFromMatrix(countData = d_one, colData = design_one, design=~condition)
dds_one <- DESeq(dds_one, fitType="local")
res <- results(dds_one, contrast=c("condition","D7piMinus_R","D7piPlus_R"), cooksCutoff=FALSE)
res2 <- as.data.frame(res)
gids <- rownames(res2)
length(gids)
#alldata <- read.table(normcoutfile, sep="\t", header=TRUE, row.names=1)
#alldata <- alldata[gids,c(g1s,g2s)]
subsetraw <- rawcount[gids,design$condition %in% c("D7piPlus_R", "D7piMinus_R")]
colnames(subsetraw) <- paste(colnames(subsetraw), "_raw_count", sep="")
subsetnorm <- normalizedcount[gids,design$condition %in% c("D7piPlus_R", "D7piMinus_R")]
colnames(subsetnorm) <- paste(colnames(subsetnorm), "_normalized_count", sep="")
subsetFPKM <- fpkm[rownames(subsetnorm),design$condition %in% c("D7piPlus_R", "D7piMinus_R")]
colnames(subsetFPKM) <- paste(colnames(subsetFPKM), "_TPM", sep="")
#res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("Fe_vs_Control_normalized_count_log2FoldChange", "Pvalue")])
res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("log2FoldChange", "pvalue")])
n <- length(colnames(res3))
colnames(res3)[(n-1):n] <- c("D7piMinus_vs_D7piPlus_Log2FoldChange", "Pvalue")
fdrs <- p.adjust(res3$Pvalue , method ="fdr")
res3 <- data.frame(res3, FDR=fdrs)
res3 <- res3[order(res3$Pvalue, decreasing=FALSE), ]
diff_res3 <- res3[res3$FDR <= 0.05 & abs(res3$D7piMinus_vs_D7piPlus_normalized_count_log2FoldChange) >= 1,]
write.table(data.frame(Gene_locus_id=rownames(diff_res3), diff_res3), "D7piMinus_vs_D7piPlus_R_gene_exp_significant.xls", sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(Gene_locus_id=rownames(res3), res3), "D7piMinus_vs_D7piPlus_R_gene_exp.xls", sep="\t", row.names=FALSE, quote=FALSE)


##### D21piMinus_vs_D21piPlus
d_one <- d[,cons1 %in% c("D21piPlus_R", "D21piMinus_R")]
cons_one <- c(rep("D21piPlus_R",3), rep("D21piMinus_R", 3))
design_one <- data.frame(row.names = colnames(d_one), condition=c(cons_one))
rsum_one1 <- apply(d_one[,1:3], 1, sum)
rsum_one2 <- apply(d_one[,4:6], 1, sum)
d_one <- d_one[rsum_one1 >= 30 | rsum_one2 >= 30,]
dds_one <- DESeqDataSetFromMatrix(countData = d_one, colData = design_one, design=~condition)
dds_one <- DESeq(dds_one, fitType="local")
res <- results(dds_one, contrast=c("condition","D21piMinus_R","D21piPlus_R"), cooksCutoff=FALSE)
res2 <- as.data.frame(res)
gids <- rownames(res2)
length(gids)
#alldata <- read.table(normcoutfile, sep="\t", header=TRUE, row.names=1)
#alldata <- alldata[gids,c(g1s,g2s)]
subsetraw <- rawcount[gids,design$condition %in% c("D21piPlus_R", "D21piMinus_R")]
colnames(subsetraw) <- paste(colnames(subsetraw), "_raw_count", sep="")
subsetnorm <- normalizedcount[gids,design$condition %in% c("D21piPlus_R", "D21piMinus_R")]
colnames(subsetnorm) <- paste(colnames(subsetnorm), "_normalized_count", sep="")
subsetFPKM <- fpkm[rownames(subsetnorm),design$condition %in% c("D21piPlus_R", "D21piMinus_R")]
colnames(subsetFPKM) <- paste(colnames(subsetFPKM), "_TPM", sep="")
#res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("Fe_vs_Control_normalized_count_log2FoldChange", "Pvalue")])
res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("log2FoldChange", "pvalue")])
n <- length(colnames(res3))
colnames(res3)[(n-1):n] <- c("D21piMinus_vs_D21piPlus_Log2FoldChange", "Pvalue")
fdrs <- p.adjust(res3$Pvalue , method ="fdr")
res3 <- data.frame(res3, FDR=fdrs)
res3 <- res3[order(res3$Pvalue, decreasing=FALSE), ]
diff_res3 <- res3[res3$FDR <= 0.05 & abs(res3$D21piMinus_vs_D21piPlus_normalized_count_log2FoldChange) >= 1,]
write.table(data.frame(Gene_locus_id=rownames(diff_res3), diff_res3), "D21piMinus_vs_D21piPlus_R_gene_exp_significant.xls", sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(Gene_locus_id=rownames(res3), res3), "D21piMinus_vs_D21piPlus_R_gene_exp.xls", sep="\t", row.names=FALSE, quote=FALSE)


##### D21H1piMinus_vs_D21H1piPlus
d_one <- d[,cons1 %in% c("D21H1piPlus_R", "D21H1piMinus_R")]
cons_one <- c(rep("D21H1piPlus_R",3), rep("D21H1piMinus_R", 3))
design_one <- data.frame(row.names = colnames(d_one), condition=c(cons_one))
rsum_one1 <- apply(d_one[,1:3], 1, sum)
rsum_one2 <- apply(d_one[,4:6], 1, sum)
d_one <- d_one[rsum_one1 >= 30 | rsum_one2 >= 30,]
dds_one <- DESeqDataSetFromMatrix(countData = d_one, colData = design_one, design=~condition)
dds_one <- DESeq(dds_one, fitType="local")
res <- results(dds_one, contrast=c("condition","D21H1piMinus_R","D21H1piPlus_R"), cooksCutoff=FALSE)
res2 <- as.data.frame(res)
gids <- rownames(res2)
length(gids)
#alldata <- read.table(normcoutfile, sep="\t", header=TRUE, row.names=1)
#alldata <- alldata[gids,c(g1s,g2s)]
subsetraw <- rawcount[gids,design$condition %in% c("D21H1piPlus_R", "D21H1piMinus_R")]
colnames(subsetraw) <- paste(colnames(subsetraw), "_raw_count", sep="")
subsetnorm <- normalizedcount[gids,design$condition %in% c("D21H1piPlus_R", "D21H1piMinus_R")]
colnames(subsetnorm) <- paste(colnames(subsetnorm), "_normalized_count", sep="")
subsetFPKM <- fpkm[rownames(subsetnorm),design$condition %in% c("D21H1piPlus_R", "D21H1piMinus_R")]
colnames(subsetFPKM) <- paste(colnames(subsetFPKM), "_TPM", sep="")
#res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("Fe_vs_Control_normalized_count_log2FoldChange", "Pvalue")])
res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("log2FoldChange", "pvalue")])
n <- length(colnames(res3))
colnames(res3)[(n-1):n] <- c("D21H1piMinus_vs_D21H1piPlus_Log2FoldChange", "Pvalue")
fdrs <- p.adjust(res3$Pvalue , method ="fdr")
res3 <- data.frame(res3, FDR=fdrs)
res3 <- res3[order(res3$Pvalue, decreasing=FALSE), ]
diff_res3 <- res3[res3$FDR <= 0.05 & abs(res3$D21H1piMinus_vs_D21H1piPlus_normalized_count_log2FoldChange) >= 1,]
write.table(data.frame(Gene_locus_id=rownames(diff_res3), diff_res3), "D21H1piMinus_vs_D21H1piPlus_R_gene_exp_significant.xls", sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(Gene_locus_id=rownames(res3), res3), "D21H1piMinus_vs_D21H1piPlus_R_gene_exp.xls", sep="\t", row.names=FALSE, quote=FALSE)


##### D21H6piMinus_vs_D21H6piPlus
d_one <- d[,cons1 %in% c("D21H6piPlus_R", "D21H6piMinus_R")]
cons_one <- c(rep("D21H6piPlus_R",3), rep("D21H6piMinus_R", 3))
design_one <- data.frame(row.names = colnames(d_one), condition=c(cons_one))
rsum_one1 <- apply(d_one[,1:3], 1, sum)
rsum_one2 <- apply(d_one[,4:6], 1, sum)
d_one <- d_one[rsum_one1 >= 30 | rsum_one2 >= 30,]
dds_one <- DESeqDataSetFromMatrix(countData = d_one, colData = design_one, design=~condition)
dds_one <- DESeq(dds_one, fitType="local")
res <- results(dds_one, contrast=c("condition","D21H6piMinus_R","D21H6piPlus_R"), cooksCutoff=FALSE)
res2 <- as.data.frame(res)
gids <- rownames(res2)
length(gids)
#alldata <- read.table(normcoutfile, sep="\t", header=TRUE, row.names=1)
#alldata <- alldata[gids,c(g1s,g2s)]
subsetraw <- rawcount[gids,design$condition %in% c("D21H6piPlus_R", "D21H6piMinus_R")]
colnames(subsetraw) <- paste(colnames(subsetraw), "_raw_count", sep="")
subsetnorm <- normalizedcount[gids,design$condition %in% c("D21H6piPlus_R", "D21H6piMinus_R")]
colnames(subsetnorm) <- paste(colnames(subsetnorm), "_normalized_count", sep="")
subsetFPKM <- fpkm[rownames(subsetnorm),design$condition %in% c("D21H6piPlus_R", "D21H6piMinus_R")]
colnames(subsetFPKM) <- paste(colnames(subsetFPKM), "_TPM", sep="")
#res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("Fe_vs_Control_normalized_count_log2FoldChange", "Pvalue")])
res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("log2FoldChange", "pvalue")])
n <- length(colnames(res3))
colnames(res3)[(n-1):n] <- c("D21H6piMinus_vs_D21H6piPlus_Log2FoldChange", "Pvalue")
fdrs <- p.adjust(res3$Pvalue , method ="fdr")
res3 <- data.frame(res3, FDR=fdrs)
res3 <- res3[order(res3$Pvalue, decreasing=FALSE), ]
diff_res3 <- res3[res3$FDR <= 0.05 & abs(res3$D21H6piMinus_vs_D21H6piPlus_normalized_count_log2FoldChange) >= 1,]
write.table(data.frame(Gene_locus_id=rownames(diff_res3), diff_res3), "D21H6piMinus_vs_D21H6piPlus_R_gene_exp_significant.xls", sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(Gene_locus_id=rownames(res3), res3), "D21H6piMinus_vs_D21H6piPlus_R_gene_exp.xls", sep="\t", row.names=FALSE, quote=FALSE)

##### D21H24piMinus_vs_D21H24piPlus
d_one <- d[,cons1 %in% c("D21H24piPlus_R", "D21H24piMinus_R")]
cons_one <- c(rep("D21H24piPlus_R",3), rep("D21H24piMinus_R", 3))
design_one <- data.frame(row.names = colnames(d_one), condition=c(cons_one))
rsum_one1 <- apply(d_one[,1:3], 1, sum)
rsum_one2 <- apply(d_one[,4:6], 1, sum)
d_one <- d_one[rsum_one1 >= 30 | rsum_one2 >= 30,]
dds_one <- DESeqDataSetFromMatrix(countData = d_one, colData = design_one, design=~condition)
dds_one <- DESeq(dds_one, fitType="local")
res <- results(dds_one, contrast=c("condition","D21H24piMinus_R","D21H24piPlus_R"), cooksCutoff=FALSE)
res2 <- as.data.frame(res)
gids <- rownames(res2)
length(gids)
#alldata <- read.table(normcoutfile, sep="\t", header=TRUE, row.names=1)
#alldata <- alldata[gids,c(g1s,g2s)]
subsetraw <- rawcount[gids,design$condition %in% c("D21H24piPlus_R", "D21H24piMinus_R")]
colnames(subsetraw) <- paste(colnames(subsetraw), "_raw_count", sep="")
subsetnorm <- normalizedcount[gids,design$condition %in% c("D21H24piPlus_R", "D21H24piMinus_R")]
colnames(subsetnorm) <- paste(colnames(subsetnorm), "_normalized_count", sep="")
subsetFPKM <- fpkm[rownames(subsetnorm),design$condition %in% c("D21H24piPlus_R", "D21H24piMinus_R")]
colnames(subsetFPKM) <- paste(colnames(subsetFPKM), "_TPM", sep="")
#res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("Fe_vs_Control_normalized_count_log2FoldChange", "Pvalue")])
res3 <- data.frame(subsetraw, subsetnorm, subsetFPKM, res2[,c("log2FoldChange", "pvalue")])
n <- length(colnames(res3))
colnames(res3)[(n-1):n] <- c("D21H24piMinus_vs_D21H24piPlus_Log2FoldChange", "Pvalue")
fdrs <- p.adjust(res3$Pvalue , method ="fdr")
res3 <- data.frame(res3, FDR=fdrs)
res3 <- res3[order(res3$Pvalue, decreasing=FALSE), ]
diff_res3 <- res3[res3$FDR <= 0.05 & abs(res3$D21H24piMinus_vs_D21H24piPlus_normalized_count_log2FoldChange) >= 1,]
write.table(data.frame(Gene_locus_id=rownames(diff_res3), diff_res3), "D21H24piMinus_vs_D21H24piPlus_R_gene_exp_significant.xls", sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(Gene_locus_id=rownames(res3), res3), "D21H24piMinus_vs_D21H24piPlus_R_gene_exp.xls", sep="\t", row.names=FALSE, quote=FALSE)

