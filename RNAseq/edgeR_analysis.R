library(RnaSeqGeneEdgeRQL)
library(edgeR)
library(statmod)
library(org.Hs.eg.db)

pheno = read.table("phenoTable.txt",header=T)
data1 = read.table("genewiseCounts_entrezid.txt",header=T)
### all pure sample comparison ###
sample_sel=which(pheno$include_in_pure_tumor_analysis == "y")
group = pheno[sample_sel,]$genomic_group
group <-factor(group)
# process data #
y <- DGEList(as.matrix(data1[,sample_sel]), group=group,genes=data1[,1,drop=FALSE])
y$samples
y$genes$Symbol <- mapIds(org.Hs.eg.db, rownames(y),keytype="ENTREZID", column="SYMBOL")
head(y$genes$Symbol)
design=model.matrix(~0+group)
colnames(design) <- levels(group)
keep <- filterByExpr(y, design)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
AveLogCPM <- aveLogCPM(y)
hist(AveLogCPM)
y <- calcNormFactors(y)
y$samples
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)
# end #
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)
plotQLDisp(fit)
logCPM <- cpm(y, prior.count=2, log=TRUE)
rownames(logCPM) <- y$genes$Symbol
pdf("coolmap_allgenes_20samples.pdf")
coolmap(logCPM, margins=c(7,7), lhei=c(2,10), lwid=c(2,3),cexRow=0.5)
dev.off()

###Pathway analysis###
sample_sel=which(pheno$include_in_pure_tumor_analysis == "y" & 
                   pheno$genomic_group=="G1" | pheno$genomic_group=="G3")
group = pheno[sample_sel,]$genomic_group
group <-factor(group)
# process data #
y <- DGEList(as.matrix(data1[,sample_sel]), group=group,genes=data1[,1,drop=FALSE])
y$samples
y$genes$Symbol <- mapIds(org.Hs.eg.db, rownames(y),keytype="ENTREZID", column="SYMBOL")
head(y$genes$Symbol)
design=model.matrix(~0+group)
colnames(design) <- levels(group)
keep <- filterByExpr(y, design)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
AveLogCPM <- aveLogCPM(y)
hist(AveLogCPM)
y <- calcNormFactors(y)
y$samples
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)

Hs.c2 <- readRDS(url("http://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.c2.all.v7.1.entrez.rds"))
idx <- ids2indices(Hs.c2,id=row.names(y))
g1vsg3 <- makeContrasts(
  G1sG3=G1-G3,
  levels=design)

cam <- camera(y, idx, design, contrast=g1vsg3, inter.gene.cor=0.01)
options(digits=2)
head(cam,5)
write.table(head(cam,50),"MsigDB_g1vsg3.txt",quote=F,sep="\t")
res <- glmQLFTest(fit, contrast=g1vsg3)
barcodeplot(res$table$logFC,index=idx[["NIKOLSKY_BREAST_CANCER_17Q21_Q25_AMPLICON"]],index2=idx[["NIKOLSKY_BREAST_CANCER_8Q23_Q24_AMPLICON"]],
            labels=c("G1","G3"),
            main="G1 vs. G3",
            alpha=1)

##################### Lung sample cluster ###########
sample_sel=which(pheno$include_in_pure_tumor_analysis == "y" & pheno$tissue == "Ln")
group = pheno[sample_sel,]$genomic_group
group <-factor(group)
y <- DGEList(as.matrix(data1[,sample_sel]), group=group,genes=data1[,1,drop=FALSE])
y$samples
y$genes$Symbol <- mapIds(org.Hs.eg.db, rownames(y),keytype="ENTREZID", column="SYMBOL")
head(y$genes$Symbol)
design=model.matrix(~0+group)
colnames(design) <- levels(group)
keep <- filterByExpr(y, design)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
AveLogCPM <- aveLogCPM(y)
hist(AveLogCPM)
y <- calcNormFactors(y)
y$samples
plotMD(y, column=5)
abline(h=0, col="red", lty=2, lwd=2)
y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)
plotQLDisp(fit)
logCPM <- cpm(y, prior.count=2, log=TRUE)
rownames(logCPM) <- y$genes$Symbol
pdf("coolmap_allgenes_Lnsamples.pdf")
coolmap(logCPM, margins=c(7,7), lhei=c(2,10), lwid=c(2,3),cexRow=0.5)
dev.off()
##################### liver sample cluster ###########
sample_sel=which(pheno$include_in_pure_tumor_analysis == "y" & pheno$tissue == "Lv")
group = pheno[sample_sel,]$genomic_group
group <-factor(group)
y <- DGEList(as.matrix(data1[,sample_sel]), group=group,genes=data1[,1,drop=FALSE])
y$samples
y$genes$Symbol <- mapIds(org.Hs.eg.db, rownames(y),keytype="ENTREZID", column="SYMBOL")
head(y$genes$Symbol)
design=model.matrix(~0+group)
colnames(design) <- levels(group)
keep <- filterByExpr(y, design)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
AveLogCPM <- aveLogCPM(y)
hist(AveLogCPM)
y <- calcNormFactors(y)
y$samples
plotMD(y, column=5)
abline(h=0, col="red", lty=2, lwd=2)
y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)
plotQLDisp(fit)
logCPM <- cpm(y, prior.count=2, log=TRUE)
rownames(logCPM) <- y$genes$Symbol
pdf("coolmap_allgenes_Lvsamples.pdf")
coolmap(logCPM, margins=c(7,7), lhei=c(2,10), lwid=c(2,3),cexRow=0.5)
dev.off()


