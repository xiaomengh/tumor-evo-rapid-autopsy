######### plot CN from new_data ###########
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(MASS)
library(pvclust)
new_data_raw = as.data.frame(read.table("RNAseq_normalized101gene.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1))
allsamplenames = colnames(new_data_raw)
goodsamples= c(16,35,19,32,
               2,20,
               1,5,10,17,18,36,21,15,34,23,
               24,39,
               25,3)

badsamples = c(13,40,12,22,38,7,14,11)

normalsamples = c(4,6,8,9,26,27,28,29,30,31,33,37,41,42)

inMatrix = new_data_raw[,c(goodsamples,badsamples,normalsamples,43)]
pdf("Rnaseq_cna_allsamples_groups_clustered.pdf")
Heatmap(inMatrix[,-ncol(inMatrix)], cluster_rows = F, cluster_columns = F, show_row_names = FALSE, split=as.numeric(inMatrix[,ncol(inMatrix)]), 
        col = colorRamp2(c(-0.5,0,0.5), c("blue", "wheat1", "red")), column_names_gp = gpar(fontsize = 9), row_title_gp = gpar(fontsize = 9))
dev.off()
### plot PCA #####
pheno = read.table("phenoTable.txt",header=T)
for (i in 1:nrow(pheno)){
  if (pheno$sample_type[i]=='n'){
    pheno$group_color[i]= "gray"
  } else {
    if (pheno$tumor_genomic_group[i] == 'G1') {
      pheno$group_color[i] = "green"
    } else if (pheno$tumor_genomic_group[i] == 'G2') {
      pheno$group_color[i] = "purple"
    } else if (pheno$tumor_genomic_group[i] == 'G3') {
      pheno$group_color[i] = "orange"
    } else if (pheno$tumor_genomic_group[i] == 'G4') {
      pheno$group_color[i] = "pink"
    } else if (pheno$tumor_genomic_group[i] == 'Ln7') {
      pheno$group_color[i] = "blue"
    } else if (pheno$tumor_genomic_group[i] == 'Ln9') {
      pheno$group_color[i] = "dark blue"
    } else if (pheno$tumor_genomic_group[i] == 'Ln1') {
      pheno$group_color[i] = "black"
    } else if (pheno$tumor_genomic_group[i] == 'G5') {
      pheno$group_color[i] = "brown"
    }
  }
}

plotPCAColoredByGroup=function(data, selector){
  svd2 = svd(data[,selector])
  plot(svd2$d)
  plot(svd2$v[,1],svd2$v[,2], xlab="1st PC", ylab="2nd PC", col="white")
  text(svd2$v[,1],svd2$v[,2], labels=colnames(data[,selector]), col=pheno[selector,]$group_color,cex=1)
}
all_samples = c(goodsamples,badsamples,normalsamples)
plotPCAColoredByGroup(new_data_raw,all_samples)

### Hierarchical clustering plot #####
plotHclust = function(cluster_data) {
  d = dist(as.matrix(cluster_data), method="euclidean")
  fit = hclust(d, method="ward.D")
  plot(fit)
}
plotHclust(t(new_data_raw[,goodsamples]))



