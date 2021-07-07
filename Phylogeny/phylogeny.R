require(VariantAnnotation);
require(reshape);
library(reshape2)
library(ggplot2)
library(ggdendro)


vcf=readVcf("strictSomatic.vcf")
AO = as.data.frame(geno(vcf)$AO);
RO = geno(vcf)$RO;
DP = geno(vcf)$DP;

AF=matrix(rep(0, NROW(AO)*NCOL(AO)), ncol=NCOL(AO), byrow=T);
for(i in 1:NROW(AO)) {
  for(j in 1:NCOL(AO)) {
    AF[i, j] = AO[[i,j]] / (AO[[i,j]] + RO[[i, j]]);
  }
}
colnames(AF)=colnames(AO)
rownames(AF)=rownames(AO)

################### Tumor Purity ##################
tp = c(0.51,0.865,0.85,0.85,0.6,0.92,
       0.91,0.95,0.94,0.79,0.88,0.9,0.9,0.51,0.81,
       0.95,0.93,0.91,0.8,0.85,0.73,0.94,0.64,0.66,
       0.83,0.84,0.89,0.67)

correctAF=AF
for(j in 1:length(tp)){
  correctAF[,j]=AF[,j]/tp[j]
}

#################### ggplot of heatmap ###############
correctAF_bi = 0+(correctAF>0.1)
mydist=function(d){
  mat=matrix(,nrow=ncol(d), ncol=ncol(d))
  rownames(mat)=colnames(d)
  colnames(mat)=colnames(d)
  for(i in 1:ncol(d)){
    for(j in 1:ncol(d)){
      a=as.vector(d[,i])
      b=as.vector(d[,j])
      mat[i,j]=length(a)-sum(a-b==0)
    }
  }
  return(as.dist(mat))
}
myhclust=function(d){return(hclust(d, 'ave'))}
sample_dendro=as.dendrogram(myhclust(d=mydist(correctAF_bi)))
sampledendroplot=ggdendrogram(data=sample_dendro)
print(sampledendroplot)
sample_order = order.dendrogram(sample_dendro)
###### cluster variants #####
heatmap_data4=correctAF_bi[apply(correctAF[,c(1:26)]>0.2,1,sum)>0,]
variant_dendro=as.dendrogram(myhclust(d=mydist(t(heatmap_data4))))
variantdendroplot=ggdendrogram(data=variant_dendro,rotate=TRUE)
print(variantdendroplot)
variant_order=order.dendrogram(variant_dendro)

heatmap_data3 = correctAF[apply(correctAF[,c(1:26)]>0.2,1,sum)>0,sample_order]
heatmap_data5 = heatmap_data3[variant_order,]
df_heatmap=melt(heatmap_data5,id.vars="variants")
ggplot(df_heatmap, aes(Var2, Var1)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  ylab("variants") +
  xlab("samples") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Allele Frequency")

###########  UPGMA Clustering Method based on hamming distance  ###########
plot_tree=function(data){
  mat=matrix(,nrow=ncol(data), ncol=ncol(data))
  rownames(mat)=colnames(data)
  colnames(mat)=colnames(data)
  for(i in 1:ncol(data)){
    for(j in 1:ncol(data)){
      a=as.vector(data[,i])
      b=as.vector(data[,j])
      mat[i,j]=length(a)-sum(a-b==0)
    }
  }
  plot(hclust(as.dist(mat), "ave"))
}
## since we know BrP, BrM and Ln1 are samples mixtured with different evolution lineage. We exclude them from this analysis.
plot_tree(correctAF_bi[,c(2:26)])
