
###########  UPGMA Clustering Method  ###########
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

