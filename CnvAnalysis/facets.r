library(facets)
set.seed(1234)
sample ="tumor"
rcmat = readSnpMatrix(paste(sample,'.snppileup.csv',sep=''))
xx = preProcSample(rcmat, ndepth=15)
oo = procSample(xx, min.nhet=5)
fit = emcncf(oo)

print("ploidy")
print(fit$ploidy)

print("purity")
print(fit$purity)

rdatafile = paste(sample,'.RData',sep='')
save(xx, oo, fit, file=rdatafile)

fitfile = paste(sample,'.cncf.txt',sep='')
write.table(fit$cncf, fitfile, sep="\t", row.names=F, col.names=T, quote=F)
pdf(paste(sample,".facets.pdf",sep=''))
plotSample(x=oo, emfit=fit)
dev.off()
