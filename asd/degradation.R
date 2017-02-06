# make degradation matrix
library(sva)
dir = '/dcl01/lieber/ajaffe/Brady/asd/'
load('rdas/pheno.rda')
pd = cbind(pd,read.csv('/dcl01/lieber/ajaffe/Brady/asd/annotated_pd.csv'))
pd$degCov = paste0(dir,'degradation/',pd$SAMPLE_ID,'_degradeStats_Ribozero.txt')
table(file.exists(pd$degCov))

degMatList = parallel::mclapply(pd$degCov, function(x) {
  cat(".")
  read.delim(pipe(paste("cut -f10", x)), as.is=TRUE,skip=1)[,1]
}, mc.cores=8)
degMat = do.call("cbind", degMatList)/50 # sum reads/read length
colnames(degMat) = pd$SAMPLE_ID
qSVs = qsva(degMat) # 5 PCs

save(degMat,qSVs,file = 'rdas/qSVAs-geschwind_asd.rdas' )