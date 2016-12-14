# analyze human down syndrome brain RNAseq
source('/users/ajaffe/Lieber/Projects/RnaQualPaper/OLD/qsva_functions.R')

########################
## load degradation data
load("/dcl01/lieber/ajaffe/Brady/down/down_human_brains_samples.rda")
ind = which(!pd$Region %in% c('FC','HIP','IPC', 'MFC','S1C','STC'))
pd = pd[ind,]
p = exonCounts[,ind]

tIndexes=split(seq(along=pd$Region), pd$Region)
sapply(tIndexes, function(ii) {
  xx = tapply(pd$RIN[ii],pd$Case[ii], mean)
  xx[-1] - xx[1]
})

### make qSVA 
load("/users/ajaffe/Lieber/Projects/RnaQualPaper/OLD/ABRF/ABRF_RNase_RIN_genes.rda")
getPcaVars = function(pca)  signif(((pca$sdev)^2)/(sum((pca$sdev)^2)),3)*100
digIndex=which(map$Gene %in% theGenes$Symbol)

coefMatQsva = lapply(tIndexes, function(ii) {
  mod = model.matrix(~Case+AgeGroup,data = pd[ii,])
  colnames(mod)[2] = levels(pd$Case)[-1]
  pp = p[,ii]
  pca = prcomp(t(pp[digIndex,]))
  nsv = num.sv(pp[digIndex,], mod)
  qual = pca$x[,1:nsv]
  lmFit(pp, cbind(mod,qual))$coef[,2]
})
coefMatQsva = do.call("cbind", coefMatQsva)
colnames(coefMatQsva) = paste0(names(tIndexes),":", colnames(coefMatQsva))


### overall RIN effect
coefMat = lapply(tIndexes, function(ii) {
  mod = model.matrix(~pd$Case[ii])
  colnames(mod)[2] = levels(pd$Case)[-1]
  lmFit(p[,ii], mod)$coef[,2]
})
coefMat = do.call("cbind", coefMat)
colnames(coefMat) = paste0(names(tIndexes),":", colnames(coefMat))

tMat = lapply(tIndexes, function(ii) {
  mod = model.matrix(~pd$Case[ii])
  colnames(mod)[2] = levels(pd$Case)[-1]
  ebayes(lmFit(p[,ii], mod))$t[,2]
})
tMat = do.call("cbind", tMat)
colnames(tMat) = paste0(names(tIndexes),":", colnames(tMat))

coefMatRin = sapply(tIndexes, function(ii) {
  modRin = model.matrix(~pd$RIN[ii])
  lmFit(p[,ii], modRin)$coef[,2]
})
tMatRin = sapply(tIndexes, function(ii) {
  modRin = model.matrix(~pd$RIN[ii])
  ebayes(lmFit(p[,ii], modRin))$t[,2]
})
pMatRin = sapply(tIndexes, function(ii) {
  modRin = model.matrix(~pd$RIN[ii])
  ebayes(lmFit(p[,ii], modRin))$p[,2]
})

modRin = model.matrix(~pd$RIN + pd$Region)
fitRin = lmFit(p,modRin)
ebRin = ebayes(fitRin)


cor(coefMatQsva, fitRin$coef[,2])
cor(coefMatQsva, coefMatRin)
tMatQsva = lapply(tIndexes, function(ii) {
  mod = model.matrix(~pd$Case[ii])
  colnames(mod)[2] = levels(pd$Case)[-1]
  pp = p[,ii]
  pca = prcomp(t(pp[digIndex,]))
  nsv = num.sv(pp[digIndex,], mod)
  qual = pca$x[,1:nsv]
  ebayes(lmFit(pp, cbind(mod,qual)))$t[,2]
})
tMatQsva = do.call("cbind", tMatQsva)
colnames(tMatQsva) = paste0(names(tIndexes),":", colnames(tMatQsva))

cor(tMatQsva, tMatRin)
cor(tMatQsva, ebRin$t[,2])
