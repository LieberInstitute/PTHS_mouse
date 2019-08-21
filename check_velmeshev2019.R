##
library(lmerTest)

meta = read.delim("asd/tables/meta_velmeshev2019.tsv",as.is=TRUE)

## table
tt = table(meta$sample, meta$cluster)
tt = prop.table(tt,1)
class(tt) = "numeric"

## reduce pheno
pd = meta[match(rownames(tt), meta$sample),3:12]
pd$diagnosis = factor(pd$diagnosis, levels = c("Control", "ASD"))
colnames(pd)[9:10] = c("PMI", "RIN")

## analysis
pd = cbind(pd, as.matrix(tt))

## univariate
lmObjs = lapply(pd[,11:ncol(pd)], function(x) summary(lm(x ~ pd$diagnosis + pd$region))$coef)
dxEffect = t(sapply(lmObjs, function(x) x[2,]))
signif(dxEffect,3)

dxEffect = dxEffect[order(dxEffect[,4]),]
colnames(dxEffect) = c("PropDiff", "SE", "T-stat","P-value")
write.csv(dxEffect, "asd/tables/table_s7_velmeshev2019.csv")