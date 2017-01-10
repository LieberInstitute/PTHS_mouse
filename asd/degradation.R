# make degradation matrix
library(sva)

load('rdas/pheno.rda')
pd = cbind(pd,read.delim())