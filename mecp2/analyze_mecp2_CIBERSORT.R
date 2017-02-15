library(reshape2)
library(ggplot2)
library(nnet)

options(stringsAsFactors = F)
#######################################
# load the cibersort and phenotype data
load('rdas/pheno.rda')

dat = read.delim('tables/celltype_mecp2_CIBERSORT.txt',row.names = 1)[,1:7]
id_vars = c('Genotype')

###############################
# reformat into long data table
datLong = melt(cbind(dat,Genotype = pd[match(rownames(dat),pd$SampleID),'Genotype']),
               id.vars = id_vars,variable.name = 'Celltype',
               value.name = "Fraction")
datLong$Type = paste0(datLong$Celltype,'.',datLong$Region)

##############################################
# fit multinomial model for shift in p1 brains
tmp = summary(lm(formula = Fraction~Genotype,data = datLong))
coefs = tmp$coefficients[2,'Estimate']
pvals = tmp$coefficients[2,'Pr(>|t|)']

pvals[pvals<0.05]
coefs[pvals<0.05]