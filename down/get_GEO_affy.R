# extract Down syndrome microarray data for DE, GSE59630
library(GEOquery)
library(jaffelab)
options(stringsAsFactors = F)

###########################
# get unfiltered annotation
dat = getGEO("GSE59630")[[1]]
exonCounts = exprs(dat)

############################
# restructure phenotype data
pd = pData(dat)
pd$FileID = as.character(pd$geo_accession)
pd$SampleID = ss(as.character(pd$characteristics_ch1),": ",2)
pd$Region = ss(as.character(pd$characteristics_ch1.1),": ",2)
pd$Case = factor(ss(as.character(pd$characteristics_ch1.2),": ",2),
                 levels = c('CTL','DS'))
pd$Sex = ss(as.character(pd$characteristics_ch1.3),": ",2)
pd$Age = ss(as.character(pd$characteristics_ch1.4),": ",2)
pd$AgeGroup = ifelse(grepl('wg',pd$Age),'Fetal',
                     ifelse(grepl('mo',pd$Age),'Infant',
                            ifelse(as.numeric(ss(pd$Age,'yr'))<10,'Child',
                                   ifelse(as.numeric(ss(pd$Age,'yr'))<20,'Teen','Adult'))))
pd$AgeGroup = factor(pd$AgeGroup,levels = 
                       c('Fetal','Infant','Child','Teen','Adult'))
pd$Stage = as.numeric(ss(as.character(pd$characteristics_ch1.5),": ",2))
pd$PMI = as.numeric(ss(as.character(pd$characteristics_ch1.6),": ",2))
pd$Race = ss(as.character(pd$characteristics_ch1.7),": ",2)
pd$RIN = as.numeric(ss(as.character(pd$characteristics_ch1.8),": ",2))

######################
# clean up annotations
map = fData(dat)
map$Gene = ss(as.character(map$gene_assignment),' // ',2)

#####################
# save the data
save(exonCounts,map,pd,file="/dcl01/lieber/ajaffe/Brady/down/down_human_brains_samples.rda",compress=TRUE)





