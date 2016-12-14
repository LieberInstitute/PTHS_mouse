# make phenotype data for parkinson mice, MPTP lession RNAseq PRJEB6957
library(jaffelab)

dir = "/dcl01/lieber/ajaffe/Brady/parkinson/"
pd = read.delim("./tables/SraRunTable.txt",stringsAsFactors = FALSE)
pd$Condition = ifelse(grepl('Unlesioned',pd$Sample_Name_s),'Unlesioned','Lesioned')
pd$Condition = factor(pd$Condition,c('Unlesioned','Lesioned'))
pd$SampleID = pd$Run_s
save(pd,file = './rdas/pheno.rda')
#write.table(pd$SampleID,file = paste0(dir,'/SAMPLE_IDs.txt'),col.names = F,row.names = F,quote = F)
