# make phenotype data for mef2c mice RNAseq SRP090281
library(jaffelab)

dir = "/dcl01/lieber/ajaffe/Brady/mef2c/"
pd = read.delim("./tables/SraRunTable.txt",stringsAsFactors = FALSE)
pd$Genotype = factor(pd$genotype_variation_s, levels=c("control (Mef2cfl/fl)", "Mef2c cKO (Mef2cfl/fl; Emx1-Cre)"))
pd$SampleID = pd$Run_s
save(pd,file = './rdas/pheno.rda')
write.table(pd$SampleID,file = paste0(dir,'/SAMPLE_IDs.txt'),col.names = F,row.names = F,quote = F)
