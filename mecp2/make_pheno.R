# make phenotype data for mecp2 mice RNAseq
library(jaffelab)

dir = "/dcl01/lieber/ajaffe/Brady/mouseRNAseq/Rett"
pd = read.csv("./tables/SraRunInfo.csv",stringsAsFactors = FALSE)
pd$Genotype = factor(pd$Geno, levels=c("wild type", "MeCP2 KO"))
save(pd,file = './rdas/pheno.rda')


write.table(pd$SampleID,file = paste0(dir,'/SAMPLE_IDs.txt'),col.names = F,row.names = F,quote = F)
