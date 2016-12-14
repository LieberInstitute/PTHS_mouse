# make phenotype data for human huntington's disease GSE64810
library(jaffelab)

dir = "/dcl01/lieber/ajaffe/Brady/huntington/"
pd = read.delim("./tables/SraRunTable.txt",stringsAsFactors = FALSE)
pd$SampleID = pd$Run_s
pd$Case = factor(ifelse(grepl('normal',pd$diagnosis_s),'CTL','HNT'))
pd$RIN = as.numeric(pd$rin_s)
pd$PMI = as.numeric(pd$pmi_s)
pd$Age = as.numeric(pd$age_of_death_s)

save(pd,file = './rdas/pheno.rda')
write.table(pd$SampleID,file = paste0(dir,'/SAMPLE_IDs.txt'),col.names = F,row.names = F,quote = F)


#############
# re-write 
r1 = paste0(dir,'FASTQ/',pd$SampleID,'_1.fastq.gz')
r1_new = paste0(dir,'FASTQ/',pd$SampleID,'_R1_001.fastq.gz')

r2 = paste0(dir,'FASTQ/',pd$SampleID,'_2.fastq.gz')
r2_new = paste0(dir,'FASTQ/',pd$SampleID,'_R2_001.fastq.gz')

# file.exists(r1) #TRUE
#parallel::mclapply(paste('mv',r1,r1_new),system,mc.cores = parallel::detectCores())
#parallel::mclapply(paste('mv',r2,r2_new),system,mc.cores = parallel::detectCores())