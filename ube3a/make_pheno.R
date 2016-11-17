library(jaffelab)
library(readxl)
library(parallel)

dir = '/dcl01/lieber/ajaffe/Brady/ube3a'
pd = read_excel('../philpot/tables/Thaxton_Philpot_BGI_RNaseq_file_info_Maher_modified .xlsx')

###################################
# recode phenotype data
pd$Line = ss(pd$"File rename/Unique identifier","-")
pd$FlowCell = ss(pd$"BGI file name", "_")
pd$Case = factor(pd$'General genotype', levels=c("Control", "Mutant"))
pd$Age = ifelse(as.numeric(ss(pd$Age,"P",2)) > 30, "Adult","P1")
pd$Age = factor(pd$Age, levels = c("P1", "Adult"))
pd$Sex = factor(pd$Sex, levels = c("M","F"))
pd$Extract = factor(pd$"Initial RNA extraction/ homogenization (Date)")
pd$Repurify = as.character(pd$"re-purification with Zymo (Date")
pd$Repurify = pd$Repurify == "2016-03-22"
pd = pd[pd$Line=='Ube3a',]
pd$Genotype = factor(pd$Genotype)
pd$FASTQ = paste0('/dcl01/lieber/ajaffe/Brady/philpot/FASTQ/',pd$"BGI file name")
#file.exists(pd$FASTQ) #TRUE

save(pd,file = './rdas/pheno.rda')
cat(gsub('.fq.gz','',pd$"BGI file name"),file = paste0(dir,'/SAMPLE_IDs.txt'),sep = '\n')