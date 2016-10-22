library(jaffelab)
library(readxl)
library(parallel)

dir = '/dcl01/lieber/ajaffe/Brady/philpot'
pd = read_excel('./tables/Thaxton_Philpot_BGI_RNaseq_file_info_Maher_modified .xlsx')

###################################
# put all reads in FASTQ folder
#fns = list.files(paste0(dir,'/Reads/Clean'),recursive = T,pattern = '.gz',full.names = T)
#fns = fns[unlist(sapply(pd$"BGI file name",grep,fns))]
thepath = paste0(dir,'/FASTQ/')
pd$FASTQ = paste0(thepath,pd$"BGI file name")
#thecall = paste('cp',fns,pd$FASTQ)
# mclapply(thecall,system,mc.cores = 18)
#write.table(gsub('.fq.gz','',pd$"BGI file name"),file = paste0(dir,'/SAMPLE_IDs.txt'),
#            col.names = F,row.names = F,quote = F)

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
pd = pd[pd$Line!='Ube3a',]
pd$Genotype = ifelse(pd$'General genotype'=='Control','WT',pd$Line)
pd$Genotype = factor(pd$Genotype,levels = c('WT','Act','Nest','Del','R579W'))

save(pd,file = './rdas/pheno.rda')
#write.table(gsub('.fq.gz','',pd$"BGI file name"),file = paste0(dir,'/SAMPLE_IDs.txt'),
 #           col.names = F,row.names = F,quote = F)