# phenotype data for Psych ENCODE ASD
library(jaffelab)

load('/dcl01/lieber/ajaffe/psychENCODE_Data/UCLA_R01MH094714/rdas/annotated_phenotype_asd.rda')
dat1 = '/dcl01/lieber/ajaffe/psychENCODE_Data/UCLA_R01MH094714/FASTQ/'
dat2 = '/dcl01/lieber/ajaffe/Brady/asd/FASTQ/'

##########################################
# remove 15q duplication syndrome patients
pd = pd[pd$Detailed.Diagnosis!="Chromosome 15q Duplication Syndrome",]

####################
# rename fastq files
pd$SampleID  = ss(pd$Fastq.file.names,'.R1/2.fastq.gz')
old = c(paste0(dat1,pd$SampleID,'.R1.fastq.gz'), paste0(dat1,pd$SampleID,'.R2.fastq.gz'))
new = c(paste0(dat2,pd$SampleID,'_R1_001.fastq.gz'), paste0(dat2,pd$SampleID,'_R2_001.fastq.gz'))
thecall = paste('cp',old,new)
#parallel::mclapply(thecall, system, mc.cores = parallel::detectCores())

cat(paste0(dat2,pd$SampleID),file = '/dcl01/lieber/ajaffe/Brady/asd/SAMPLE_IDs.txt',sep = '\n')
save(pd,file = 'rdas/pheno.rda')