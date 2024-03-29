# make phenotype data for sweatt mice data

library(jaffelab)
library(parallel)
dir = '/dcl01/lieber/ajaffe/Brady/sweatt'

##########################################
# make bam into fastq to align with HISAT2
bamFile = list.files(path = "/dcl01/lieber/ajaffe/Brady/sweatt/BAM",
                     pattern = "].bam$",full.names = T)
FileID = ss(ss(ss(bamFile,"\\[",2),"\\.",1),"\\]",1)

#resortedBam = paste0(dir,'/BAM/',FileID,'.bam')
#thecall = paste('samtools sort -n',bamFile, paste0(dir,'/BAM/',FileID))
#mclapply(thecall,system,mc.cores = 16)

R1 = paste0(dir,'/FASTQ/',FileID,'.fastq')
R2 = paste0(dir,'/FASTQ/',FileID,'_R2.fastq')
thecall = paste('qsub -V -l mf=10G,h_vmem=15G,h_stack=256M -cwd -b y bedtools bamtofastq -i',bamFile,'-fq',R1)
# mclapply(thecall,system, mc.cores= 18)
# mclapply(paste('gzip',R1),system, mc.cores= 18)

#####################
# make phenotype data
pd = data.frame(SampleID = FileID, FileID = FileID)
pd$Genotype = factor(ifelse(grepl("PTHS",pd$FileID),"HT","WT"))
pd$Hemi = factor(ifelse(grepl("L",pd$FileID),"Left","Right"))
save(pd,file = './rdas/pheno.rda')
write.table(FileID,file = paste0(dir,'/SAMPLE_IDs.txt'),col.names = F,row.names = F,quote = F)
