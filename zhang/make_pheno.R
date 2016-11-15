# make pheno for Zhang brain cell types RNA seq data
library(jaffelab)

dir = "/dcl01/lieber/ajaffe/Brady/zhang/"
pd = read.delim("./tables/SraRunTable-Zhang.txt",stringsAsFactors = FALSE)
pd$Cell = factor(pd$cell_type_s)
pd$SampleID = pd$Run_s
pd = pd[,apply(pd,2,function(x)!all(x=='<not provided>'))] #rm unnecessary columns

save(pd,file = './rdas/pheno.rda')
write.table(pd$SampleID,file = paste0(dir,'/SAMPLE_IDs.txt'),col.names = F,row.names = F,quote = F)

Dir = "/dcl01/lieber/ajaffe/Brady/zhang/FASTQ"
reads = list.files(Dir, pattern="fastq.gz$", recur=TRUE, full.names=TRUE)

R1_old = paste0(Dir,'/',pd$SampleID,'_1.fastq.gz')
R2_old = paste0(Dir,'/',pd$SampleID,'_2.fastq.gz')

R1 = paste0(Dir,'/',pd$SampleID,'_R1.fastq.gz')
R2 = paste0(Dir,'/',pd$SampleID,'_R2.fastq.gz')

thecall1 = paste('cp',R1_old,R1)
thecall2 = paste('cp',R2_old,R2)
# parallel::mclapply(c(thecall1,thecall2),system,mc.cores = 18)

thecall1 = paste('rm',R1_old)
thecall2 = paste('rm',R2_old)
# parallel::mclapply(c(thecall1,thecall2),system,mc.cores = 18)
