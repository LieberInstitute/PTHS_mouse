# make the phenotype table for chd8 het mouse
library(readxl)

dir = '/dcl01/lieber/ajaffe/Brady/chd8'
pd = read_excel('tables/DRA003116_pheno.xlsx')
pd = pd[pd$Expr_Type =='TRANSCRIPTOMIC',apply(pd,2,function(x) !all(is.na(x)))]
pd$Genotype = factor(pd$Genotype, levels = c('WildType','CHD8hetero'))
pd$Age = factor(pd$Age,levels = c('E10.5','E12.5','E14.5','E16.5','E18.5','Adult'))
pd$FASTQ.bz2 = paste0(dir,'/FASTQ/',pd$Experiment_s,'/',pd$Run_s,'.fastq.bz2')
pd$FASTQ.gz = paste0(dir,'/FASTQ/',pd$Experiment_s,'/',pd$Run_s,'.fastq.gz')

table(file.exists(pd$FASTQ.bz2))

with(pd[!duplicated(pd$Sample_s),],table(Age,Genotype))
# N of 1 for each age/genotype

# change bz2 to gz
thecall = with(pd, paste('bunzip2 -c <',FASTQ.bz2,'| gzip -c > ',FASTQ.gz))
#parallel::mclapply(thecall,system, mc.cores = 36)
table(file.exists(pd$FASTQ.gz))

#merged fastq files

pd = pd[!duplicated(pd$Sample_s),c("TITLE","Sample_s",
                                   "Age", "Genotype", "Expr_Type")]
pd$FASTQ = paste0(dir,'/merged_fastq/',pd$Sample_s,'.fastq.gz')
table(file.exists(pd$FASTQ))

save(pd,file = 'rdas/pheno_chd8.rda')

#make manifest file

manifest = with(pd,paste(FASTQ, 0, Sample_s,sep = '\t'))

cat(manifest, file = paste0(dir,'/samples.manifest'),sep = '\n')

##########################
# run the RNA seq pipeline
# /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/rnaseq-run-all.sh --experiment chd8 --prefix feb7 --reference mm10 --stranded "TRUE" --ercc "FALSE" --cores 8 --large "TRUE" --fullcov "FALSE"