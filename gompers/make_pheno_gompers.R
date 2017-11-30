# make the phenotype table for gompers het mouse
library(jaffelab)

dir = '/dcl01/lieber/ajaffe/Brady/gompers'
pd = read.delim('tables/gompers_SraRunTable.txt')
rownames(pd) = pd$Run_s
pd = pd[pd$LibrarySource_s =='TRANSCRIPTOMIC',apply(pd,2,function(x) !all(is.na(x)))]
pd$Genotype = factor(pd$genotype_s, levels = c('WT','Het'))
pd$Age = factor(pd$Stage_s,levels = c('e12.5','e14.5','e17.5','P0','P56'))
pd$FASTQ.gz = paste0(dir,'/FASTQ/',pd$Run_s,'.fastq.gz')

# get the paired end fastq files
whichPE = which(pd$LibraryLayout_s=='PAIRED')
pd$FASTQ.gz[whichPE] = paste0(dir,'/FASTQ/',pd$Run_s[whichPE],'_1.fastq.gz')
pd$FASTQ.gz2 = NA
pd$FASTQ.gz2[whichPE] = paste0(dir,'/FASTQ/',pd$Run_s[whichPE],'_2.fastq.gz')

# see how many samples per condition
table(file.exists(c(pd$FASTQ.gz,pd$FASTQ.gz2)))
with(pd[!duplicated(pd$Run_s),],table(Age,Genotype))

#merged fastq files
save(pd,file = 'rdas/pheno_gompers.rda')

#make manifest file for single end reads
whichInd = which(pd$LibraryLayout_s=='SINGLE')
manifest = paste(pd$FASTQ.gz[whichInd], 0, pd$Run_s[whichInd],sep = '\t')
cat(manifest, file = paste0(dir,'/samples.manifest'),sep = '\n')

#make manifest file for paired end reads
whichInd = which(pd$LibraryLayout_s=='PAIRED')
manifest = paste(pd$FASTQ.gz[whichInd], 0,pd$FASTQ.gz2[whichInd],0, pd$Run_s[whichInd],sep = '\t')
cat(manifest, file = paste0('/dcl01/lieber/ajaffe/Brady/gompers_PE/samples.manifest'),sep = '\n')



#Andrew's code 2017NOV30
library(jaffelab)

## reads
pd = read.delim("gompers_chd8_SraRunTable.txt", as.is=TRUE)
fqPath = "/dcl01/lieber/ajaffe/Brady/chd8_het/FASTQ/"

## only P0 and P56
pd = pd[pd$Stage_s %in% c("P0", "P56"),]

dat = data.frame(read = paste0(fqPath, pd$Run_s, ".fastq.gz"),
                 MD5 = 0, SampleID = pd$Run_s,stringsAsFactors=FALSE)
all(file.exists(dat$read))

write.table(dat, file = "samples.manifest",
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
prep_samples.R (END) 




##########################
# run the RNA seq pipeline
# /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/rnaseq-run-all.sh --experiment gompers --prefix SE --reference mm10 --stranded "reverse" --ercc "FALSE" --cores 8 --large "TRUE" --fullcov "FALSE"
# /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/rnaseq-run-all.sh --experiment gompers --prefix PE --reference mm10 --stranded "reverse" --ercc "FALSE" --cores 8 --large "TRUE" --fullcov "FALSE"